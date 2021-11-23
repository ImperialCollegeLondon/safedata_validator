from collections import Counter
from typing import Union, Optional
import dataclasses
from io import StringIO
from itertools import groupby
import sqlite3
import requests
from enforce_typing import enforce_types

from safedata_validator.logger import LOGGER, FORMATTER, CH, log_and_raise
from safedata_validator.validators import (GetDataFrame, IsLower, IsNotBlank,
                                           IsNotPadded, IsString, HasDuplicates)

"""
This module describes classes and methods used to compile taxonomic data from
datasets and to validate taxonomy against the GBIF backbone database.

The Taxon dataclass is used to store data about a taxon entry in a dataset. It
is initialised with user data and then the taxon Validator classes can be used
to update a Taxon object with the result of GBIF validation. The two validation
classes use either the online API (`RemoteGBIFValidator`) or faster validation
against a local copy of the GBIF backbone (`LocalGBIFValidator`).

The dataset 'Taxa' worksheet provides a set of taxonomic entries and the Taxa
class is used to load and collate the set of taxonomic entries from a dataset.

Note that we explicitly exclude form and variety from the set of GBIF backbone
taxonomic levels because they cannot be matched into the backbone hierarchy 
without extra API calls.
"""

BACKBONE_RANKS = ['kingdom', 'phylum', 'order', 'class', 'family',
                  'genus', 'species', 'subspecies']

# TODO - A lot of complexity could be lost here if the two validation
#        sources had more similar structure. Could have a row return method
#        that is used within a generic id_lookup and search class.
#        Notably there us much less information in the row return from the
#        local SQL db (e.g. merging acceptedKey and parentKey, but also not
#        providing parent hierarcy taxon names.  Could defining some SQL
#        queries to add some of that information into the local database?


class GBIFError(Exception):
    """Exception class for remote GBIF errors

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message="GBIF ID not found"):
        self.message = message
        super().__init__(self.message)


@enforce_types
@dataclasses.dataclass
class Taxon:
    """Holds taxonomic information from a user, which can be populated using GBIF validation.

    There are 4 class properties that can be used to create an instance:
        * name
        * rank
        * gbif_id
        * ignore_gbif: a GBIF id match to reject
    The remaining properties are populated by processing functions not when
    an instance is created.
        * is_backbone: the taxon is at a taxonomic level included in the GBIF backbone
        * is_canon: the taxon is considered canon in GBIF
        * lookup_status: the outcome of the lookup with one of the following values:
          found, no_match, validation_fail, unknown_id, id_mismatch
        * taxon_status: the taxonomic status of the taxon with one of the following values:
          accepted, doubtful, synonym etc. etc.
        * parent_id: a GBIF id for the accepted parent taxon.
        * canon_usage: a Taxon instance holding the canonical usage for the taxon
        * note: a string of any extra information provided by the search
        * hierarchy: a list of 2-tuples of rank and GBIF ID for the taxonomic hierarchy
    """
    
    # Init properties
    name: str
    rank: str
    gbif_id: Optional[Union[int, float]] = None
    ignore_gbif: Optional[Union[int, float]] = None
    is_backbone: bool = dataclasses.field(init=False)
    is_canon: bool = dataclasses.field(init=False)
    canon_usage: 'Taxon' = dataclasses.field(init=False)  # https://stackoverflow.com/questions/33533148
    parent_id: int = dataclasses.field(init=False)
    taxon_status: str = dataclasses.field(init=False)
    lookup_status: str = dataclasses.field(init=False)
    note: str = dataclasses.field(init=False)
    hierarchy: list = dataclasses.field(init=False)

    def __post_init__(self):
        """Sets the defaults for the post-init properties and checks inputs
        """

        if self.gbif_id is not None:
            if isinstance(self.gbif_id, float) and not self.gbif_id.is_integer():
                raise ValueError('GBIF Id is not an integer')
            self.gbif_id = int(self.gbif_id)

        if self.ignore_gbif is not None:
            if isinstance(self.ignore_gbif, float) and not self.ignore_gbif.is_integer():
                raise ValueError('Ignore GBIF is not an integer')
            self.ignore_gbif = int(self.ignore_gbif)

        self.rank = self.rank.lower()
        self.is_backbone = self.rank in BACKBONE_RANKS
        self.is_canon = False
        self.canon_usage = None
        self.parent_id = None
        self.taxon_status = None
        self.lookup_status = 'unvalidated'
        self.note = None
        self.hierarchy = []

    def __repr__(self):

        if not self.is_backbone:
            return f"{self.name} (not of GBIF backbone rank)"
        elif self.found:
            if self.is_canon:
                return f"{self.name}"
            else:
                return f"{self.name} ({self.taxon_status})"
        else:
            return f"{self.name} (not found: {self.lookup_status})"

    @property
    def found(self):
        """Convenience boolean property for processing"""
        return self.is_backbone and self.lookup_status == 'found'


class LocalGBIFValidator:

    def __init__(self, resources):

        conn = sqlite3.connect(resources.gbif_database)
        conn.row_factory = sqlite3.Row
        self.gbif_conn = conn

    def __del__(self):

        self.gbif_conn.close()

    def search(self, taxon: Taxon):

        """
        Looks for a taxon in the GBIF database using name and rank and
        an optional GBIF ID for disambiguation. The ignore_gbif optional
        parameter allows a user to override a particular match.

        Args:
            taxon: A Taxon instance

        Returns:
            A Taxon instance
        """

        if not taxon.is_backbone:
            taxon.lookup_status = 'Cannot validate non-backbone taxa'
            return taxon

        if taxon.gbif_id is not None:

            # get the record associated with the provided ID
            id_taxon = self.id_lookup(taxon.gbif_id)

            # Check that name and rank are congruent with id
            if (id_taxon.name != taxon.name) or (id_taxon.rank != taxon.rank):
                taxon.lookup_status = 'ID does not match name and rank'
                return taxon

            return id_taxon

        else:
            # get the set of records associated with the taxon and rank

            sql = (f"select * from backbone where canonical_name ='{taxon.name}' "
                   f"and rank= '{taxon.rank.upper()}';")

            taxon_rows = self.gbif_conn.execute(sql).fetchall()
            selected_row = None

            if len(taxon_rows) == 0:
                # No matching rows
                taxon.lookup_status = 'No match found'
                return taxon
            elif len(taxon_rows) == 1:
                # one matching row - extract it from the list
                selected_row = taxon_rows[0]
            else:
                # More than one row - try to mimic the preferred hits reported
                # by the GBIF API to select a single hit by looking at the counts
                # of the different statuses.

                # First, get the taxon statuses
                tx_status = [tx['status'].lower() for tx in taxon_rows]
                tx_counts = Counter(tx_status)

                if 'accepted' in tx_counts.keys() and tx_counts['accepted'] == 1:
                    # Single accepted hits are first preference
                    selected_row = taxon_rows[tx_status.index('accepted')]
                elif 'doubtful' in tx_counts.keys() and tx_counts['doubtful'] == 1:
                    # Doubtful hits get next preference - not quite sure about this!
                    selected_row = taxon_rows[tx_status.index('doubtful')]
                else:
                    # Rows now contain only synonyms (of varying kinds) and
                    # misapplied. Both of these types have accepted usage
                    # values, so look for a unique accepted usage, trapping the
                    # edge case of kingdoms, which have no parent_key.
                    tx_acc = {tx['parent_key'] for tx in taxon_rows
                              if tx['parent_key'] is not None}

                    if len(tx_acc) == 1:
                        # A single accepted usage - pick the first row to index
                        selected_row = taxon_rows[0]

            if selected_row is None:
                # No single row has been accepted as the best, so return no
                # match and a note, as the API interface does.
                taxon.lookup_status = f'Multiple equal matches for {taxon.name}'
                return taxon

            # Should now have a single row for the preferred hit, which can be
            # extracted from the database
            return self.id_lookup(selected_row['id'])

    def id_lookup(self, gbif_id: int):
        """Method to return a Taxon directly from a GBIF ID

        Params:
            gbif_id: An integer

        Returns:
            A Taxon object.
        """

        if not isinstance(gbif_id, int):
            raise TypeError()

        if not gbif_id > 0:
            raise ValueError()

        # get the record associated with the provided ID
        sql = f"select * from backbone where id = {gbif_id}"
        taxon_row = self.gbif_conn.execute(sql).fetchone()

        # check there is a result and that it is congruent with any
        # provided taxon or rank information
        if taxon_row is None:
            raise GBIFError()

        # Create and populate taxon
        taxon = Taxon(name=taxon_row['canonical_name'],
                      rank=taxon_row['rank'].lower(),
                      gbif_id=taxon_row['id'])
        taxon.lookup_status = 'found'
        taxon.taxon_status = taxon_row['status'].lower()
        taxon.parent_id = taxon_row['parent_key']

        # Add the taxonomic hierarchy, using a mapping of backbone ranks (except
        # subspecies) to backbone table fields.
        taxon.hierarchy = [(rk, taxon_row[ky])
                           for rk, ky in [(r, r + '_key') for r in BACKBONE_RANKS[:-1]]
                           if ky in taxon_row.keys()]

        # parent key in the local database has the odd property that the parent
        # tax_gbif['parent_key'] does dual duty: points up to parent for canon
        # taxa and 'up' to canon for non-canon taxa, so need to look through both
        # to get the canon and parent populated.
        if taxon.taxon_status in ['accepted', 'doubtful']:
            taxon.is_canon = True
        else:
            taxon.is_canon = False
            taxon.canon_usage = self.id_lookup(taxon.parent_id)
            taxon.parent_id = taxon.canon_usage.parent_id

        return taxon


@enforce_types
class RemoteGBIFValidator:
    """This provides a validate method for a Taxon using the online GBIF
    API. Unlike the LocalGBIFValidator, this doesn't need an __init__ method
    and just contains methods, but duplicates the structure so that
    the two Validators are interchangeable.
    """
    def search(self, taxon: Taxon):

        """
        Validates a taxon against the GBIF web API. It uses the API endpoint
        species/match?name=XXX&rank=YYY&strict=true
        endpoint to



        # safe to assume that the next least nested taxonomic level is the parent.
        # So, populate GBIF ID using species/match and then use species/{id} to populate the
        # return values

        Args:
            taxon (Taxon): A Taxon to be validated

        Returns:
            An updated Taxon object.
        """

        if not taxon.is_backbone:
            raise ValueError('Cannot validate non-backbone taxa')

        if taxon.gbif_id is not None:

            # get the record associated with the provided ID
            id_taxon = self.id_lookup(taxon.gbif_id)

            # Check that name and rank are congruent with id
            if (id_taxon.name != taxon.name) or (id_taxon.rank != taxon.rank):
                taxon.lookup_status = 'ID does not match name and rank'
                return taxon

            return id_taxon
        else:

            # If no ID is provided then use the species/match?name={} endpoint to
            # try and find an ID for the combination
            url = (f'http://api.gbif.org/v1/species/match?name={taxon.name}'
                   f'&rank={taxon.rank}&strict=true')
            taxon_row = requests.get(url)

            # catch errors
            if taxon_row.status_code != 200:
                raise GBIFError('Connection error to remote server')

            response = taxon_row.json()

            # check the response status
            if response['matchType'] == u'NONE':
                # No match found - look for explanatory notes
                if 'note' in response:
                    taxon.lookup_status = response['note']
                    return taxon
                else:
                    taxon.lookup_status = 'No match found'
                    return taxon

            # Return the match
            return self.id_lookup(response['usageKey'])

    def id_lookup(self, gbif_id: int):
        """Method to return a Taxon directly from a GBIF ID

        Params:
            gbif_id: An integer

        Returns:
            A Taxon object.
        """

        if not isinstance(gbif_id, int):
            raise TypeError()

        if not gbif_id > 0:
            raise ValueError()

        # Use the species/{id} endpoint
        taxon_row = requests.get(f"http://api.gbif.org/v1/species/{gbif_id}")

        # unknown ID numbers return a 404 error
        if taxon_row.status_code == 404:
            raise GBIFError()
        elif taxon_row.status_code != 200:
            raise GBIFError('Connection error to remote server')

        # Extract the response
        response = taxon_row.json()

        # Create and populate taxon
        taxon = Taxon(name=response['canonicalName'],
                      rank=response['rank'].lower(),
                      gbif_id=response['nubKey'])

        # First, set the parent key - in the GBIF API, this is always provided
        # and always points to the true parent taxon, unlike parent_key in
        # the simple local database. Note that parentKey is not included in the
        # response for Kingdom level taxa, hence get().
        taxon.parent_id = response.get('parentKey')
        taxon.taxon_status = response['taxonomicStatus'].lower()
        taxon.lookup_status = 'found'

        # Add the taxonomic hierarchy from the accepted usage - these are tuples
        # to be used to extend a set for the taxonomic hierarchy
        taxon.hierarchy = [(rk, response[ky])
                           for rk, ky in [(r, r + 'Key') for r in BACKBONE_RANKS[:-1]]
                           if ky in response]

        # Now populate the canon details, which requires another look up if the
        # user details are for a synonym etc.
        if taxon.taxon_status in ('accepted', 'doubtful'):
            taxon.is_canon = True
        else:
            taxon.is_canon = False
            # acceptedKey is not provided in the response for canon taxa.
            taxon.canon_usage = self.id_lookup(response['acceptedKey'])

        return taxon


@dataclasses.dataclass()
class Taxa:
    """
           Attempts to load and check the content of the Taxa worksheet. The
        method checks that all taxa have a local name and then validates
        taxon names and parent names against the GBIF backbone. It populates
        two things:
            i)  the taxon_names attribute of the dataset, which is just a set of
                names used as a validation list for taxon names used in data worksheets.
            ii) the taxon_index attribute of the dataset, which contains a set
                of lists recording the full hierarchy of the taxa in the dataset
                for use in dataset searching, so including not just the named taxa,
                but all higher taxa needed to complete the backbone.

                Each list consists of:
                    [worksheet_name (str),
                     gbif_id (int),
                     gbif_parent_id (int),
                     canonical_name (str),
                     taxonomic_rank (str),
                     status (str)]

                Where a taxon is not accepted or doubtful on GBIF, two entries are
                inserted for the taxon, one under the canon name and one under the
                provided name. They will share the same worksheet name and so can
                be paired back up for description generation. The worksheet name
                for parent taxa and deeper taxonomic hierarchy is set to None.

            The index can then be used:
            a) to generate the taxonomic coverage section of the dataset description, and
            b) as the basis of a SQL dataset_taxa table to index the taxonomic coverage
               of datasets.
    """

    taxon_index = []
    taxon_names = set()

    def load(self, worksheet, resources):
        """

        Args:
            worksheet:
            resources:

        Returns:
            Updates the taxon_names and taxon_index attributes of the class instance.

        """

        start_errors = CH.counters['ERROR']

        # Get the field headers from the first row
        LOGGER.info("Reading taxa headers")
        FORMATTER.push()
        dframe = GetDataFrame(worksheet)

        # Dupe headers likely cause serious issues, so stop
        if 'duplicated' in dframe.bad_headers:
            LOGGER.error('Cannot parse taxa with duplicated headers')
            return

        # Get the headers
        headers = IsLower(dframe.headers).values

        # Field cleaning
        core_fields = {'name', 'taxon name', 'taxon type'}
        missing_core = core_fields.difference(headers)

        if missing_core:
            # core names are not found so can't continue
            LOGGER.error('Missing core fields: ', extra={'join': missing_core})
            return

        # Any duplication in names
        dupl_taxon_names = HasDuplicates([dframe.data_columns[headers.index('name')]])

        if dupl_taxon_names:
            LOGGER.error('Duplicated names found: ',
                         extra={'join': dupl_taxon_names.duplicated})

        # Clean string fields of padding
        str_fields = ['name', 'taxon name', 'taxon type', 'parent name', 'parent type']

        for each_field in set(headers).intersection(str_fields):
            idx = headers.index(each_field)
            field_data = IsNotPadded(dframe.data_columns[idx])

            if not field_data:
                LOGGER.error(f'Field {each_field} contains padded strings: ',
                             extra={'join': field_data.failed})
                dframe.data_columns[idx] = field_data.values

        # get dictionaries of the taxa
        taxa = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]

        # check number of taxa found
        if len(taxa) == 0:
            LOGGER.info('No taxon rows found')
            return

        # Standardise to the expected fields, filling in None for any
        # completely missing fields (parent fields could be missing).
        tx_fields = ['name', 'taxon name', 'taxon type', 'taxon id', 'ignore id',
                     'parent name', 'parent type', 'parent id']
        taxa = [{fld: tx.get(fld) for fld in tx_fields} for tx in taxa]

        # Standardize the taxon representation into pairs of taxon and parent dicts.
        # Note that parent tuples cannot have an ignore id.
        #     (name,
        #       (taxon name, taxon type, taxon id, ignore id),
        #       (parent name, parent type, parent id, None))

        taxon_tuples = []

        for row in taxa:
            m_tuple = (row['taxon name'], row['taxon type'], row['taxon id'], row['ignore id'])
            p_tuple = (row['parent name'], row['parent type'], row['parent id'], None)

            # If there is no parent information, replace the parent tuple with None
            if p_tuple == (None, None, None, None):
                p_tuple = None

            self.taxon_names.update(row['name'])
            taxon_tuples.append((row['name'], m_tuple, p_tuple))

        # Get a validator instance
        if resources.use_local_gbif:
            validator = LocalGBIFValidator(resources)
        else:
            validator = RemoteGBIFValidator()

        # Build a dictionary keyed on taxon tuples that stores Taxon
        # objects for parents and children. This is then used to look up
        # combinations of parent and taxon in cross validation

        parent_dict = {}
        taxon_hier = set()

        parents = set([tx[2] for tx in taxon_tuples if tx[2] is not None])

        if len(parents):
            LOGGER.info(f'Checking {len(parents)} parent taxa')
            FORMATTER.push()

            for p_tuple in parents:

                # Create a taxon object
                p_taxon = Taxon(name=p_tuple[0], rank=p_tuple[1],
                                gbif_id=p_tuple[2], ignore_gbif=p_tuple[3])

                # Parents _must_ be a backbone type
                if not p_taxon.is_backbone:
                    LOGGER.error(f'{p_taxon}: is not of a backbone rank')
                    continue

                # Look for a match
                p_taxon = validator.search(p_taxon)

                if not p_taxon.found:
                    LOGGER.error(f'{p_taxon}: {p_taxon.lookup_status}')

                # Store the search results keyed by parent tuple and update hierarchy and index
                parent_dict[p_tuple] = p_taxon
                taxon_hier.update([rw for rw in p_taxon.hierarchy if rw[1] is not None])
                self.taxon_index.append([None, p_taxon.gbif_id, p_taxon.parent_id,
                                         p_taxon.name, p_taxon.rank, p_taxon.taxon_status])

                # If not canon usage, issue a warning using the canon usage
                # stored in the Taxon object and update the index
                if not p_taxon.is_canon:
                    LOGGER.warning(f'{p_taxon.name}: considered a {p_taxon.taxon_status}'
                                   f' of {p_taxon.canon_usage.name} in GBIF backbone')
                    taxon_hier.update([rw for rw in p_taxon.canon_usage.hierarchy if rw[1] is not None])
                    self.taxon_index.append([None,
                                             p_taxon.canon_usage.gbif_id,
                                             p_taxon.canon_usage.parent_id,
                                             p_taxon.canon_usage.name,
                                             p_taxon.canon_usage.rank,
                                             p_taxon.canon_usage.taxon_status])

            FORMATTER.pop()

        # Now check main taxa
        #
        # The taxon index is now populated with parent Taxon objects keyed by
        # data tuples, so now loop over taxon_data to validate the named taxa
        # and then the combinations of taxon and parent status.
        #
        # The combinations are shown below. The taxon row is valid (O) for: a
        # found taxon (with or without a valid parent); a non-matching taxon
        # with a valid parent; and a non-backbone taxon type with a valid
        # parent.
        #
        # Everything else is invalid (X), possibly including a found taxon with
        # a valid parent that isn't actually a parent of the child taxon (case c)
        #
        #                | None  | pr_inv | pr_val |
        # tx_found       |  O a) |  X b)  |  ? c)  |
        # tx_nomatch     | [    X d)    ] |  O e)  |
        # tx_nonbackbone | [    X f)    ] |  O g)  |

        LOGGER.info(f'Validating {len(taxon_tuples)} taxa')
        FORMATTER.push()

        for idx, (m_name, m_tuple, p_tuple) in enumerate(taxon_tuples):

            # Get a row label and retrieve the parent taxon (which is often None)
            rw_ref = f'Row {idx + 2} ({m_name})'
            p_taxon = parent_dict.get(p_tuple)

            # Create the taxon instance
            m_taxon = Taxon(name=m_tuple[0], rank=m_tuple[1],
                            gbif_id=m_tuple[2], ignore_gbif=m_tuple[3])

            # Handle non-backbone cases first - just needs a valid parent.
            if not m_taxon.is_backbone:

                if p_taxon is None or not p_taxon.found:
                    # f) Non backbone with no or bad parent information
                    LOGGER.error(f'{rw_ref}: taxa of type {m_taxon.rank} must '
                                 'have valid parent information.')
                else:
                    # g) Taxon is a non backbone type with good parent info
                    LOGGER.info(f'{rw_ref}: {m_taxon.rank} with valid parent information ')
                    # Update index - no taxon hierarchy except for parent
                    self.taxon_index.append([m_name, -1, p_taxon.gbif_id,
                                             m_taxon.name, m_taxon.rank,
                                             'user'])

                continue

            # Otherwise try and validate backbone taxon

            m_taxon = validator.search(m_taxon)

            # Flag ignored matches
            if ((not m_taxon.is_canon) and
                    (m_taxon.ignore_gbif is not None) and
                    (m_taxon.ignore_gbif == m_taxon.canon_usage.gbif_id)):
                m_taxon.ignore_match = True
            else:
                m_taxon.ignore_match = False

            if m_taxon.found and p_taxon is None:

                # Add the index entry and update hierarchy
                self.taxon_index.append([m_name, m_taxon.gbif_id, m_taxon.parent_id,
                                         m_taxon.name, m_taxon.rank,
                                         m_taxon.taxon_status])
                taxon_hier.update([rw for rw in m_taxon.hierarchy if rw[1] is not None])

                # a) Good backbone with no parent, provide info on taxon status
                if m_taxon.is_canon:
                    LOGGER.info(f'{rw_ref}: in GBIF backbone ({m_taxon.taxon_status})')
                elif m_taxon.ignore_match:
                    LOGGER.error(f'{rw_ref}: ignoring match to canon usage '
                                 f'({m_taxon.canon_usage.name}) but no parent provided')
                else:
                    LOGGER.warn(f'{rw_ref}: considered a {m_taxon.taxon_status} '
                                f'of {m_taxon.canon_usage.name} in GBIF backbone')

                    # Add the canon index entry and update hierarchy
                    self.taxon_index.append([m_name, m_taxon.canon_usage.gbif_id,
                                             m_taxon.canon_usage.parent_id,
                                             m_taxon.canon_usage.name,
                                             m_taxon.canon_usage.rank,
                                             m_taxon.canon_usage.taxon_status])
                    taxon_hier.update([rw for rw in m_taxon.canon_usage.hierarchy if rw[1] is not None])

            elif m_taxon.found and p_taxon is not None:

                if m_taxon.ignore_match and p_taxon.found:

                    # Add the index entry and update hierarchy - record the GBIF
                    # ID but use _specified parent id and hierarchy not canon_
                    self.taxon_index.append([m_name, m_taxon.gbif_id, p_taxon.gbif_id,
                                             m_taxon.name, m_taxon.rank,
                                             'user'])
                    taxon_hier.update([rw for rw in p_taxon.hierarchy if rw[1] is not None])

                    LOGGER.info(f'{rw_ref}: ignoring match to ({m_taxon.canon_usage.name})')

                elif p_taxon.found:
                    # c) Good backbone with good parent - are they compatible? Check if all
                    #    entries in the parent hierarchy appear in the taxon hierarchy
                    parent_hier = set(p_taxon.hierarchy)
                    taxon_hier = set(m_taxon.hierarchy)

                    if not set(parent_hier).issubset(taxon_hier):
                        LOGGER.error(f'{rw_ref}: additional parent information is incompatible')
                    else:
                        LOGGER.info(f'{rw_ref}: in GBIF backbone ({m_taxon.taxon_status}) with '
                                    'compatible parent info')

                    # Add to index and hierarchy
                    self.taxon_index.append([m_name, m_taxon.gbif_id, m_taxon.parent_id,
                                             m_taxon.name, m_taxon.rank,
                                             m_taxon.taxon_status])
                    taxon_hier.update([rw for rw in m_taxon.hierarchy if rw[1] is not None])

                else:
                    # b) Good backbone with bad parent
                    LOGGER.error(f'{rw_ref}: additional parent information is not valid.')

                    # Add to index and hierarchy
                    self.taxon_index.append([m_name, m_taxon.gbif_id, m_taxon.parent_id,
                                             m_taxon.name, m_taxon.rank,
                                             m_taxon.taxon_status])
                    taxon_hier.update([rw for rw in m_taxon.hierarchy if rw[1] is not None])

            elif not m_taxon.found and p_taxon is not None:

                if p_taxon.found:

                    if m_taxon.is_backbone:
                        # e) Taxon is a backbone type that is not present in GBIF
                        #    but the user has provided a valid set of parent taxon
                        #    information.
                        LOGGER.info(f'{rw_ref}: not found in GBIF but has valid parent information')
                    else:
                        # g) Taxon is a non backbone type with good parent info
                        LOGGER.info(f'{rw_ref}: {m_taxon.rank} with valid parent information ')

                    # Add to index  - parent already in hierarchy so nothing to add
                    self.taxon_index.append([m_name, -1, p_taxon.gbif_id,
                                             m_taxon.name, m_taxon.rank,
                                             'user'])

                elif m_taxon.is_backbone and p_taxon is None:
                    # d1) Taxon is a backbone type but is not found in GBIF.
                    if m_taxon.note is None:
                        LOGGER.error(f'{rw_ref}: name and rank combination not found')
                    else:
                        LOGGER.error(f'{rw_ref}: {m_taxon.note}')
                elif m_taxon.is_backbone and not p_taxon.found:
                    # d2) Taxon is a backbone type not found in GBIF and the
                    #     provided parent isn't valid
                    LOGGER.error(f'{rw_ref}): not found in GBIF and has invalid parent information.')
                elif not m_taxon.is_backbone and (p_taxon is None or not p_taxon.found):
                    # f) Non backbone with no or bad parent information
                    LOGGER.error(f'{rw_ref}: taxa of type {m_taxon.rank} must have valid '
                                 'parent information.')

        # Use the taxon hierarchy entries to add higher taxa
        # - drop taxa with a GBIF ID already in the index

        known = [tx[1] for tx in self.taxon_index if tx[1] != -1]
        to_add = [tx for tx in taxon_hier if tx[1] not in known]
        to_add.sort(key=lambda val: BACKBONE_RANKS.index(val[0]))
        LOGGER.info('Indexing taxonomic hierarchy')
        FORMATTER.push()

        # Look up the taxonomic hierarchy
        for tx_lev, tx_id in to_add:
            higher_taxon = validator.id_lookup(tx_id)
            self.taxon_index.append([None,  higher_taxon.gbif_id,
                                     higher_taxon.parent_id,
                                     higher_taxon.name,
                                     higher_taxon.rank,
                                     higher_taxon.taxon_status])
            LOGGER.info(f'Added {tx_lev} {higher_taxon}')

        FORMATTER.pop()

        # summary of processing
        n_errors = CH.counters['ERROR'] - start_errors
        if n_errors > 0:
            LOGGER.info('Taxa contains {} errors'.format(n_errors))
        else:
            LOGGER.info('{} taxa loaded correctly'.format(len(self.taxon_names)))

        FORMATTER.pop()


def taxon_index_to_text(taxa, html=False, depth=4):
    """
    Turns the taxon index from a Taxa instance into a text representation
    of the taxonomic hierarchy used in the dataset.
    """

    lbr = '<br>' if html else '\n'

    def indent(n, html=html):

        ind = '&ensp;' * depth if html else ' ' * depth
        return ind * n

    def format_name(tx, html=html):

        # format the canonical name
        if tx[4] in ['genus', 'species', 'subspecies']:
            if html:
                return f'<i>{tx[3]}</i>'
            else:
                return f'_{tx[3]}_'
        elif tx[4] in ['morphospecies', 'functional group']:
            return f'[{tx[0]}]'
        else:
            return tx[3]

    # Container to hold the output
    html = StringIO()

    # group by parent taxon id in position 2, subsitituting 0 for None
    taxa.sort(key=lambda x: x[2] or 0)
    grouped = {k: list(v) for k, v in groupby(taxa, lambda x: x[2])}

    # start the stack with the kingdoms - these taxa will have None as a parent
    stack = [{'current': grouped[None][0], 'next': grouped[None][1:]}]

    while stack:

        # Handle the current top of the stack: format the canonical name
        current = stack[-1]['current']
        canon_name = format_name(current)

        # Look for a non-None entry in next that shares the same worksheet name
        next_ws_names = [tx[0] for tx in stack[-1]['next']
                         if tx[0] is not None]

        if current[0] in next_ws_names:
            # pop out the matching entry and find which is 'accepted'
            name_pair = stack[-1]['next'].pop(next_ws_names.index(current[0]))
            if current[5] == 'accepted':
                as_name = format_name(name_pair)
                as_status = name_pair[5]
            else:
                as_name = canon_name
                as_status = current[5]
                canon_name = format_name(name_pair)

            txt = f'{indent(len(stack))} {canon_name} (as {as_status}: {as_name}){lbr}'
        else:
            txt = f'{indent(len(stack))} {canon_name}{lbr}'

        html.write(txt)

        # Is this taxon a parent for other taxa - if so add that taxon to the top of
        # the stack, otherwise start looking for a next taxon to push onto the stack.
        # If there is none at the top, pop and look down.
        parent_id = current[1]
        if parent_id in grouped:
            stack.append({'current': grouped[parent_id][0], 'next': grouped[parent_id][1:]})
        else:
            while stack:
                push = stack.pop()
                if push['next']:
                    stack.append({'current': push['next'][0], 'next': push['next'][1:]})
                    break

    return html.getvalue()

import dataclasses
import sqlite3
from collections import Counter
from io import StringIO
from itertools import groupby
from logging import Formatter
from typing import Optional, Union

import requests
from enforce_typing import enforce_types

from safedata_validator.logger import (CH, FORMATTER, LOGGER,
                                       loggerinfo_push_pop)
from safedata_validator.validators import (GetDataFrame, HasDuplicates,
                                           IsLower, IsNotPadded)

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
#        providing parent hierarchy taxon names.  Could define an SQL
#        document to structure the local database to match remote?


class GBIFError(Exception):
    """Exception class for remote GBIF errors

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message="GBIF ID not found"):
        self.message = message
        super().__init__(self.message)


@enforce_types  # However, the checking code should handle bad inputs elegantly.
@dataclasses.dataclass
class Taxon:
    """Holds taxonomic information from a user, which can be populated using GBIF validation.

    There are 3 class properties that can be used to create an instance:
        * name
        * rank
        * gbif_id
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
    is_backbone: bool = dataclasses.field(init=False)
    is_canon: bool = dataclasses.field(init=False)
    canon_usage: 'Taxon' = dataclasses.field(init=False)  # https://stackoverflow.com/questions/33533148
    parent_id: int = dataclasses.field(init=False)
    taxon_status: str = dataclasses.field(init=False)
    lookup_status: str = dataclasses.field(init=False)
    hierarchy: list = dataclasses.field(init=False)

    def __post_init__(self):
        """Sets the defaults for the post-init properties and checks inputs
        """

        if self.gbif_id is not None:
            if isinstance(self.gbif_id, float) and not self.gbif_id.is_integer():
                raise ValueError('GBIF Id is not an integer')
            self.gbif_id = int(self.gbif_id)

        self.rank = self.rank.lower()
        self.is_backbone = self.rank in BACKBONE_RANKS
        self.is_canon = False
        self.canon_usage = None
        self.parent_id = None
        self.taxon_status = None
        self.lookup_status = 'unvalidated'
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

    # TODO - could make _eq_ method and allow Taxon to be used directly
    #        instead of using taxon tuples as data and parent index.


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
        an optional GBIF ID for disambiguation. 

        Args:
            taxon: A Taxon instance

        Returns:
            A Taxon instance
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
                           if ky in taxon_row.keys() and taxon_row[ky] is not None]

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
                      gbif_id=response['key'])

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
                           if ky in response and response[ky] is not None]

        # Now populate the canon details, which requires another look up if the
        # user details are for a synonym etc.
        if taxon.taxon_status in ('accepted', 'doubtful'):
            taxon.is_canon = True
        else:
            taxon.is_canon = False
            # acceptedKey is not provided in the response for canon taxa.
            taxon.canon_usage = self.id_lookup(response['acceptedKey'])

        return taxon


class Taxa:

    def __init__(self, resources):
        """A class to hold a list of taxon names and a validated taxonomic 
        index for those taxa and their taxonomic hierarchy. The validate_taxon
        method checks that taxon details and their optional parent taxon can be
        matched into the the GBIF backbone and populates two things:

        i)  the taxon_names attribute of the dataset, which is just a set of
            names used as a validation list for taxon names used in data worksheets.
        ii) the taxon_index attribute of the dataset, which contains a set
            of lists structured as:

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
        
        The index_higher_taxa method can be used to extend the taxon_index to
        include all of the higher taxa linking the validated taxa.

        The index can then be used:

        a) to generate the taxonomic coverage section of the dataset description, and
        b) as the basis of a SQL dataset_taxa table to index the taxonomic coverage
        of datasets.
        """

        self.taxon_index = []
        self.taxon_names = set()
        self.parents = dict()
        self.hierarchy = set()
        self.n_errors = None
        self.taxon_names_used = set()

        # Get a validator instance
        if resources.use_local_gbif:
            self.validator = LocalGBIFValidator(resources)
        else:
            self.validator = RemoteGBIFValidator()
    
    @loggerinfo_push_pop('Loading Taxa worksheet')
    def load(self, worksheet):
        """Loads a set of taxa from the rows of a SAFE formatted Taxa worksheet and 
        then adds the higher taxa for those rows.

        Args:
            worksheet: An openpyxl worksheet instance following the Taxa formatting

        Returns:
            Updates the taxon_names and taxon_index attributes of the class instance
            using the data in the worksheet.
        """

        start_errors = CH.counters['ERROR']

        # Get the data read in.
        LOGGER.info("Reading taxa data")
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

        # get dictionaries of the taxa
        taxa = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]
        FORMATTER.pop()

        # check number of taxa found
        if len(taxa) == 0:
            LOGGER.info('No taxon rows found')
            return

        # Standardise to the expected fields, filling in None for any
        # completely missing fields (parent fields could be missing).
        tx_fields = ['name', 'taxon name', 'taxon type', 'taxon id', 'ignore id',
                     'parent name', 'parent type', 'parent id']
        taxa = [{fld: tx.get(fld) for fld in tx_fields} for tx in taxa]

        # Standardize the taxon representation into lists of taxon and parent data
        # Note that parent tuples cannot have an ignore id.
        #     [name,
        #       [taxon name, taxon type, taxon id, ignore id],
        #       [parent name, parent type, parent id]]

        for idx, row in enumerate(taxa):
            taxon_info = [row['taxon name'], row['taxon type'], row['taxon id'], row['ignore id']]
            parent_info = [row['parent name'], row['parent type'], row['parent id']]

            # If there is no parent information, replace the parent tuple with None
            if parent_info == [None, None, None]:
                parent_info = None

            self.taxon_names.update([row['name']])
            LOGGER.info(f"Validating row {idx + 1}: {row['name']}")
            FORMATTER.push()
            self.validate_and_add_taxon((row['name'], taxon_info, parent_info))
            FORMATTER.pop()
        
        # Add the higher taxa
        self.index_higher_taxa()

        # summary of processing
        self.n_errors = CH.counters['ERROR'] - start_errors
        if self.n_errors > 0:
            LOGGER.info('Taxa contains {} errors'.format(self.n_errors))
        else:
            LOGGER.info('{} taxa loaded correctly'.format(len(self.taxon_names)))

        FORMATTER.pop()

    # TODO - would be nice to use the decorator, but more complex than
    #        I had anticipated: https://stackoverflow.com/questions/11731136/
    #        Could do this via e.g. @loggerinfo_push_pop(f'Validating
    #        {self._row_description}') but this implementation ties
    #        validate_and_add_taxon() to needing that property populated
    
    def validate_and_add_taxon(self, taxon_input):
        """ Takes user information on a taxon and optionally a parent taxon, 
        validates it and updates the Taxa instance to  include the new details. 
        This is principally used to process rows found in a Taxa worksheet, but
        is deliberately separated out so that a Taxa instance can be populated
        independently of an Excel dataset.

        The taxon_input has the form:

        ['worksheet_name', 
         ['taxon name', 'taxon type', 'taxon id', 'ignore id'],
         ['parent name', 'parent type', 'parent id']]

        If there is no parent information, the structure is:
        ['worksheet_name', 
         ['taxon name', 'taxon type', 'taxon id', 'ignore id'],
          None]

        Args:
            taxon_input: Taxon information in standard form as above

        Returns:
            Updates the taxon_names and taxon_index attributes of the class instance.
        """

        m_name, taxon_info, parent_info = taxon_input

        # Sanitise worksheet names for taxa - only keep unpadded strings.
        if m_name is None or not isinstance(m_name, str) or m_name.isspace():
            LOGGER.error('Worksheet name missing, whitespace only or not text')
        elif m_name != m_name.strip():
            LOGGER.error(f"Worksheet name has whitespace padding: {repr(m_name)}")
            m_name = m_name.strip()
            self.taxon_names.add(m_name)
        else:
            self.taxon_names.add(m_name)
        
        # Check the parent details
        p_fail = False
        if parent_info is not None:
            # Name and rank must be unpadded strings - can still check cleaned padded strings
            for idx, idx_name in ((0, 'Parent name'), (1, 'Parent rank')):
                val = parent_info[idx]

                if val is None or not isinstance(val, str):
                    LOGGER.error(f'{idx_name} missing or not text')
                    p_fail = True
                elif val != val.strip():
                    LOGGER.error(f"{idx_name} has whitespace padding: {repr(val)}")
                    parent_info[idx] = val.strip()

            # ID can be None or an integer (openpyxl loads all values as float)
            if not(parent_info[2] is None or 
                   (isinstance(parent_info[2], float) and parent_info[2].is_integer()) or 
                   isinstance(parent_info[2], int)) :
                LOGGER.error('Parent GBIF ID contains value that is not an integer')
                p_fail = True
        
        # Check the main taxon details
        mfail = False

        # Name and rank must be unpadded strings - can still check cleaned padded strings
        for idx, idx_name in ((0, 'Taxon name'), (1, 'Taxon rank')):
            val = taxon_info[idx]

            if val is None or not isinstance(val, str) or val.isspace():
                LOGGER.error(f'{idx_name} missing, whitespace only or not text')
                mfail = True
            elif val != val.strip():
                LOGGER.error(f"{idx_name} has whitespace padding: {repr(val)}")
                taxon_info[idx] = val.strip()

        # GBIF ID and Ignore ID can be None or an integer (openpyxl loads all values as float)
        for idx, idx_name in ((2, 'GBIF ID'), (3, 'Ignore ID')):
            val = taxon_info[idx]
        
            if not(val is None or 
                   (isinstance(val, float) and val.is_integer()) or 
                   isinstance(val, int)) :
                LOGGER.error(f'{idx_name} contains value that is not an integer: {val}')
                mfail = True
        
        if p_fail:
            LOGGER.error(f'Parent taxon details not properly formatted, cannot validate')
        
        if mfail:
            LOGGER.error(f'Taxon details not properly formatted, cannot validate')
        
        if mfail or p_fail:
            return

        # Now that inputs are sanitised, continue with checking...
        # Parent taxon checking - can be None, already processed with a previous
        # information and stored in the parent index using a tuple of the parent
        # as a key, or be new and need processing.
        if parent_info is None:
            p_taxon = None
        elif tuple(parent_info) in self.parents:
            p_taxon = self.parents[tuple(parent_info)]
        else:
            # Create a taxon object
            p_taxon = Taxon(name=parent_info[0], rank=parent_info[1], gbif_id=parent_info[2])

            # Look for a match
            if p_taxon.is_backbone:
                p_taxon = self.validator.search(p_taxon)

                # Update the hierarchy and index with the search results
                self.hierarchy.update([rw for rw in p_taxon.hierarchy if rw[1] is not None])
                self.taxon_index.append([None, p_taxon.gbif_id, p_taxon.parent_id,
                                            p_taxon.name, p_taxon.rank, p_taxon.taxon_status])

                if p_taxon.is_backbone and p_taxon.found and not p_taxon.is_canon:
                    self.hierarchy.update([rw for rw in p_taxon.canon_usage.hierarchy if rw[1] is not None])
                    self.taxon_index.append([None,
                                            p_taxon.canon_usage.gbif_id,
                                            p_taxon.canon_usage.parent_id,
                                            p_taxon.canon_usage.name,
                                            p_taxon.canon_usage.rank,
                                            p_taxon.canon_usage.taxon_status])
            
            # Store the parent taxon keyed by parent information (needs tuple)
            self.parents[tuple(parent_info)] = p_taxon

        # Report on the parent information
        if p_taxon is not None:
            if not p_taxon.is_backbone:
                LOGGER.error(f'Parent taxon ({p_taxon.name}) is not of a backbone rank')
                
            elif not p_taxon.found:
                LOGGER.error(f'Parent taxon ({p_taxon.name}) {p_taxon.lookup_status}')

            elif not p_taxon.is_canon:
                LOGGER.warning(f'Parent taxon ({p_taxon.name}) considered a {p_taxon.taxon_status}'
                               f' of {p_taxon.canon_usage.name} in GBIF backbone')
            else:
                LOGGER.info(f'Parent taxon ({p_taxon.name}) accepted')
        else:
                LOGGER.info('No parent taxon provided')

        # Now check main taxa
        #
        # The parent list is now populated with parent Taxon objects keyed by
        # data tuples, so now loop over taxon_data to validate the named taxa
        # and then the combinations of taxon and parent status.
        #
        # The combinations are shown below. The taxon row is valid (O) for: a
        # found taxon (with or without a valid parent); a non-matching taxon
        # with a valid parent; a non-backbone taxon type with a valid
        # parent; and a backbone taxon set to ignore the match with a valid 
        # parent.
        #
        # Everything else is invalid (X), possibly including a found taxon with
        # a valid parent that isn't actually a parent of the child taxon 
        #
        #                | None  | pr_inv | pr_val |
        # tx_ignore      |  X    |  X     |  O     |
        # tx_found       |  O    |  X     |  ?     |
        # tx_nomatch     |  X    |  X     |  O     |
        # tx_nonbackbone |  X    |  X     |  O     |

        # Create the taxon instance
        m_taxon = Taxon(name=taxon_info[0], rank=taxon_info[1], gbif_id=taxon_info[2])
        ignore_gbif = taxon_info[3]

        if ignore_gbif is not None:
            # Handle ignored matches first

            # The taxon must be a backbone taxon - can't ignore impossible matches
            if not m_taxon.is_backbone:
                LOGGER.error('Ignore ID can only be used with GBIF backbone taxon ranks')
            else:
                # It should also be found and the ignore ID should match to the actual usage
                # or canon usage.
                m_taxon = self.validator.search(m_taxon)

                if not m_taxon.found:
                    LOGGER.error('Taxon with Ignore ID not found in GBIF backbone')
                elif m_taxon.is_canon and (m_taxon.gbif_id != ignore_gbif):
                    LOGGER.error(f'Ignore ID does not match the canon GBIF usage ({m_taxon.gbif_id})')
                elif not m_taxon.is_canon and (m_taxon.canon_usage.gbif_id != ignore_gbif):
                    LOGGER.error(f'Taxon is non-canon and Ignore ID does not match the canon GBIF usage ({m_taxon.canon_usage.gbif_id})')
                else:
                    LOGGER.info('Canon GBIF usage ignored')

            # It must also have a valid parent.
            if p_taxon is None:
                LOGGER.error('Taxa with Ignore ID must provide parent information.')
            elif not p_taxon.found:
                LOGGER.error('Taxon with Ignore ID has invalid parent information.')
            else:
                LOGGER.info('Taxon with ignored canon usage has valid parent information.')
                # Update index - no taxon hierarchy except for parent
                self.taxon_index.append([m_name, -1, p_taxon.gbif_id,
                                            m_taxon.name, m_taxon.rank,
                                            'user'])

        elif not m_taxon.is_backbone:

            # Now handle non-backbone cases - just needs a valid parent.
            if p_taxon is None:
                LOGGER.error(f'Taxon of type {m_taxon.rank} must provide parent information.')
            elif not p_taxon.found:
                # Non backbone with bad parent information
                LOGGER.error(f'Taxon of type {m_taxon.rank} has invalid parent information.')
            else:
                # Non backbone with with good parent info
                LOGGER.info(f'Taxon of type {m_taxon.rank} has valid parent information')
                # Update index - no taxon hierarchy except for parent
                self.taxon_index.append([m_name, -1, p_taxon.gbif_id,
                                            m_taxon.name, m_taxon.rank,
                                            'user'])
        
        else:
            # Otherwise try and validate backbone taxon
            m_taxon = self.validator.search(m_taxon)

            if m_taxon.found and p_taxon is None:

                # Add the index entry and update hierarchy
                self.taxon_index.append([m_name, m_taxon.gbif_id, m_taxon.parent_id,
                                         m_taxon.name, m_taxon.rank,
                                         m_taxon.taxon_status])
                
                self.hierarchy.update([rw for rw in m_taxon.hierarchy if rw[1] is not None])

                # Good backbone with no parent, provide info on taxon status
                if m_taxon.is_canon:
                    LOGGER.info(f'Taxon found in GBIF backbone ({m_taxon.taxon_status})')
                else:
                    LOGGER.warning(f'Taxon considered a {m_taxon.taxon_status} '
                                f'of {m_taxon.canon_usage.name} in GBIF backbone')

                    # Add the canon index entry and update hierarchy
                    self.taxon_index.append([m_name, m_taxon.canon_usage.gbif_id,
                                             m_taxon.canon_usage.parent_id,
                                             m_taxon.canon_usage.name,
                                             m_taxon.canon_usage.rank,
                                             m_taxon.canon_usage.taxon_status])
                    self.hierarchy.update([rw for rw in m_taxon.canon_usage.hierarchy if rw[1] is not None])

            elif m_taxon.found and p_taxon is not None:

                if p_taxon.found:

                    # Good backbone with good parent - are they compatible? Check if all
                    # entries in the parent hierarchy appear in the taxon hierarchy
                    if not set(p_taxon.hierarchy).issubset(m_taxon.hierarchy):
                        LOGGER.error(f'Taxon in GBIF backbone ({m_taxon.taxon_status}) with incompatible parent information')
                    else:
                        LOGGER.info(f'Taxon in GBIF backbone ({m_taxon.taxon_status}) with compatible parent information')

                else:
                    # Good backbone with bad parent
                    LOGGER.error(f'Taxon in GBIF backbone ({m_taxon.taxon_status}) but with invalid parent information.')

                # Add to index and hierarchy
                self.taxon_index.append([m_name, m_taxon.gbif_id, m_taxon.parent_id,
                                            m_taxon.name, m_taxon.rank,
                                            m_taxon.taxon_status])
                self.hierarchy.update([rw for rw in m_taxon.hierarchy if rw[1] is not None])

            elif not m_taxon.found:

                if p_taxon is None:
                    # Taxon is a backbone type but is not found in GBIF and has no parent info
                    if m_taxon.lookup_status == 'No match found':
                        LOGGER.error('Taxon name and rank combination not found')
                    else:
                        LOGGER.error(f'GBIF issue: {m_taxon.lookup_status}')

                elif not p_taxon.found:
                    # Taxon is a backbone type but not found and parent not found either
                    LOGGER.error(f'Taxon not found in GBIF and has invalid parent information.')
                else:
                    # Taxon is a backbone type but not found but does have valid parent info
                    LOGGER.info('Taxon not found in GBIF but has valid parent information')
                    
                    # Add to index  - parent already in hierarchy so nothing to add
                    self.taxon_index.append([m_name, -1, p_taxon.gbif_id,
                                             m_taxon.name, m_taxon.rank,
                                             'user'])

    @loggerinfo_push_pop('Indexing taxonomic hierarchy')
    def index_higher_taxa(self):

        # Use the taxon hierarchy entries to add higher taxa
        # - drop taxa with a GBIF ID already in the index

        known = [tx[1] for tx in self.taxon_index if tx[1] != -1]
        to_add = [tx for tx in self.hierarchy if tx[1] not in known]
        to_add.sort(key=lambda val: BACKBONE_RANKS.index(val[0]))

        # Look up the taxonomic hierarchy
        for tx_lev, tx_id in to_add:
            higher_taxon = self.validator.id_lookup(tx_id)
            self.taxon_index.append([None,  
                                     higher_taxon.gbif_id,
                                     higher_taxon.parent_id,
                                     higher_taxon.name,
                                     higher_taxon.rank,
                                     higher_taxon.taxon_status])
            LOGGER.info(f'Added {tx_lev} {higher_taxon}')

    @property
    def is_empty(self):
        return len(self.taxon_names) == 0

def taxon_index_to_text(taxon_index, html=False, indent_width=4):
    """
    Turns the taxon index from a Taxa instance into a text representation
    of the taxonomic hierarchy used in the dataset.
    """

    lbr = '<br>' if html else '\n'

    def indent(n, use_html=html):

        ind = '&ensp;' if use_html else ' '
        return ind * indent_width * (n - 1)

    def format_name(tx, use_html=html):

        # format the canonical name
        if tx[4] in ['genus', 'species', 'subspecies']:
            if use_html:
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
    taxon_index.sort(key=lambda x: x[2] or 0)
    grouped = {k: list(v) for k, v in groupby(taxon_index, lambda x: x[2])}

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

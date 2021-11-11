import os
import sqlite3
import requests
from collections import Counter
import xlrd
from typing import Union, Optional
import dataclasses
from enforce_typing import enforce_types

from .logger import LOGGER
from .validators import (ToLower, IsUnique, BlankToNone, CellToValue, 
                         IsPadded, IsNotBlank, AllNone)

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
#  sources had more similar structure. Could have a row return method
#  that is used within a generic id_lookup and search class.
#  Notably there us much less information in the row return from the
#  local SQL db (e.g. merging acceptedKey and parentKey, but also not
#  providing parent hierarcy taxon names.  Could defining some SQL
#  queries to add some of that information into the local database?

class GBIFError(Exception):
    """Exception class for GBIF errors

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
    canon_id: int = dataclasses.field(init=False)
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
        self.canon_id = None
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

    def __init__(self, gbif_database: str):
        # TODO - remove use of LOGGER within these functions. More portable

        LOGGER.info('Validating local GBIF database: ' + gbif_database)

        # Does the provided path exist and is it a functional SQLite database
        # with a backbone table? Because sqlite3 can connect to any existing path,
        # use a query attempt to reveal exceptions

        if not os.path.exists(gbif_database):
            LOGGER.critical('Local GBIF database not found')
            raise IOError('Local GBIF database not found')
        try:
            conn = sqlite3.connect(gbif_database)
            _ = conn.execute('select count(*) from backbone;')
        except sqlite3.OperationalError:
            LOGGER.critical('Local GBIF database does not contain the backbone table')
            raise IOError('Local GBIF database does not contain the backbone table')
        except sqlite3.DatabaseError:
            LOGGER.critical('Local SQLite database not valid')
            raise IOError('Local SQLite database not valid')
        else:
            self.gbif_conn = conn
            self.gbif_conn.row_factory = sqlite3.Row

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
            raise GBIFError('Cannot validate non-backbone taxa')

        if taxon.gbif_id is not None:

            # get the record associated with the provided ID
            id_taxon = self.id_lookup(taxon.gbif_id)

            # Check that name and rank are congruent with id
            if (id_taxon.name != taxon.name) or (id_taxon.rank != taxon.rank):
                raise GBIFError('ID does not match name and rank')

            return id_taxon

        else:
            # get the set of records associated with the taxon and rank

            sql = (f"select * from backbone where canonical_name ='{taxon.name}' "
                   f"and rank= '{taxon.rank.upper()}';")
            taxon_rows = self.gbif_conn.execute(sql).fetchall()
            selected_row = None

            if len(taxon_rows) == 0:
                # No matching rows
                raise GBIFError('No match found')
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
                raise GBIFError(f'Multiple equal matches for {taxon.name}')

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

        # Add the taxonomic hierarchy, using a mapping of backbone ranks (except
        # subspecies) to backbone table fields.
        taxon.hierarchy = [(rk, taxon_row[ky])
                           for rk, ky in [(r, r + '_key') for r in BACKBONE_RANKS[:-1]]
                           if ky in taxon_row.keys()]

        # parent key in the local database has the odd property that the parent
        # tax_gbif['parent_key']dual duty: points up to parent for canon taxa
        # and 'up' to canon for non-canon taxa.
        if taxon.taxon_status in ['accepted', 'doubtful']:
            taxon.is_canon = True
            taxon.parent_id = taxon_row['parent_key']
        else:
            taxon.is_canon = False
            taxon.canon_id = taxon_row['parent_key']

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
                raise GBIFError('ID does not match name and rank')

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
                    raise GBIFError(response['note'])
                else:
                    raise GBIFError('No match found')

            return self.id_lookup(response['usageKey'])

    @staticmethod
    def id_lookup(gbif_id: int):
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
            taxon.canon_id = response['acceptedKey']

        return taxon


# class Taxa:
#
#     """
#     Validates a list of the taxa defined in a dataset and provides an
#     interface to the taxon names and taxon index
#
#     :return:
#     """
#
#     def __init__(self):
#         """
#         Creates a Taxa instance from a Dataset config
#         # - If gbif_database isn't provided to __init__ or in the config then use the
#         #   online API. Otherwise try and use the provided file, preferring the __init__
#         #   value over the config value
#
#
#         :param config:
#         """
#
#         self.taxon_data = []
#         self.taxon_index = {}
#         self.taxon_names = set()
#         self.taxon_hierarchy = set()
#
#     def load_taxa(self, workbook):
#         """
#         Loads the taxa from a taxon worksheet and populates the taxa_dict attribute.
#
#         :param workbook: An xlrd.Workbook instance
#         :return:
#         """
#
#         try:
#             worksheet = workbook.sheet_by_name('Taxa')
#         except xlrd.XLRDError:
#             # This might mean that the study doesn't have any taxa, so return an empty
#             # set. If the datasets then contain taxonomic names, it'll fail gracefully.
#             LOGGER.warn("No taxa worksheet found - assuming no taxa in data for now!")
#             return
#
#         # Get and check the headers
#         worksheet_rows = worksheet.get_rows()
#         headers = ToLower(CellToValue(next(worksheet_rows)))
#
#         # check which fields are found
#         if 'name' not in headers:
#             # dataset names are not found so can't continue
#             LOGGER.error('No name column found - no further checking')
#             return
#
#         # Duplicated headers are a problem in that it will cause values in the taxon
#         # dictionaries to be overwritten. Depending on what gets overwritten, this can
#         # produce really unpredictable bugs, so just stop here.
#         headers_unique = IsUnique(headers)
#         if not headers_unique:
#             LOGGER.error('Duplicated column headers in Taxa worksheet: ',
#                          extra={'join': headers_unique.invalid})
#             return
#
#         # Load dictionaries of the taxa, strip out any rows that consist
#         # of nothing but empty cells and replace empty cells with None
#         taxa = list(worksheet_rows)
#         taxa = [BlankToNone(CellToValue(row, none_ok=True)) for row in taxa]
#         taxa = [dict(zip(headers, row)) for row in taxa]
#
#         # check number of taxa found and standardise names
#         if len(taxa) == 0:
#             LOGGER.info('No taxon rows found')
#             return
#
#         # Standardise to the core fields, filling in None for any
#         # completely missing fields (parent fields could be missing).
#         tx_fields = ['name', 'taxon name', 'taxon type', 'taxon id', 'ignore id',
#                      'parent name', 'parent type', 'parent id']
#         taxa = [{fld: tx.get(fld) for fld in tx_fields} for tx in taxa]
#
#         # Clean the contents of whitespace padding - TODO fix here - not stripping...
#         for idx, tx in enumerate(taxa):
#             tx_pad = IsPadded(list(tx.values()), none_ok=True)
#             if tx_pad:
#                 LOGGER.error(f'Whitespace padding in row {idx + 2}: {",".join(tx_pad.valid)}')
#                 tx = {tx: vl.strip() for tx, vl in tx}
#
#         # Validate the names used within the data
#         taxon_names = IsNotBlank([tx['name'] for tx in taxa])
#         if taxon_names.invalid:
#             LOGGER.error('Blank entries in name column')
#
#         # Any duplication in cleaned names
#         taxon_names = IsUnique(taxon_names.valid)
#         if taxon_names.invalid:
#             LOGGER.error('Duplicated names found: ', extra={'join': taxon_names.invalid})
#
#         # set the unique taxon names for the dataset instance
#         self.taxon_names = set(taxon_names.valid)
#
#         # Standardize the taxon representation into pairs of taxon and parent dicts.
#         # Note that parent tuples cannot have an ignore id.
#         #     (name,
#         #       (taxon name, taxon type, taxon id, ignore id),
#         #       (parent name, parent type, parent id, None))
#         for row in taxa:
#             m_tuple = (row['taxon name'], row['taxon type'], row['taxon id'], row['ignore id'])
#             p_tuple = (row['parent name'], row['parent type'], row['parent id'], None)
#
#             # If there is no parent information, replace the parent tuple with None
#             if AllNone(p_tuple, none_ok=True):
#                 p_tuple = None
#
#             self.taxon_data.append((row['name'], m_tuple, p_tuple))
#
#     def validate(self, validator):
#         """
#         Validates the taxon data loaded from the worksheet using an instance
#         of one of the two GBIFValidator classes.
#
#         :param validator: A GBIFValidator instance
#         :return:
#         """
#
#         # Build a dictionary keyed on taxon tuples that stores validated Taxon
#         # objects for each taxon. Start with the parents.
#         parents = set([tx[2] for tx in self.taxon_data if tx[2] is not None])
#
#         if len(parents):
#             LOGGER.info(f'Checking {len(parents)} parent taxa',
#                         extra={'indent_after': 2})
#
#         for p_tuple in parents:
#
#             p_taxon = Taxon(name=p_tuple[0], rank=p_tuple[1],
#                             gbif_id=p_tuple[2], ignore_gbif=p_tuple[3])
#
#             # Parents _must_ be a backbone type
#             if not p_taxon.is_backbone:
#                 LOGGER.error(f'{p_taxon}: is not of a backbone rank')
#             else:
#                 try:
#                     p_taxon = validator.search(p_taxon)
#                 except GBIFError as err:
#                     LOGGER.error(f'{p_taxon}: {err}')
#
#                 self.taxon_index[p_tuple] = p_taxon
#
#             # Non canon
#             if not p_taxon.is_canon:
#                 c_taxon = validator.id_lookup(p_taxon.canon_id)
#                 LOGGER.warn(f'{p_taxon.name}: considered a {p_taxon.taxon_status}'
#                             f' of {c_taxon.name} in GBIF backbone')
#                 self.taxon_index[p_tuple + tuple(['canon'])] = p_taxon
#
#         # Now check main taxa
#         #
#         # The taxon index is now populated with parent Taxon objects keyed by
#         # data tuples, so now loop over taxon_data to validate the named taxa
#         # and then the combinations of taxon and parent status.
#         #
#         # The combinations are shown below. The taxon row is valid (O) for: a
#         # found taxon (with or without a valid parent); a non-matching taxon
#         # with a valid parent; and a non-backbone taxon type with a valid
#         # parent.
#         #
#         # Everything else is invalid (X), possibly including a found taxon with
#         # a valid parent that isn't actually a parent of the child taxon (case c)
#         #
#         #                | None  | pr_inv | pr_val |
#         # tx_found       |  O a) |  X b)  |  ? c)  |
#         # tx_nomatch     | [    X d)    ] |  O e)  |
#         # tx_nonbackbone | [    X f)    ] |  O g)  |
#
#         LOGGER.info(f'Validating {len(self.taxon_data)} taxa',
#                     extra={'indent_before': 1, 'indent_after': 2})
#
#         for idx, (m_name, m_tuple, p_tuple) in enumerate(self.taxon_data):
#
#             # Get a row label and retrieve the parent taxon (which is often None)
#             rw_ref = f'Row {idx + 2} ({m_name})'
#             p_taxon = self.taxon_index.get(p_tuple)
#
#             # Create the taxon
#             m_taxon = Taxon(name=m_tuple[0], rank=m_tuple[1],
#                             gbif_id=m_tuple[2], ignore_gbif=m_tuple[3])
#
#             if not m_taxon.is_backbone:
#                 # Handle non-backbone cases first - just needs a valid parent.
#
#                 if p_taxon is None or not p_taxon.found:
#                     # f) Non backbone with no or bad parent information
#                     LOGGER.error(f'{rw_ref}: taxa of type {m_taxon.user.rank} must have valid '
#                                   'parent information.')
#                 else:
#                     # g) Taxon is a non backbone type with good parent info
#                     LOGGER.info(f'{rw_ref}: {m_taxon.rank} with valid parent information ')
#             else:
#                 # Try and validate backbone taxon, capturing any validation errors
#                 search_error = None
#                 c_taxon = None
#
#                 try:
#                     m_taxon = validator.search(m_taxon)
#                     if not m_taxon.is_canon:
#                         c_taxon = validator.id_lookup(m_taxon.canon_id)
#                 except GBIFError as err:
#                     search_error = err
#
#                 if m_taxon.found and p_taxon is None:
#
#                     # a) Good backbone with no parent, provide info on taxon status
#                     if m_taxon.is_canon:
#                         LOGGER.info(f'{rw_ref}: in GBIF backbone ({m_taxon.taxon_status})')
#                     else:
#                         LOGGER.warn(f'{rw_ref}: considered a {m_taxon.taxon_status} '
#                                     f'of {c_taxon.name} in GBIF backbone')
#
#                     if p_taxon is not None:
#                         if p_taxon.found:
#                             # c) Good backbone with good parent - are they compatible? Check if all
#                             #    entries in the parent hierarchy appear in the taxon hierarchy
#                             parent_hier = set(p_taxon.hierarchy)
#                             taxon_hier = set(m_taxon.hierarchy)
#                             if not set(parent_hier).issubset(taxon_hier):
#                                 LOGGER.error(f'{rw_ref}: additional parent information is incompatible')
#                         else:
#                             # b) Good backbone with bad parent
#                             LOGGER.error(f'{rw_ref}: additional parent information is not valid.')
#
#                 # update hierarchy and index
#                 self.taxon_index[m_tuple] = m_taxon
#                 self.taxon_hierarchy.update(m_taxon.hierarchy)
#
#             elif not m_taxon.found and p_taxon is not None:
#
#                 if p_taxon.found:
#
#                     if m_taxon.is_backbone:
#                         # e) Taxon is a backbone type that is not present in GBIF
#                         #    but the user has provided a valid set of parent taxon
#                         #    information.
#                         LOGGER.info(f'{rw_ref}: not found in GBIF but has valid parent information')
#                     else:
#                         # g) Taxon is a non backbone type with good parent info
#                         LOGGER.info(f'{rw_ref}: {m_taxon.rank} with valid parent information ')
#
#                     # update hierarchy
#                     self.taxon_hierarchy.update(p_taxon.hierarchy)
#
#                 elif m_taxon.is_backbone and p_taxon is None:
#                     # d1) Taxon is a backbone type but is not found in GBIF.
#                     if m_taxon.note is None:
#                         LOGGER.error(f'{rw_ref}: name and rank combination not found')
#                     else:
#                         LOGGER.error(f'{rw_ref}: {m_taxon.note}')
#                 elif m_taxon.is_backbone and not p_taxon.found:
#                     # d2) Taxon is a backbone type not found in GBIF and the
#                     # provided parent isn't valid
#                     LOGGER.error(f'{rw_ref}): not found in GBIF and has invalid parent information.')
#                 elif not m_taxon.is_backbone and (p_taxon is None or not p_taxon.found):
#                     # f) Non backbone with no or bad parent information
#                     LOGGER.error(f'{rw_ref}: taxa of type {m_taxon.user.rank} must have valid '
#                                  'parent information.')
#
#     def build_index(self):
#
#         index = []
#
#         for taxon in self.taxon_index.values():
#             if taxon.found:
#                 if taxon.is_canon:
#                         index.extend([[taxon.user.gbif_id, taxon.parent_id, taxon.user.name,
#                                        taxon.user.rank, taxon.user.taxon_status]])
#                 else:
#                     index.extend([[taxon.user.gbif_id, taxon.parent_id, taxon.user.name,
#                                    taxon.user.rank, taxon.user.taxon_status],
#                                   [taxon.canon.gbif_id, taxon.parent_id, taxon.canon.name,
#                                    taxon.canon.rank, taxon.canon.taxon_status]])
#             else:
#                 index.extend([[-1, None, taxon.user.name, taxon.user.rank, 'user']])
#
#         # Look up the unique taxon hierarchy entries
#         # - drop taxa with a GBIF ID already in the index
#         indexed = [tx[1] for tx in self.taxon_index.values()]
#         taxon_hierarchy = {tx for tx in self.taxon_hierarchy if tx[1] not in indexed}
#
#         # - sort into ascending taxonomic order
#         taxon_hierarchy = list(taxon_hierarchy)
#         taxon_hierarchy.sort(key=lambda val: BACKBONE_RANKS.index(val[0]))
#         LOGGER.info('Indexing taxonomic hierarchy', extra={'indent_before': 1, 'indent_after': 2})
#
#         # Look up the taxonomic hierarchy
#         for tx_lev, tx_id in taxon_hierarchy:
#             taxon = Taxon(gbif_id=tx_id)
#             validator.validate(taxon)
#             self.taxon_index.append([None] + canon)
#             LOGGER.info('Added {}'.format(canon[2]))
#
#         # summary of processing
#         n_errors = CH.counters['ERROR'] - start_errors
#         if n_errors > 0:
#             LOGGER.info('Taxa contains {} errors'.format(n_errors),
#                         extra={'indent_before': 1})
#         else:
#             LOGGER.info('{} taxa loaded correctly'.format(len(self.taxon_names)),
#                         extra={'indent_before': 1})

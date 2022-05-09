"""## The taxa submodule
This module describes classes and methods used to compile taxonomic data from
datasets and to validate taxonomy against the GBIF backbone database and/or the
NCBI taxonomy database.

The two parallel Taxon dataclasses (GBIFTaxon and NCBITaxon) are used to store
data about a taxon entry in a dataset. They are initialised with user data and
then the relevant taxon Validator classes (GBIF or NCBI) can be used to update a
Taxon object with the result of validaton against a particular database. Online
GBIF validation can be performed using the online API (`RemoteGBIFValidator`),
whereas online NCBI validation makes use of the Entrez Eutils (`RemoteNCBIValidator`).
Faster validation may be performed using local copies of the databases
(`LocalGBIFValidator` and `LocalNCBIValidator`).

Parallel 'Taxa' worksheets (GBIFTaxa and NCBITaxa) are defined, which are used to
load and collate the set of taxonomic entries from a dataset.

Note that we explicitly exclude form and variety from the set of GBIF backbone
taxonomic levels because they cannot be matched into the backbone hierarchy
without extra API calls.

When validating against the NCBI database supplied taxa of any rank (i.e. strain
or clade) which can be successfully validated will be recorded. However, associated
higher taxa will only be recorded if their ranks are either a GBIF backbone rank
or superkingdom.
"""

import dataclasses
import sqlite3
from collections import Counter
from io import StringIO
from itertools import groupby
from logging import Formatter
from typing import Optional, Union
import time
from lxml import etree

import requests
from enforce_typing import enforce_types

from safedata_validator.resources import Resources
from safedata_validator.logger import (COUNTER_HANDLER, FORMATTER, LOGGER,
                                       loggerinfo_push_pop)
from safedata_validator.validators import (GetDataFrame, HasDuplicates,
                                           IsLower, blank_value)

BACKBONE_RANKS = ['kingdom', 'phylum', 'order', 'class', 'family',
                  'genus', 'species', 'subspecies']

# Extended version of backbone ranks to capture superkingdoms
BACKBONE_RANKS_EX = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                    'family', 'genus', 'species', 'subspecies']

# TODO - Modify the resource file to ask the user to provide an email address
# This should only be done if the user actually wants to use this module as it
# isn't need elsewhere (as far as I know)
# Also should ask for api key

# TODO - Work out how this best integrates with taxa
# Lot of different potential approaches (e.g. making making NCBI and GBIF validators
# sub-classes of a validator class). Need to establish what the most extensible
# and stable option is

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

class NCBIError(Exception):
    """Exception class for remote NCBI errors

    Attributes:
        message: explanation of the error
    """

    def __init__(self, message="No entry found for ID, probably bad ID"):
        self.message = message
        super().__init__(self.message)

@enforce_types  # However, the checking code should handle bad inputs elegantly.
@dataclasses.dataclass
class GBIFTaxon:
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
        * canon_usage: a GBIFTaxon instance holding the canonical usage for the taxon
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

    # TODO - could make _eq_ method and allow GBIFTaxon to be used directly
    #        instead of using taxon tuples as data and parent index.

@enforce_types
@dataclasses.dataclass
class NCBITaxon:
    """Holds taxonomic information from a user on a microbial taxa. This can be
    populated using NCBI validation.

    There are 3 class properties that can be used to create an instance:
        * name
        * taxa_hier: Dictionary of valid taxonomic hierarchy (with ID's)
        * ncbi_id: NCBI ID for full taxa (i.e. including non-backbone ranks)
    The remaining properties are populated by processing functions not when an
    instance is created.
        * superseed: is supplied taxon name/ID still the accepted usage
        * orig: Records rank of the orginal taxa if it was not found
    """

    # Init properties
    name: str
    rank: str
    taxa_hier: dict
    ncbi_id: Optional[Union[int, float]] = None
    superseed: bool = dataclasses.field(init=False)
    orig: str = dataclasses.field(init=False)

    def __post_init__(self):
        """Sets the defaults for the post-init properties and checks inputs
        """

        if self.ncbi_id is not None:
            if isinstance(self.ncbi_id, float) and not self.ncbi_id.is_integer():
                raise TypeError('NCBI Id is not an integer')
            self.ncbi_id = int(self.ncbi_id)

        if self.taxa_hier is not None:
            if len(self.taxa_hier) == 0:
                raise ValueError('Taxa hierarchy dictonary empty')
            elif all(isinstance(x,str) for x in self.taxa_hier.keys()) == False:
                raise ValueError('Not all taxa dictionary keys are strings')
            elif all(isinstance(x,tuple) for x in self.taxa_hier.values()) == False:
                raise ValueError('Not all taxa dictionary values are tuples')
            elif all(list(map(type,x)) == [str, int, int] or list(map(type,x)) ==
                 [str, int, type(None)] for x in self.taxa_hier.values()) == False:
                 raise ValueError('Taxa tuples not all in [string integer integer] form')

        if self.rank.lower() != list(self.taxa_hier.keys())[-1]:
            raise ValueError(f'Provided rank ({self.rank.lower()}) does not match'
                             f' lowest rank in supplied hiearchy ('
                             f'{list(self.taxa_hier.keys())[-1]})')

        if self.name != (list(self.taxa_hier.values())[-1])[0]:
            raise ValueError(f'Provided taxon name ({self.name}) does not match'
                             f' name of the lowest entry in supplied hiearchy ('
                             f'{(list(self.taxa_hier.values())[-1])[0]})')

        if self.ncbi_id != None and self.ncbi_id != (list(self.taxa_hier.values())[-1])[1]:
            raise ValueError(f'Provided NCBI ID ({self.ncbi_id}) does not match'
                             f' ID of the lowest entry in supplied hiearchy ('
                             f'{(list(self.taxa_hier.values())[-1])[1]})')

        self.rank = self.rank.lower()
        self.superseed = False
        self.orig = None

    def __repr__(self):

        if self.orig != None:
            return f"{self.name} (resolved as {self.rank} rather than {self.orig})"
        elif self.superseed:
            return f"{self.name} (superseeded taxon details provided)"
        else:
            return f"{self.name}"

class LocalGBIFValidator:

    def __init__(self, resources):

        conn = sqlite3.connect(resources.gbif_database)
        conn.row_factory = sqlite3.Row
        self.gbif_conn = conn

    def __del__(self):

        self.gbif_conn.close()

    def search(self, taxon: GBIFTaxon):

        """
        Looks for a taxon in the GBIF database using name and rank and
        an optional GBIF ID for disambiguation.

        Args:
            taxon: A GBIFTaxon instance

        Returns:
            A GBIFTaxon instance
        """

        if not taxon.is_backbone:
            raise ValueError('Cannot validate non-backbone taxa')

        if taxon.gbif_id is not None:

            # get the record associated with the provided ID
            try:
                id_taxon = self.id_lookup(taxon.gbif_id)
            except GBIFError as err:
                taxon.lookup_status = f'GBIF ID problem: {err.message}'
                return taxon

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
        """Method to return a GBIFTaxon directly from a GBIF ID. It will raise
        a GBIFError if the provided ID cannot be found.

        Params:
            gbif_id: An integer

        Returns:
            A GBIFTaxon object.
        """

        if not isinstance(gbif_id, int):
            raise ValueError('Non-integer GBIF code')

        if not gbif_id >= 0:
            # 0 is kingdom placeholder for incertae sedis
            raise ValueError('Negative GBIF code')

        # get the record associated with the provided ID
        sql = f"select * from backbone where id = {gbif_id}"
        taxon_row = self.gbif_conn.execute(sql).fetchone()

        # check there is a result and that it is congruent with any
        # provided taxon or rank information
        if taxon_row is None:
            raise GBIFError()

        # Create and populate taxon
        taxon = GBIFTaxon(name=taxon_row['canonical_name'],
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
    """This provides a validate method for a GBIFTaxon using the online GBIF
    API. Unlike the LocalGBIFValidator, this doesn't need an __init__ method
    and just contains methods, but duplicates the structure so that
    the two Validators are interchangeable.
    """
    def search(self, taxon: GBIFTaxon):

        """
        Validates a taxon against the GBIF web API. It uses the API endpoint
        species/match?name=XXX&rank=YYY&strict=true
        endpoint to



        # safe to assume that the next least nested taxonomic level is the parent.
        # So, populate GBIF ID using species/match and then use species/{id} to populate the
        # return values

        Args:
            taxon (GBIFTaxon): A GBIFTaxon to be validated

        Returns:
            An updated GBIFTaxon object.
        """

        if not taxon.is_backbone:
            raise ValueError('Cannot validate non-backbone taxa')

        if taxon.gbif_id is not None:

            # get the record associated with the provided ID
            try:
                id_taxon = self.id_lookup(taxon.gbif_id)
            except GBIFError as err:
                taxon.lookup_status = f'GBIF ID problem: {err.message}'
                return taxon

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
        """Method to return a GBIFTaxon directly from a GBIF ID. It will raise
        a GBIFError if the provided ID cannot be found, or if there is a
        connection error.

        Params:
            gbif_id: An integer

        Returns:
            A GBIFTaxon object.
        """

        if not isinstance(gbif_id, int):
            raise ValueError('Non-integer GBIF code')

        if not gbif_id > 0:
            raise ValueError('Negative GBIF code')

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
        taxon = GBIFTaxon(name=response['canonicalName'],
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

class LocalNCBIValidator:
    """This provides a validate method for a NCBITaxon using SQLite queries on a
    locally downloaded copy of the NCBI database.
    """
    def __init__(self, resources):

        conn = sqlite3.connect(resources.ncbi_database)
        conn.row_factory = sqlite3.Row
        self.ncbi_conn = conn

    def __del__(self):

        self.ncbi_conn.close()

    # Functionality to find taxa information from genbank ID
    def id_lookup(self, nnme: str, ncbi_id: int):
        """Method to return full taxonomic information from a NCBI ID. It will
        raise a NCBIError if the provided ID cannot be found, or if there is a
        connection error.

        Params:
            nnme: A nickname to identify the taxon
            ncbi_id: Unique identifer for the taxon

        Returns:
            NCBITaxon: Complete taxon info
        """

        if not isinstance(ncbi_id, int):
            raise TypeError('Non-integer NCBI taxonomy ID')

        if not isinstance(nnme, str):
            raise TypeError('Non-string nickname')

        if not ncbi_id > 0:
            raise ValueError('Negative NCBI taxonomy ID')

        superseed = False

        # Look for record associated with the provided ID
        sql = f"select * from nodes where tax_id = {ncbi_id}"
        taxon_row = self.ncbi_conn.execute(sql).fetchone()

        # If nothing found check if this ID has been merged
        if taxon_row == None:
            sql = f"select * from merge where old_tax_id = {ncbi_id}"
            taxon_row = self.ncbi_conn.execute(sql).fetchone()
            # If it's not found then give bad ID error
            if taxon_row == None:
                raise NCBIError()
            else:
                sql = f"select * from nodes where tax_id = {taxon_row['new_tax_id']}"
                taxon_row = self.ncbi_conn.execute(sql).fetchone()
                superseed = True
                # Warn user that they've given a superseeded taxa ID
                LOGGER.warning(f"NCBI ID {ncbi_id} has been superseeded by ID "
                               f"{taxon_row['tax_id']}")

        # Extract relevant info from the taxon row
        t_rank = taxon_row['rank']
        good_id = taxon_row['tax_id']
        # Then use to find and store name
        sql = f"select * from names where tax_id = {good_id} and "
        sql += "name_class = 'scientific name'"
        name_row = self.ncbi_conn.execute(sql).fetchone()
        t_name = name_row["name_txt"]

        # Then setup loop to find the whole lineage
        lin_fnd = False
        linx = []

        while lin_fnd == False:
            tmp_dic = {}
            # Find node and name of the parent taxon
            sql = f"select * from nodes where tax_id = {taxon_row['parent_tax_id']}"
            taxon_row = self.ncbi_conn.execute(sql).fetchone()
            sql = f"select * from names where tax_id = {taxon_row['tax_id']} and "
            sql += "name_class = 'scientific name'"
            name_row = self.ncbi_conn.execute(sql).fetchone()
            # Store all relevant info
            tmp_dic['TaxID'] = taxon_row['tax_id']
            tmp_dic['ScientificName'] = name_row['name_txt']
            tmp_dic['Rank'] = taxon_row['rank']
            # And add all of it to the Lineage
            linx.append(tmp_dic)
            # End this when the parent taxon is root (ID=1)
            if taxon_row['parent_tax_id'] == 1:
                lin_fnd = True

        # Reverse the order of the lineage
        linx = list(reversed(linx))
        # Find number of taxonomic ranks
        tx_len = len(linx)
        # Extract parent taxa ranks
        rnks = [linx[i]["Rank"] for i in range(tx_len)]

        # Check that the taxon rank provided is a backbone rank
        if t_rank in BACKBONE_RANKS_EX:
            # In this case use provided rank
            rnk = BACKBONE_RANKS_EX.index(t_rank)

        # Filter out ID's with non-backbone ranks (e.g. strains)
        else:
            # Set as not a valid taxa
            vld_tax = False
            # Find lowest index
            r_ID = len(BACKBONE_RANKS_EX) - 1

            # While loop that runs until valid taxa is found
            while vld_tax == False:
                # Check if taxa id is found
                if any([rnks[i] == f"{BACKBONE_RANKS_EX[r_ID]}" for i in range(tx_len)]):
                    # Close loop and store rank number
                    vld_tax = True
                    # Add 1 to the rank as only including lineage in this case
                    rnk = r_ID + 1
                    # Warn user that non-backbone rank has been supplied
                    LOGGER.warning(f'{nnme} of non-backbone rank: {t_rank}')
                # Raise error once backbone ranks have been exhausted
                elif r_ID < 1:
                    LOGGER.error(f'Taxon hierarchy for {nnme} contains no backbone ranks')
                    return
                else:
                    r_ID -= 1

        # Make list of backbone ranks we are looking for
        actual_bb_rnks = BACKBONE_RANKS_EX[0:rnk]

        # Number of missing ranks initialised to zero
        m_rnk = 0


        # Check that all desired backbone ranks appear in the lineage
        if all(item in rnks for item in actual_bb_rnks) == False:
            # Find all missing ranks
            miss = list(set(actual_bb_rnks).difference(rnks))
            # Count missing ranks
            m_rnk = len(miss)
            # Remove missing ranks from our list of desired ranks
            for i in range(0,m_rnk):
                actual_bb_rnks.remove(miss[i])

        # Find valid indices (e.g. those corresponding to backbone ranks)
        vinds = [idx for idx, element in enumerate(rnks) if element in actual_bb_rnks]

        # Create dictonary of valid taxa lineage using a list
        if len(actual_bb_rnks) != 0:
            red_taxa = {f"{actual_bb_rnks[0]}":(str(linx[vinds[0]]["ScientificName"]),
                        int(linx[vinds[0]]["TaxID"]),None)}

            # Requersively add all the hierarchy data in
            for i in range(1,rnk-m_rnk):
                red_taxa[f"{actual_bb_rnks[i]}"] = (str(linx[vinds[i]]["ScientificName"]),
                                                    int(linx[vinds[i]]["TaxID"]),
                                                    int(linx[vinds[i-1]]["TaxID"]))

            # Then add taxa information as a final entry
            red_taxa[f"{t_rank}"] = (str(t_name),
                                     int(good_id),
                                     int(linx[vinds[-1]]["TaxID"]))
        else:
            # Just make taxa with individual rank if no valid lineage provided
            red_taxa = {f"{t_rank}":((t_name),int(good_id),None)}

        # Create and populate microbial taxon
        mtaxon = NCBITaxon(name=t_name, rank=t_rank, ncbi_id=int(good_id),
                           taxa_hier=red_taxa)

        # Record whether a superseeded NCBI ID has been provided
        mtaxon.superseed = superseed

        return mtaxon

    # New function to read in taxa information
    def taxa_search(self, nnme: str, taxah: dict):
        """Method that takes in taxonomic information, and finds the corresponding
        NCBI ID. This NCBI ID is then used to generate a NCBITaxon object,
        which is returned. This function also makes use of parent taxa information
        to distinguish between ambigious taxa names.

        Params:
            nnme: A nickname to identify the taxon
            taxah: A dictonary containing taxonomic information

        Returns:
            NCBITaxon: Complete taxon info
        """

        if isinstance(taxah, dict) == False:
            raise TypeError('Taxa hierarchy should be a dictonary')
        elif all(isinstance(x,str) for x in taxah.keys()) == False:
            raise ValueError('Not all taxa dictionary keys are strings')
        elif all(isinstance(x,str) for x in taxah.values()) == False:
            raise ValueError('Not all taxa dictionary values are strings')

        # Warn user if no backbone ranks have been provided
        if bool(set(taxah.keys()) & set(BACKBONE_RANKS_EX)) == False:
            LOGGER.warning(f"No backbone ranks provided in {nnme}'s taxa hierarchy")

        # Find last dictionary key
        f_key = list(taxah.keys())[-1]

        # Then find corresponding entry as a search term
        s_term = taxah[f_key]

        # Search for name in the local names database
        sql = f"select * from names where name_txt = '%s'" % s_term
        taxon_row = self.ncbi_conn.execute(sql).fetchall()

        # Store count of the number of rows found
        c = len(taxon_row)

        # Check that a singular record has been provided before proceeding
        if c == 1:
            # Find taxa ID as single entry in the list
            taxon_row = self.ncbi_conn.execute(sql).fetchone()
            tID = taxon_row['tax_id']
            # Use ID lookup function to generate as a NCBITaxon object
            mtaxon = self.id_lookup(nnme,tID)

        # Catch cases where no record is found
        elif c == 0:
            # Check if there actually is any higher taxonomy provided
            if len(taxah.keys()) == 1:
                LOGGER.error(f'Taxa {nnme} cannot be found and its higher '
                             f'taxonomic hierarchy is absent')
                return

            # If there is then set up a loop over it
            fnshd = False
            cnt = 1
            while fnshd == False:
                # Increment counter
                cnt += 1
                # Use to find higher taxonomic level to use in search
                f_key = list(taxah.keys())[-cnt]
                new_s_term = taxah[f_key]

                # Search for name in the local names database
                sql = f"select * from names where name_txt = '%s'" % new_s_term
                taxon_row = self.ncbi_conn.execute(sql).fetchall()

                # Process the response
                c = len(taxon_row)

                # Check how many records have been found
                if c == 1:
                    # Find taxa ID as single entry in the list
                    taxon_row = self.ncbi_conn.execute(sql).fetchone()
                    tID = taxon_row['tax_id']
                    # Use ID lookup function to generate as a NCBITaxon object
                    mtaxon = self.id_lookup(nnme,tID)
                    # Store orginally supplied rank
                    mtaxon.orig = list(taxah.keys())[-1]
                    # Warn the user that a higher taxonomic rank is being used
                    LOGGER.warning(f'{s_term} not registered with NCBI, but '
                                   f'higher level taxon {new_s_term} is')
                    fnshd = True
                elif c > 1: # Not going to handle ambiguities in this case
                    LOGGER.error(f'Taxa {nnme} cannot be found and its higher '
                                 f'taxonomic hierarchy is ambigious')
                    return
                # Catch when all the provided hierachy has been exhausted
                elif cnt == len(taxah.keys()):
                    fnshd = True

                # If valid higher taxon not found print error
                if 'mtaxon' not in locals():
                    LOGGER.error(f'Taxa {nnme} cannot be found and neither can '
                                 f'its higher taxonomic hierarchy')
                    return

        # Case where only one rank has been provided
        elif len(taxah) == 1:
            # Preallocate container for the ranks
            t_ranks = []
            # First check if multiple taxa have the same rank
            for i in range(c):
                # Find taxa ID as single entry in the list
                tID = taxon_row[i]['tax_id']
                # Use ID lookup function to generate a tempoary taxon
                temp_taxon = self.id_lookup(nnme,tID)
                # Add rank to list
                t_ranks.append(list(temp_taxon.taxa_hier.keys())[-1])

            # Record whether ranks match expected rank
            if f_key == "kingdom":
                mtch = [True if x == f_key or x == "superkingdom" else False
                        for x in t_ranks]
            else:
                mtch = [True if x == f_key else False for x in t_ranks]

            # If there's only one match, we have determined the correct entry
            if sum(mtch) == 1:
                # Find relevant index
                ind = ([i for i, x in enumerate(mtch) if x])[0]
                # Find taxa ID as single entry in the list
                tID = taxon_row[ind]['tax_id']
                # Use ID lookup function to generate as a NCBITaxon object
                mtaxon = self.id_lookup(nnme,tID)
            # Then check whether multiple taxonomic levels have been provided
            else:
                # If not raise an error
                LOGGER.error(f'Taxa {nnme} cannot be found using only one '
                             f'taxonomic level, more should be provided')
                return
        # Higher ranks provided
        else:
            # Find second from last dictonary key
            f_key = list(taxah.keys())[-2]
            # Then find corresponding entry as a search term
            s_term = taxah[f_key]

            # Search for name in the local names database
            sql = f"select * from names where name_txt = '%s'" % s_term
            p_taxon_row = self.ncbi_conn.execute(sql).fetchall()

            # Store count of the number of records found
            pc = len(p_taxon_row)

            # Check that single parent taxa exists in records
            if pc == 0:
                LOGGER.error(f'Provided parent taxa for {nnme} not found')
                return
            elif pc > 1:
                LOGGER.error(f'More than one possible parent taxa for {nnme} found')
                return
            else:
                # Find parent taxa ID as single entry in the list
                p_taxon_row = self.ncbi_conn.execute(sql).fetchone()
                pID = p_taxon_row['tax_id']
                # Then use ID lookup function to find generate as a NCBITaxon object
                ptaxon = self.id_lookup("parent",pID)

            # Save parent taxa rank and name
            p_key = list(ptaxon.taxa_hier.keys())[-1]
            p_val = ptaxon.taxa_hier[p_key]

            # Use list comprehension to make list of potential taxa
            potents = [self.id_lookup(f"c{i}",taxon_row[i]['tax_id'])
                        for i in range(0,c)]

            # Check if relevant rank exists in the child taxa
            child = [p_key in potents[i].taxa_hier for i in range(0,c)]

            # Then look for matching rank name pairs
            for i in range(0,c):
                # Only check cases where rank exists in the child taxa
                if child[i] == True:
                    child[i] = (potents[i].taxa_hier[f"{p_key}"] == p_val)

            # Check for errors relating to finding too many or few child taxa
            if sum(child) == 0:
                LOGGER.error(f'Parent taxa not actually a valid parent of {nnme}')
                return
            elif sum(child) > 1:
                LOGGER.error(f'Parent taxa for {nnme} refers to multiple '
                             f'possible child taxa')
                return
            else:
                # Find index corresponding to correct child taxa
                tID = int(taxon_row[child.index(True)]['tax_id'])
                # Use ID lookup function to find generate as a NCBITaxon object
                mtaxon = self.id_lookup(nnme,tID)

        # Find last dictonary key
        f_key = list(mtaxon.taxa_hier.keys())[-1]

        # Check if taxonomic rank supplied is used
        if mtaxon.orig == None and f_key != list(taxah.keys())[-1] and f_key != "superkingdom":
            # If not raise an error
            LOGGER.error(f'{list(taxah.values())[-1]} is a {f_key}'
                         f' not a {list(taxah.keys())[-1]}')
            return
        elif list(taxah.keys())[-1] == "kingdom" and f_key == "superkingdom":
            # If not print a warning
            LOGGER.warning(f'NCBI records {(mtaxon.taxa_hier[f_key])[0]} as '
                           f'a superkingdom rather than a kingdom')
        # Then check whether orginally supplied name is still used
        elif mtaxon.orig == None and taxah[f_key] != (mtaxon.taxa_hier[f_key])[0]:
            # If not print a warning
            LOGGER.warning(f'{taxah[f_key]} not accepted usage should be '
                           f'{(mtaxon.taxa_hier[f_key])[0]} instead')
            # And record the fact that usage is superseeded
            mtaxon.superseed = True

        return mtaxon


@enforce_types
class RemoteNCBIValidator:
    """This provides a validate method for a NCBITaxon using the online NCBI Entrez
    tools. This validator duplicates the structure of LocalNCBIValidator so that
    the two Validators are interchangeable.
    """
    # COULD SET THIS UP TO POPULATE FROM THE RESOURCE FILE
    # HOWEVER NCBI SAY THAT THIS SHOULD BE DEVELOPER DEFINED RATHER THAN USER DEFINED
    # (SEE LINK BELOW) "https://www.ncbi.nlm.nih.gov/books/NBK25497/"
    # NEED TO THINK ABOUT TOOL NAME AND EMAIL IF SO
    def __init__(self):

        # Question as to whether I need email and tool here
        # THINK WE NEED TO DEFINE THEM AND REGISTER WITH NCBI, BUT NOT NEEDED BEFORE THEN
        # THINK I NEED TO ADD TOOL AND EMAIL TO THE
        self.email = "jacobcook1995@gmail.com"
        self.tool = "SAFE_data_validator"
        self.api_key = "1738fe86eba2d8fc287ff0d1dcbfeda44a0a"

    def taxonomy_efetch(self, ncbi_id: int):
        """A function that uses the online NCBI eutils function efetch to fetch
        the entry for a specific NCBI ID. This function checks for connection
        errors, and rate limits to ensure that only 10 requests per second are
        made. If this is successful the xml output is stored as an element tree.

        Params:
            ncbi_id: The NCBI ID to fetch the record for

        Returns:
            lxml.etree._Element: Output XML stored as an element tree
        """
        # Construct url
        url = (f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db='
               f'taxonomy&id={ncbi_id}&api_key={self.api_key}')

        # Set up while loop to make the request up to 5 times if neccessary
        success = False
        att = 0
        while success == False and att < 5:
            # Increment counter and make the request
            att += 1
            recrd = requests.get(url)

            # Wait a 10th of a second after each request
            time.sleep(0.1)

            # exit loop if a proper response is recived
            if recrd.status_code == 200:
                success = True

        # raise error if a successful response hasn't been obtaines
        if success == False:
            raise NCBIError('Connection error to remote server')

        # Parse the xml
        root = etree.fromstring(recrd.content)

        # Check to see if the xml contains any information
        if len(root) == 0:
            # If not delete it
            root = None

        return root

    def taxonomy_esearch(self, t_name: str):
        """A function that uses the online NCBI eutils function esearch to search
        for a particular taxon name, and to return information on all matching records.
        This function checks for connection errors, and rate limits to ensure that
        only 10 requests per second are made. If this is successful the xml output
        is stored as an element tree.

        Params:
            t_name: taxon name to search for

        Returns:
            lxml.etree._Element: Output XML stored as an element tree
        """
        # Construct url
        url = (f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db='
               f'taxonomy&term={t_name}&api_key={self.api_key}')

        # Set up while loop to make the request up to 5 times if neccessary
        success = False
        att = 0
        while success == False and att < 5:
            # Increment counter and make the request
            att += 1
            recrd = requests.get(url)

            # Wait a 10th of a second after each request
            time.sleep(0.1)

            # exit loop if a proper response is recived
            if recrd.status_code == 200:
                success = True

        # raise error if a successful response hasn't been obtaines
        if success == False:
            raise NCBIError('Connection error to remote server')

        # Parse the xml
        root = etree.fromstring(recrd.content)

        return root

    # Functionality to find taxa information from genbank ID
    def id_lookup(self, nnme: str, ncbi_id: int):
        """Method to return full taxonomic information from a NCBI ID. It will
        raise a NCBIError if the provided ID cannot be found, or if there is a
        connection error.

        Params:
            nnme: A nickname to identify the taxon
            ncbi_id: Unique identifer for the taxon

        Returns:
            NCBITaxon: Complete taxon info
        """

        if not isinstance(ncbi_id, int):
            raise TypeError('Non-integer NCBI taxonomy ID')

        if not isinstance(nnme, str):
            raise TypeError('Non-string nickname')

        if not ncbi_id > 0:
            raise ValueError('Negative NCBI taxonomy ID')

        # Use efetch to find taxonomy details based on the index
        taxon_row = self.taxonomy_efetch(ncbi_id)

        # Catch case where no entry was found
        if taxon_row == None:
            raise NCBIError()

        # Track down parent ID
        PID = taxon_row.findall("./Taxon/ParentTaxId")

        # Use this to decide whether to track down the lineage
        if PID[0].text != "1":
            # Preallocate list to store lineage
            linx = []
            # Find full lineage
            Lin = taxon_row.findall("./Taxon/LineageEx/Taxon")
            # Then loop over elements in the lineage
            for ind in range (0,len(Lin)):
                tmp_dic = {}
                # Store all relevant info
                tmp_dic['TaxID'] = Lin[ind][0].text
                tmp_dic['ScientificName'] = Lin[ind][1].text
                tmp_dic['Rank'] = Lin[ind][2].text
                # And add all of it to the Lineage
                linx.append(tmp_dic)
        # If no lineage raise an error
        else:
            LOGGER.error(f'Taxon hierarchy for {nnme} contains no backbone ranks')
            return

        # Find number of taxonomic ranks
        tx_len = len(linx)
        # Extract parent taxa ranks
        rnks = [linx[i]["Rank"] for i in range(tx_len)]

        # Find rank of the taxon
        RNK = taxon_row.findall("./Taxon/Rank")

        # Check that the taxon rank provided is a backbone rank
        if RNK[0].text in BACKBONE_RANKS_EX:
            # In this case use provided rank
            rnk = BACKBONE_RANKS_EX.index(RNK[0].text)

        # Filter out ID's with non-backbone ranks (e.g. strains)
        else:
            # Set as not a valid taxa
            vld_tax = False
            # Find lowest index
            r_ID = len(BACKBONE_RANKS_EX) - 1

            # While loop that runs until valid taxa is found
            while vld_tax == False:
                # Check if taxa id is found
                if any([rnks[i] == f"{BACKBONE_RANKS_EX[r_ID]}" for i in range(tx_len)]):
                    # Close loop and store rank number
                    vld_tax = True
                    # Add 1 to the rank as only including lineage in this case
                    rnk = r_ID + 1
                    # Warn user that non-backbone rank has been supplied
                    LOGGER.warning(f'{nnme} of non-backbone rank: {RNK[0].text}')
                # Raise error once backbone ranks have been exhausted
                elif r_ID < 1:
                    LOGGER.error(f'Taxon hierarchy for {nnme} contains no backbone ranks')
                    return
                else:
                    r_ID -= 1

        # Make list of backbone ranks we are looking for
        actual_bb_rnks = BACKBONE_RANKS_EX[0:rnk]

        # Number of missing ranks initialised to zero
        m_rnk = 0

        # Check that all desired backbone ranks appear in the lineage
        if all(item in rnks for item in actual_bb_rnks) == False:
            # Find all missing ranks
            miss = list(set(actual_bb_rnks).difference(rnks))
            # Count missing ranks
            m_rnk = len(miss)
            # Remove missing ranks from our list of desired ranks
            for i in range(0,m_rnk):
                actual_bb_rnks.remove(miss[i])

        # Find valid indices (e.g. those corresponding to backbone ranks)
        vinds = [idx for idx, element in enumerate(rnks) if element in actual_bb_rnks]

        # Find scientific name and ID of the taxon
        SNME = taxon_row.findall("./Taxon/ScientificName")
        TID = taxon_row.findall("./Taxon/TaxId")

        # Create dictonary of valid taxa lineage using a list
        if len(actual_bb_rnks) != 0:
            red_taxa = {f"{actual_bb_rnks[0]}":(str(linx[vinds[0]]["ScientificName"]),
                        int(linx[vinds[0]]["TaxID"]),None)}

            # Requersively add all the hierarchy data in
            for i in range(1,rnk-m_rnk):
                red_taxa[f"{actual_bb_rnks[i]}"] = (str(linx[vinds[i]]["ScientificName"]),
                                                    int(linx[vinds[i]]["TaxID"]),
                                                    int(linx[vinds[i-1]]["TaxID"]))

            # Then add taxa information as a final entry
            red_taxa[f"{RNK[0].text}"] = (str(SNME[0].text),int(TID[0].text),
                                          int(linx[vinds[-1]]["TaxID"]))
        else:
            # Just make taxa with individual rank if no valid lineage provided
            red_taxa = {f"{RNK[0].text}":(str(SNME[0].text),
                        int(TID[0].text),None)}

        # Create and populate microbial taxon
        mtaxon = NCBITaxon(name=SNME[0].text,rank=RNK[0].text,ncbi_id=int(TID[0].text),
                           taxa_hier=red_taxa)

        # Check if AkaTaxIds exists in taxonomic information
        AKA = taxon_row.findall("./Taxon/AkaTaxIds")
        if len(AKA) != 0:
            # Warn user that they've given a superseeded taxa ID
            LOGGER.warning(f"NCBI ID {ncbi_id} has been superseeded by ID "
                            f"{TID[0].text}")
            # Record that a superseeded NCBI ID has been provided
            mtaxon.superseed = True

        return mtaxon

    # New function to read in taxa information
    def taxa_search(self, nnme: str, taxah: dict):
        """Method that takes in taxonomic information, and finds the corresponding
        NCBI ID. This NCBI ID is then used to generate a NCBITaxon object,
        which is returned. This function also makes use of parent taxa information
        to distinguish between ambigious taxa names.

        Params:
            nnme: A nickname to identify the taxon
            taxah: A dictonary containing taxonomic information

        Returns:
            NCBITaxon: Complete taxon info
        """

        if isinstance(taxah, dict) == False:
            raise TypeError('Taxa hierarchy should be a dictonary')
        elif all(isinstance(x,str) for x in taxah.keys()) == False:
            raise ValueError('Not all taxa dictionary keys are strings')
        elif all(isinstance(x,str) for x in taxah.values()) == False:
            raise ValueError('Not all taxa dictionary values are strings')

        # Warn user if no backbone ranks have been provided
        if bool(set(taxah.keys()) & set(BACKBONE_RANKS_EX)) == False:
            LOGGER.warning(f"No backbone ranks provided in {nnme}'s taxa hierarchy")

        # Find last dictonary key
        f_key = list(taxah.keys())[-1]

        # Then find corresponding entry as a search term
        s_term = taxah[f_key]

        # Search the online database
        rcrds = self.taxonomy_esearch(s_term)

        # Track down record count
        CNT = rcrds.findall("./Count")
        c = int(CNT[0].text)

        # Check that a singular record has been provided before proceeding
        if c == 1:
            # Find taxa ID as single entry in the list
            ID = rcrds.findall("./IdList/Id")
            tID = int(ID[0].text)
            # Use ID lookup function to generate as a NCBITaxon object
            mtaxon = self.id_lookup(nnme,tID)
        # Catch cases where no record is found
        elif c == 0:
            # Check if there actually is any higher taxonomy provided
            if len(taxah.keys()) == 1:
                LOGGER.error(f'Taxa {nnme} cannot be found and its higher '
                             f'taxonomic hierarchy is absent')
                return

            # If there is then set up a loop over it
            fnshd = False
            cnt = 1
            while fnshd == False:
                # Increment counter
                cnt += 1
                # Use to find higher taxonomic level to use in search
                f_key = list(taxah.keys())[-cnt]
                new_s_term = taxah[f_key]

                # Search the online database
                rcrds = self.taxonomy_esearch(new_s_term)

                # Track down record count
                CNT = rcrds.findall("./Count")
                c = int(CNT[0].text)

                # Check how many records have been found
                if c == 1:
                    # Find taxa ID as single entry in the list
                    ID = rcrds.findall("./IdList/Id")
                    tID = int(ID[0].text)
                    # Use ID lookup function to generate as a NCBITaxon object
                    mtaxon = self.id_lookup(nnme,tID)
                    # Store orginally supplied rank
                    mtaxon.orig = list(taxah.keys())[-1]
                    # Warn the user that a higher taxonomic rank is being used
                    LOGGER.warning(f'{s_term} not registered with NCBI, but '
                                   f'higher level taxon {new_s_term} is')
                    fnshd = True
                elif c > 1: # Not going to handle ambiguities in this case
                    LOGGER.error(f'Taxa {nnme} cannot be found and its higher '
                                 f'taxonomic hierarchy is ambigious')
                    return
                # Catch when all the provided hierachy has been exhausted
                elif cnt == len(taxah.keys()):
                    fnshd = True

                # If valid higher taxon not found print error
                if 'mtaxon' not in locals():
                    LOGGER.error(f'Taxa {nnme} cannot be found and neither can '
                                 f'its higher taxonomic hierarchy')
                    return

        # Case where only one rank has been provided
        elif len(taxah) == 1:
            # Preallocate container for the ranks
            t_ranks = []
            # Find all the taxa ID's
            ID = rcrds.findall("./IdList/Id")
            # First check if multiple taxa have the same rank
            for i in range(c):
                # Find taxa ID as specific entry in the list
                tID = int(ID[i].text)
                # Use ID lookup function to generate a tempoary taxon
                temp_taxon = self.id_lookup(nnme,tID)
                # Add rank to list
                t_ranks.append(list(temp_taxon.taxa_hier.keys())[-1])

            # Record whether ranks match expected rank
            if f_key == "kingdom":
                mtch = [True if x == f_key or x == "superkingdom" else False
                        for x in t_ranks]
            else:
                mtch = [True if x == f_key else False for x in t_ranks]

            # If there's only one match, we have determined the correct entry
            if sum(mtch) == 1:
                # Find relevant index
                ind = ([i for i, x in enumerate(mtch) if x])[0]
                # Find taxa ID as single entry in the list
                tID = int(ID[ind].text)
                # Use ID lookup function to generate as a NCBITaxon object
                mtaxon = self.id_lookup(nnme,tID)
            # Then check whether multiple taxonomic levels have been provided
            else:
                # If not raise an error
                LOGGER.error(f'Taxa {nnme} cannot be found using only one '
                             f'taxonomic level, more should be provided')
                return
        # Higher ranks provided
        else:
            # Find second from last dictonary key
            f_key = list(taxah.keys())[-2]
            # Then find corresponding entry as a search term
            s_term = taxah[f_key]

            # Search the online database for parent record
            p_rcrds = self.taxonomy_esearch(s_term)

            # Track down record count
            PCNT = p_rcrds.findall("./Count")
            pc = int(PCNT[0].text)

            # Check that single parent taxa exists in records
            if pc == 0:
                LOGGER.error(f'Provided parent taxa for {nnme} not found')
                return
            elif pc > 1:
                LOGGER.error(f'More than one possible parent taxa for {nnme} found')
                return
            else:
                # Find parent taxa ID as single entry in the list
                PAR = p_rcrds.findall("./IdList/Id")
                pID = int(PAR[0].text)
                # Then use ID lookup function to find generate as a NCBITaxon object
                ptaxon = self.id_lookup("parent",pID)

            # Save parent taxa rank and name
            p_key = list(ptaxon.taxa_hier.keys())[-1]
            p_val = ptaxon.taxa_hier[p_key]

            # Find ID's of potential child taxa
            ID = rcrds.findall("./IdList/Id")

            # Use list comprehension to make list of potential taxa
            potents = [self.id_lookup(f"c{i}",int(ID[i].text))
                        for i in range(0,c)]

            # Check if relevant rank exists in the child taxa
            child = [p_key in potents[i].taxa_hier for i in range(0,c)]

            # Then look for matching rank name pairs
            for i in range(0,c):
                # Only check cases where rank exists in the child taxa
                if child[i] == True:
                    child[i] = (potents[i].taxa_hier[f"{p_key}"] == p_val)

            # Check for errors relating to finding too many or few child taxa
            if sum(child) == 0:
                LOGGER.error(f'Parent taxa not actually a valid parent of {nnme}')
                return
            elif sum(child) > 1:
                LOGGER.error(f'Parent taxa for {nnme} refers to multiple '
                             f'possible child taxa')
                return
            else:
                # Find index corresponding to correct child taxa
                tID = int(ID[child.index(True)].text)
                # Use ID lookup function to find generate as a NCBITaxon object
                mtaxon = self.id_lookup(nnme,tID)

        # Find last dictonary key
        f_key = list(mtaxon.taxa_hier.keys())[-1]

        # Check if taxonomic rank supplied is used
        if mtaxon.orig == None and f_key != list(taxah.keys())[-1] and f_key != "superkingdom":
            # If not raise an error
            LOGGER.error(f'{list(taxah.values())[-1]} is a {f_key}'
                         f' not a {list(taxah.keys())[-1]}')
            return
        elif list(taxah.keys())[-1] == "kingdom" and f_key == "superkingdom":
            # If not print a warning
            LOGGER.warning(f'NCBI records {(mtaxon.taxa_hier[f_key])[0]} as '
                           f'a superkingdom rather than a kingdom')
        # Then check whether orginally supplied name is still used
        elif mtaxon.orig == None and taxah[f_key] != (mtaxon.taxa_hier[f_key])[0]:
            # If not print a warning
            LOGGER.warning(f'{taxah[f_key]} not accepted usage should be '
                           f'{(mtaxon.taxa_hier[f_key])[0]} instead')
            # And record the fact that usage is superseeded
            mtaxon.superseed = True

        return mtaxon



class GBIFTaxa:

    def __init__(self, resources: Resources):
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

    @loggerinfo_push_pop('Loading GBIFTaxa worksheet')
    def load(self, worksheet):
        """Loads a set of taxa from the rows of a SAFE formatted GBIFTaxa worksheet
        and then adds the higher taxa for those rows.

        Args:
            worksheet: An openpyxl worksheet instance following the GBIFTaxa formatting

        Returns:
            Updates the taxon_names and taxon_index attributes of the class instance
            using the data in the worksheet.
        """

        start_errors = COUNTER_HANDLER.counters['ERROR']

        # Get the data read in.
        LOGGER.info("Reading taxa data")
        FORMATTER.push()
        dframe = GetDataFrame(worksheet)

        if not dframe.data_columns:
            LOGGER.error('No data or only headers in GBIFTaxa worksheet')
            return

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

            # Standardise blank values to None
            row = {ky: None if blank_value(vl) else vl for ky, vl in row.items()}
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
        self.n_errors = COUNTER_HANDLER.counters['ERROR'] - start_errors
        if self.n_errors > 0:
            LOGGER.info('GBIFTaxa contains {} errors'.format(self.n_errors))
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
        validates it and updates the GBIFTaxa instance to include the new details.
        This is principally used to process rows found in a GBIFTaxa worksheet,
        but is deliberately separated out so that a GBIFTaxa instance can be populated
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
            taxon_input: GBIFTaxon information in standard form as above

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
            p_taxon = GBIFTaxon(name=parent_info[0], rank=parent_info[1], gbif_id=parent_info[2])

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
        # else:
        #         LOGGER.info('No parent taxon provided')

        # Now check main taxa
        #
        # The parent list is now populated with parent GBIFTaxon objects keyed by
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
        m_taxon = GBIFTaxon(name=taxon_info[0], rank=taxon_info[1], gbif_id=taxon_info[2])
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

class NCBITaxa:

    def __init__(self, resources: Resources):
        """A class to hold a list of taxon names and a validated taxonomic
        index for those taxa and their taxonomic hierarchy. The validate_taxon
        method checks that the provided taxon hierarchy and (optional) NCBI ID can be
        validated against the NCBI database (and both refer to the same taxon).
        This then populates two things:

        i)  the taxon_names attribute of the dataset, which is just a set of
            names used as a validation list for taxon names used in data worksheets.
        ii) the taxon_index attribute of the dataset, which contains a set
            of lists structured as:

                [worksheet_name (str),
                 ncbi_id (int),
                 ncbi_parent_id (int),
                 canonical_name (str),
                 taxonomic_rank (str),
                 taxon_superseeded (bool)]

        We only populate orginal_taxon_rank when the initially supplied taxon is
        not found in NCBI, but a parent (or grandparent etc) is found. In this case
        the details of the higher taxon are stored along with the orginal worksheet_name
        and with the orginal rank stored in orginal_taxon_rank. In all other cases
        orginal_taxon_rank is set to None.

        Where a taxon is superseeded in the NCBI database, two entries are
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
        self.hierarchy = set()
        self.n_errors = None

        # Get a validator instance
        if resources.use_local_ncbi:
            self.validator = LocalNCBIValidator(resources)
        else:
            self.validator = RemoteNCBIValidator()

    @loggerinfo_push_pop('Loading NCBITaxa worksheet')
    def load(self, worksheet):
        """Loads a set of taxa from the rows of a SAFE formatted NCBITaxa worksheet
        and then adds the higher taxa for those rows.

        Args:
            worksheet: An openpyxl worksheet instance following the NCBITaxa formatting

        Returns:
            None: Updates the taxon_names and taxon_index attributes of the class instance
            using the data in the worksheet.
        """

        start_errors = COUNTER_HANDLER.counters['ERROR']

        # Get the data read in.
        LOGGER.info("Reading NCBI taxa data")
        FORMATTER.push()
        dframe = GetDataFrame(worksheet)

        if not dframe.data_columns:
            LOGGER.error('No data or only headers in Taxa worksheet')
            return

        # Dupe headers likely cause serious issues, so stop
        if 'duplicated' in dframe.bad_headers:
            LOGGER.error('Cannot parse taxa with duplicated headers')
            return

        # Get the headers
        headers = IsLower(dframe.headers).values

        # Field cleaning
        core_fields = {'name', 'ncbi id'}
        missing_core = core_fields.difference(headers)

        if missing_core:
            # core names are not found so can't continue
            LOGGER.error('Missing core fields: ', extra={'join': missing_core})
            return

        # Check that at least two backbone taxa have been provided
        fnd_rnks = set(BACKBONE_RANKS_EX).intersection(headers)

        if len(fnd_rnks) < 2:
            # can't continue if less than two backbone ranks are provided
            LOGGER.error('Less than two backbone taxonomic ranks are provided')
            return

        # Find core field indices and use to isolate non core (i.e. taxonomic) headers
        core_inds = [headers.index(item) for item in core_fields]
        non_core = [element for i, element in enumerate(headers) if i not in core_inds]

        # Check to see if comments is provided
        if "comments" in non_core:
            # Check that it is the last header
            if non_core[-1] != "comments":
                LOGGER.error(f"If 'Comments' is provided as a field it must be "
                             f"the last column")
                return
            else:
                # If it's the last header go ahead and delete it
                del non_core[-1]

        # Check that backbone ranks are in the correct order
        l1 = [v for v in BACKBONE_RANKS_EX if v in fnd_rnks]
        l2 = [v for v in non_core if v in fnd_rnks]

        if l1 != l2:
            LOGGER.error('Backbone taxonomic ranks not provided in the correct order')

        # Check that if subspecies has been provided species is provided
        if 'subspecies' in headers:
            if 'species' not in headers:
                LOGGER.error('If subspecies is provided so must species')
                return

        # Same check
        if 'species' in headers:
            if 'genus' not in headers:
                LOGGER.error('If species is provided so must genus')
                return

        # Any duplication in names
        dupl_taxon_names = HasDuplicates(dframe.data_columns[headers.index('name')])

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

        # Change data format such that it is appropriate to be used to search
        # the NCBI database, i.e. convert from k__ notation, generate species
        # binomials and subspecies trinominals, this is then all collected in a
        # dictonary.

        for idx, row in enumerate(taxa):

            # Strip k__ notation so that the text is appropriate for the search
            for rnk in non_core:
                row[rnk], match = taxa_strip(row[rnk],rnk)
                if match == False:
                    print(f"Implied rank of {row[rnk]} in row {idx + 1} does not"
                          f" match rank it is assigned")

            # Standardise blank values to None
            row = {ky: None if blank_value(vl) else vl for ky, vl in row.items()}
            # Replace any NA values with None
            row = {ky: None if vl == "NA" else vl for ky, vl in row.items()}

            # Start with empty dictonary for taxonomic hierarchy
            taxa_hier = {}

            # Loop over all ranks to populate the dictonary
            for rnk in non_core:
                if row[rnk] != None:
                    if rnk == "species":
                        taxa_hier[rnk] = species_binomial(row["genus"], row[rnk])
                    elif rnk == "subspecies":
                        taxa_hier[rnk] = subspecies_trinomial(taxa_hier["species"], row[rnk])
                    else:
                        taxa_hier[rnk] = row[rnk]

            self.taxon_names.update([row['name']])
            LOGGER.info(f"Validating row {idx + 1}: {row['name']}")
            FORMATTER.push()
            self.validate_and_add_taxon((row['name'], taxa_hier, row['ncbi id']))
            FORMATTER.pop()

        # Add the higher taxa
        self.index_higher_taxa()

        # summary of processing
        self.n_errors = COUNTER_HANDLER.counters['ERROR'] - start_errors
        if self.n_errors > 0:
            LOGGER.info('Taxa contains {} errors'.format(self.n_errors))
        else:
            LOGGER.info('{} taxa loaded correctly'.format(len(self.taxon_names)))

        FORMATTER.pop()


    def validate_and_add_taxon(self, ncbi_taxon_input: list):
        """ User information is provided that names a taxon and (optionally) gives
        a NCBI taxonomy ID. This information is then validated against the NCBI database.
        This is principally used to process rows found in a NCBITaxa worksheet, but
        is deliberately separated out so that a NCBITaxa instance can be populated
        independently of an Excel dataset.

        The ncbi_taxon_input has the form:

        ['m_name', 'taxon_hier', 'ncbi_id']

        If there is no NCBI ID, the structure is:
        ['m_name', 'taxon_hier', None]

        Args:
            ncbi_taxon_input: NCBITaxon information in standard form as above

        Returns:
            None: Updates the taxon_names and taxon_index attributes of the class instance.
        """

        m_name, taxon_hier, ncbi_id = ncbi_taxon_input

        # Sanitise worksheet names for taxa - only keep unpadded strings.
        if m_name is None or not isinstance(m_name, str) or m_name.isspace() or not m_name:
            LOGGER.error('Worksheet name missing, whitespace only or not text')
            return
        elif m_name != m_name.strip():
            LOGGER.error(f"Worksheet name has whitespace padding: {repr(m_name)}")
            m_name = m_name.strip()
            self.taxon_names.add(m_name)
        else:
            self.taxon_names.add(m_name)

        # Check that ID provided is reasonable (if provided)
        i_fail = False

        # ID can be None or an integer (openpyxl loads all values as float)
        if not(ncbi_id is None or
               (isinstance(ncbi_id, float) and ncbi_id.is_integer()) or
                isinstance(ncbi_id, int)) :
            LOGGER.error('NCBI ID contains value that is not an integer')
            i_fail = True

        # Check the main taxon details
        h_fail = False

        # Check that a dictonary with at least one entry has been provided
        if not isinstance(taxon_hier,dict) or len(taxon_hier) == 0:
            LOGGER.error('Taxa hierarchy should be a (not empty) dictonary')
            h_fail = True
        # Otherwise check for padding of dictonary keys and values
        else:
            # Make a translation table
            translate = {}
            # Loop over all dictonary keys
            for idx in taxon_hier.keys():
                if not isinstance(idx,str):
                    LOGGER.error(f"Non-string dictonary key used: {repr(idx)}")
                    h_fail = True
                elif idx == "" or idx.isspace():
                    LOGGER.error('Empty dictonary key used')
                    h_fail = True
                elif idx != idx.strip():
                    LOGGER.error(f"Dictonary key has whitespace padding: {repr(idx)}")
                    # Save keys to swap to new translation table
                    translate[idx] = idx.strip()

                # Extract corresponding dictonary value
                val = taxon_hier[idx]
                # Then perform similar checks on dictonary value
                if not isinstance(val,str):
                    LOGGER.error(f"Non-string dictonary value used: {repr(val)}")
                    h_fail = True
                elif val == "" or val.isspace():
                    LOGGER.error('Empty dictonary value used')
                    h_fail = True
                elif val != val.strip():
                    LOGGER.error(f"Dictonary value has whitespace padding: {repr(val)}")
                    taxon_hier[idx] = val.strip()

            # Use translation table to replace whitespaced dictonary keys
            for old, new in translate.items():
                taxon_hier[new] = taxon_hier.pop(old)

        # Now check that the taxa hierarchy is correctly ordered
        if i_fail:
            LOGGER.error(f'Improper NCBI ID provided, cannot be validated')

        if h_fail:
            LOGGER.error(f'Taxon details not properly formatted, cannot validate')

        if h_fail or i_fail:
            return

        # Now check that dictonary containing taxa hierarchy is properly ordered
        if len(taxon_hier) > 1: # Only matters if it contains multiple entries
            # Find all keys in taxa hierarchy
            t_ord = list(taxon_hier.keys())
            # Find corresponding keys from backbone
            b_ord = [x for x in BACKBONE_RANKS_EX if x in t_ord]
            # Remove non-backbone keys from taxa order
            t_ord = [x for x in t_ord if x in b_ord]
            # Then catch cases where orders don't match
            if b_ord != t_ord:
                LOGGER.error(f'Taxon hierarchy not in correct order')
                return

        # Go straight ahead and search for the taxon
        hr_taxon = self.validator.taxa_search(m_name, taxon_hier)
        # Catch case where errors are returned rather than a taxon
        if hr_taxon == None:
            LOGGER.error(f'Search based on taxon hierarchy failed')
            return

        # Then check if a genbank ID number has been provided
        if ncbi_id != None:
            id_taxon = self.validator.id_lookup(m_name, int(ncbi_id))
            # Check if taxonomy hierarchy superseeded
            if hr_taxon.superseed == True:
                LOGGER.warning(f'Taxonomic classification superseeded for '
                f'{m_name}, using new taxonomic classification')
            elif id_taxon.superseed == True:
                LOGGER.warning(f'NCBI taxa ID superseeded for {m_name}'
                f', using new taxa ID')
            elif id_taxon != hr_taxon:
                LOGGER.error(f"The NCBI ID supplied for {m_name} does "
                f"not match hierarchy: expected {hr_taxon.ncbi_id}"
                f" got {ncbi_id}")
                return
        else:
            # Warn user that superseeded taxonomy won't be used
            if hr_taxon.superseed:
                LOGGER.warning(f'Taxonomic classification superseeded for '
                f'{m_name}, using new taxonomic classification')

        # Check that the hierachy found matches
        match = self.compare_hier(m_name,hr_taxon,taxon_hier)

        # Find parent ID
        if len(hr_taxon.taxa_hier) > 1:
            # Find taxon one level up
            f_key = list(hr_taxon.taxa_hier.keys())[-2]
            parent_id = (hr_taxon.taxa_hier[f_key])[1]
        else:
            # Set to None if hierachy is empty (i.e. top level taxa)
            parent_id = None

        # Catch cases where orginal taxon has not been found
        if hr_taxon.orig != None:
            self.taxon_index.append([m_name, -1, hr_taxon.ncbi_id,
                                     list(taxon_hier.values())[-1], hr_taxon.orig,
                                     False])
        else:
            # Then check if taxon is superseeded
            if hr_taxon.superseed == True or (ncbi_id != None and id_taxon.superseed == True):
                if ncbi_id == None or id_taxon.superseed == False:
                    superseed_id = hr_taxon.ncbi_id
                    # Find supplied name using lowest found rank
                    f_key = list(hr_taxon.taxa_hier.keys())[-1]
                    superseed_name = taxon_hier[f_key]
                else:
                    # Supplied ID is superseeded
                    superseed_id = ncbi_id
                    # Check if supplied name superseeded
                    if hr_taxon.superseed == True:
                        f_key = list(hr_taxon.taxa_hier.keys())[-1]
                        superseed_name = taxon_hier[f_key]
                    else:
                        superseed_name = hr_taxon.name

                # Add superseeded taxon to the index
                self.taxon_index.append([m_name, superseed_id, parent_id,
                                            superseed_name, hr_taxon.rank,
                                            True])

            # Then (also) add non-superseeded taxon info to the index
            self.taxon_index.append([m_name, hr_taxon.ncbi_id, parent_id,
                                        hr_taxon.name, hr_taxon.rank,
                                        False])

        self.hierarchy.update(list(hr_taxon.taxa_hier.items()))

        # Check if this has succeded without warnings or errors
        if hr_taxon.superseed == False and hr_taxon.orig == None:
            # Straight forward when there isn't a genbank id, or previously processed
            if ncbi_id == None:
                # If so inform the user of this
                LOGGER.info(f'Taxon ({m_name}) found in NCBI database')
            # Otherwise need to check for superseeded ID's
            elif id_taxon.superseed == False:
                LOGGER.info(f'Taxon ({m_name}) found in NCBI database')
        elif hr_taxon.superseed == False:
            LOGGER.info(f'Higher taxon for ({m_name}) resolved in NCBI')

    @loggerinfo_push_pop('Indexing taxonomic hierarchy')
    def index_higher_taxa(self):
        """ Function to add higher taxa to the index. The taxa are added are those
        stored in self.hierarchy that have not already been added to the taxon index.
        """

        # Use the taxon hierarchy entries to add higher taxa
        # - drop taxa with a GBIF ID already in the index

        # Only add if taxon is not already included in known
        known = [tx[1] for tx in self.taxon_index]
        to_add = [tx for tx in self.hierarchy if tx[1][1] not in known]
        to_add.sort(key=lambda val: BACKBONE_RANKS_EX.index(val[0]))

        # Look up the taxonomic hierarchy
        for tx_lev, (tx_nme, tx_id, p_id) in to_add:
            # Add all this to the taxonomy
            self.taxon_index.append([None, tx_id, p_id, tx_nme, tx_lev,
                                     False])
            LOGGER.info(f'Added {tx_lev} {tx_nme}')

    def compare_hier(self, m_name: str, mtaxon: NCBITaxon, taxon_hier: dict):
        """ Function to compare the hierachy of a taxon with the hierachy that was
        initially supplied. This function only checks that provided information
        matches, missing levels or entries are not checked for"""
        # Not worth checking hierachy in superseeded case
        if mtaxon.superseed == True:
            return

        # Find all common ranks between the hierarchies
        rnks = list(set(taxon_hier.keys()) & set(mtaxon.taxa_hier.keys()))

        # Find if all taxa names match
        t_match = [taxon_hier[r] == (mtaxon.taxa_hier[r])[0] for r in rnks]
        if all(t_match) == False:
            # Find indices of non-matching values
            inds = [i for i in range(len(t_match)) if t_match[i] == False]
            # Then loop over these indices giving appropriate warnings
            for ind in inds:
                LOGGER.warning(f'Hierarchy mismatch for {m_name} its {rnks[ind]}'
                               f' should be {(mtaxon.taxa_hier[rnks[ind]])[0]} not '
                               f'{taxon_hier[rnks[ind]]}')

    @property
    def is_empty(self):
        return len(self.taxon_names) == 0

class Taxa:

    def __init__(self, resources: Resources):
        """A class that combines information stored in seperate lower level GBIFTaxa
        and NCBITaxa worksheets. We are interested in checking that no worksheet
        names are reused when both Taxa sheets are provided, that every worksheet
        name is used somewhere in the Data worksheets, and that every taxon name
        used across the Data worksheets is defined in a Taxa worksheet.

        THIS ISN'T ACTUALLY WHAT IS DONE, UPDATE THIS
        To perform this check three lists are populated: `repeat_names` which records
        any worksheet name used in both Taxa worksheets, `unused_names` which records
        any names defined in the Taxa worksheets which are not used in any Data
        worksheet and `taxon_names_used` which records all names used in the Data
        worksheets.

        In addition, GBIFTaxa and NCBITaxa instances are generated and stored within
        this overarching class.
        """

        self.gbif_taxa = GBIFTaxa(resources)
        self.ncbi_taxa = NCBITaxa(resources)
        self.taxon_names_used = set()

    @property
    def is_empty(self):
        return (self.gbif_taxa.is_empty and self.ncbi_taxa.is_empty)

    @property
    def taxon_names(self):
        return self.gbif_taxa.taxon_names.union(self.ncbi_taxa.taxon_names)

def taxon_index_to_text(taxon_index, html=False, indent_width=4):
    """
    Turns the taxon index from a GBIFTaxa instance into a text representation
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

# New function to remove k__ notation, and to check that rank and
def taxa_strip(name: str, rank: str):
    """A function to remove k__ type notation from taxa names. The function also
    checks that the provided rank is consistent with the rank implied by its prefix.
    """
    if name == None:
        return (None, True)
    elif "__" in name:
        # Strip name down
        ind = name.rfind("_")
        s_name = name[ind+1:]
        # Check if ranks match
        match = (name[0].lower() == rank[0].lower())
        return (s_name, match)
    else:
        return (name, True)

# New function to generate species binomial
def species_binomial(genus: str, species: str):
    """A function to construct a species binomal from a genus name and a species
    name.
    """
    # First check if species is a single name
    if len(species.split()) == 1 and len(genus.split()) == 1:
        return genus.strip() + " " + species.strip()
    # Look for Candidatus formats
    elif "candidatus" in species.lower() or "candidatus" in genus.lower():
        if "candidatus" in species.lower():
            # Construct binomal with first word of species name removed
            bi = genus.strip()
            for i in range(1,len(species.split())):
                bi = bi + " " + species.split()[i]
            return bi
        else:
            return genus.strip() + " " + species.strip()
    # Then check that species name is more words than the genus name
    elif len(species.split()) > len(genus.split()):
        if genus in species:
            return species
        else:
            LOGGER.error(f'Species name ({species}) appears to be binomal but '
                         f'does not contain genus name ({genus})')
            return None
    else:
        LOGGER.error(f'Genus name ({genus}) appears to be too long')
        return None

# Equivalent function to generate subspecies trinominal
def subspecies_trinomial(species: str, subspecies: str):
    """A function to construct a subspecies trinomal from a species name and a
    subspecies name.
    """
    # First check if subspecies is a single name
    if len(subspecies.split()) == 1 and len(species.split()) == 2:
        return species.strip() + " " + subspecies.strip()
    # Look for Candidatus formats
    elif "candidatus" in subspecies.lower() or "candidatus" in species.lower():
        if "candidatus" in subspecies.lower():
            # Construct trinominal with first word of subspecies name removed
            tri = species.strip()
            for i in range(1,len(subspecies.split())):
                tri = tri + " " + subspecies.split()[i]
            return tri
        else:
            return species.strip() + " " + subspecies.strip()
    # Catch too short species name case
    elif len(species.split()) == 1:
        LOGGER.error(f'Species name ({species}) too short')
        return None
    # Then check that subspecies name is more words than the species name
    elif len(subspecies.split()) > len(species.split()):
        if species in subspecies:
            return subspecies
        else:
            LOGGER.error(f'Subspecies name ({subspecies}) appears to be trinomal'
                         f'but does not contain species name ({species})')
            return None
    else:
        LOGGER.error(f'Species name ({species}) too long')
        return None

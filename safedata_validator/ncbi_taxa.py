"""## The ncbi_taxa submodule
This module describes classes and methods used to compile taxonomic data from
datasets and to validate taxonomy against the GBIF backbone database.

The NCBITaxon dataclass is used to store data about a taxon entry in a dataset. It
is initialised with user data and then the taxon Validator classes can be used
to update a NCBITaxon object with the result of NCBI validation. The two validation
classes use either the set of online NCBI Entrez Eutils (`RemoteNCBIValidator`) or
faster validation against a local copy of the NCBI taxonomy database (`LocalNCBIValidator`).

The dataset 'NCBITaxa' worksheet provides a set of taxonomic entries and the NCBITaxa
class is used to load and collate the set of taxonomic entries from a dataset.

Supplied taxa of any rank (i.e. strain or clade) which can be successfully validated
will be recorded. However, associated higher taxa will only be recorded if their ranks
are either a GBIF backbone rank or superkingdom.
"""

from typing import Union, Optional
import dataclasses
import sqlite3
import time
from lxml import etree

import urllib
import requests
from enforce_typing import enforce_types

from safedata_validator.resources import Resources
from safedata_validator.logger import (COUNTER_HANDLER, FORMATTER, LOGGER,
                                       loggerinfo_push_pop)
from safedata_validator.validators import (GetDataFrame, HasDuplicates,
                                           IsLower, blank_value)

# TODO - Modify the resource file to ask the user to provide an email address
# This should only be done if the user actually wants to use this module as it
# isn't need elsewhere (as far as I know)
# Also should ask for api key

# TODO - Work out how this best integrates with taxa
# Lot of different potential approaches (e.g. making making NCBI and GBIF validators
# sub-classes of a validator class). Need to establish what the most extensible
# and stable option is

# Extended version of backbone ranks to capture superkingdoms
BACKBONE_RANKS_EX = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                    'family', 'genus', 'species', 'subspecies']

class NCBIError(Exception):
    """Exception class for remote NCBI errors

    Attributes:
        message: explanation of the error
    """

    def __init__(self, message="No entry found for ID, probably bad ID"):
        self.message = message
        super().__init__(self.message)

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

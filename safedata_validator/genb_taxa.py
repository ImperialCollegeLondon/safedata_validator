"""## The genb_taxa submodule
This module describes classes and methods used to validate taxonomic data against
the NCBI database and convert it to a form compatible with the GBIF backbone database.

The NCBITaxon dataclass is used to store data about a taxon entry in NCBI (GenBank)
form. It is initialised with user data and then the taxon Validator class can
be used to update a NCBITaxon object with the result of NCBI validation. At present
only online validation is possible through the `RemoteNCBIValidator` class.

FURTHER DETAILS OF STUFF DEFINED.

THIS IS STILL VERY MUCH A WORK IN PROGRESS SO MORE INFORMATION IS GOING TO BE
ENTERED HERE
"""

from typing import Union, Optional
import dataclasses
from Bio import Entrez

import requests
import urllib
from enforce_typing import enforce_types

from safedata_validator.logger import (COUNTER_HANDLER, FORMATTER, LOGGER,
                                       loggerinfo_push_pop)

# ADD TO RESOURCE FILE, AS A CHECK THAT USER HAS PROVIDED ONE
# Key should also be added to the resource file
Entrez.email = "jacobcook1995@gmail.com"
# Hard coding api key in for now
user_key = "1738fe86eba2d8fc287ff0d1dcbfeda44a0a"

# Extended version of backbone ranks to capture superkingdoms
# WILL NEED THOUGHT ABOUT THE MAPPING
BACKBONE_RANKS_EX = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                    'family', 'genus', 'species', 'subspecies']

class NCBIError(Exception):
    """Exception class for remote NCBI errors

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message="NCBI taxa ID not found"):
        self.message = message
        super().__init__(self.message)


# TODO - Correctly read in xslx data
# Lot of this has already been written by David, main thing is that I will have
# to decide on the formatting, like what data are we taking in.

# TODO - Think about data output
# Have to make sure that the indexing is compatibale with David's database
# Probably have to add in the first databasing steps as well

# POTENTIAL - Take steps to increase the speed
# Make local copy by downloading the relevant part of NCBI's database and
# running the validation locally

# TODO - Modify the resource file to ask the user to provide an email address
# This should only be done if the user actually wants to use this module as it
# isn't need elsewhere (as far as I know)
# Also should ask for api key

# QUESTIONS FOR DAVID
# SHOULD A YOU ARE NOT CONNCTED TO THE INTERNET ERROR BE SETUP?
# IS THERE A SENSIBLE WAY TO DUMMY CONNECTION ERRORS?

@enforce_types
@dataclasses.dataclass
class NCBITaxon:
    """Holds taxonomic information from a user on a microbial taxa. This can be
    populated using NCBI GenBank validation.

    There are 3 class properties that can be used to create an instance:
        * name
        * taxa_hier: Dictionary of valid taxonomic hierarchy (with ID's)
        * genbank_id: GenBank ID for full taxa (i.e. including non-backbone ranks)
    The remaining properties are populated by processing functions not when an
    instance is created.
        * diverg: does the GenBank ID and stored taxa info diverge, and if so how?
        * superseed: is supplied taxon name/ID still the accepted usage
    """

    # Init properties
    name: str
    taxa_hier: dict
    genbank_id: Optional[Union[int, float]] = None
    diverg: str = dataclasses.field(init=False)
    superseed: bool = dataclasses.field(init=False)

    def __post_init__(self):
        """Sets the defaults for the post-init properties and checks inputs
        """

        if self.genbank_id is not None:
            if isinstance(self.genbank_id, float) and not self.genbank_id.is_integer():
                raise TypeError('GenBank Id is not an integer')
            self.genbank_id = int(self.genbank_id)

        if self.taxa_hier is not None:
            if len(self.taxa_hier) == 0:
                raise ValueError('Taxa hierarchy dictonary empty')
            elif all(isinstance(x,str) for x in self.taxa_hier.keys()) == False:
                raise ValueError('Not all taxa dictionary keys are strings')
            elif all(isinstance(x,tuple) for x in self.taxa_hier.values()) == False:
                raise ValueError('Not all taxa dictionary values are tuples')
            elif all(list(map(type,x)) == [str, int] for x in
                 self.taxa_hier.values()) == False:
                 raise ValueError('Taxa tuples not all in [string integer] form')

        self.diverg = None
        self.superseed = False

@enforce_types
class RemoteNCBIValidator:
    """This provides a validate method for a MicrobeTaxon using various online
    NCBI APIs. This doesn't need an __init__ method and just contains methods.
    """
    # Functionality to find taxa information from genbank ID
    def id_lookup(self, nnme: str, genbank_id: int):
        """Method to return full taxonomic information from a GenBank ID. This
        includes details of any potential synonymus names.

        Params:
            nnme: A nickname to identify the taxon
            genbank_id: An integer

        Returns:
            A NCBITaxon object
        """

        if not isinstance(genbank_id, int):
            raise ValueError('Non-integer NCBI taxonomy ID')

        if not genbank_id > 0:
            raise ValueError('Negative NCBI taxonomy ID')

        # Set status of taxa validity as initally false
        valid = False

        # Use efetch to find taxonomy details based on the index
        try:
            handle = Entrez.efetch(db="taxonomy",id=f"{genbank_id}",
                                   retmode="xml",api_key=user_key)
            tax_dic = Entrez.read(handle)[0]
            handle.close()
        # And case where a connection can't be made to the remote server
        except urllib.error.HTTPError:
            raise NCBIError('Connection error to remote server')
        # Catch case where the index isn't found
        except IndexError:
            raise NCBIError("No entry found from ID, probably bad ID")

        # Then extract the lineage
        linx = tax_dic["LineageEx"]
        # Find number of taxonomic ranks
        tx_len = len(linx)
        # Extract parent taxa ranks
        rnks = [linx[i]["Rank"] for i in range(tx_len)]

        # Check that the taxon rank provided is a backbone rank
        if tax_dic["Rank"] in BACKBONE_RANKS_EX:
            # Set as a valid taxa
            valid = True
            # In this case use provided rank
            rnk = BACKBONE_RANKS_EX.index(tax_dic["Rank"])

        # Filter out ID's without ranks (e.g. strains)
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
                    LOGGER.warning(f'{nnme} not of backbone rank, instead resolved '
                                   f'to {BACKBONE_RANKS_EX[rnk-1]} level')
                # Raise error once backbone ranks have been exhausted
                elif r_ID < 1:
                    LOGGER.error(f'Taxon rank of {nnme} cannot be found in backbone')
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
        red_taxa = {f"{actual_bb_rnks[i]}":(str(linx[vinds[i]]["ScientificName"]),
                    int(linx[vinds[i]]["TaxId"])) for i in range(0,rnk-m_rnk)}

        # Then if taxa is valid then add taxa as final entry
        if valid == True:
            red_taxa[f"{BACKBONE_RANKS_EX[rnk]}"] = (str(tax_dic["ScientificName"]),
            int(tax_dic["TaxId"]))

        # Create and populate microbial taxon
        mtaxon = NCBITaxon(name=nnme,genbank_id=genbank_id,taxa_hier=red_taxa)

        # Check for non-backbone rank cases
        if valid == False:
            # Store (non-standard) taxon rank to explain divergence
            t = tax_dic["Rank"]
            mtaxon.diverg=f"{t}"

        # Check if AkaTaxIds exists in taxonomic information
        if 'AkaTaxIds' in tax_dic.keys():
            # Warn user that they've given a superseeded taxa ID
            LOGGER.warning(f"NCBI ID {(tax_dic['AkaTaxIds'])[0]} has been "
                            f"superseeded by ID {tax_dic['TaxId']}")
            # Record that a superseeded GenBank ID has been provided
            mtaxon.superseed = True

        return mtaxon

    # New function to read in taxa information
    def taxa_search(self, nnme: str, taxah: dict):
        """Method that takes in taxonomic information, and finds the corresponding
        genbank ID. This genbank ID is then used to generate a NCBITaxon object,
        which is returned. This function also makes use of parent taxa information
        to distinguish between ambigious taxa names.

        Params:
            taxah: A dictonary containing taxonomic information

        Returns:
            A NCBITaxon object
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
        try:
            handle = Entrez.esearch(db="taxonomy", term=s_term, api_key=user_key)
        except urllib.error.HTTPError:
            NCBIError('Connection error to remote server')

        # If it works then save the response
        record = Entrez.read(handle)
        handle.close()

        # Store count of the number of records found
        c = int(record['Count'])

        # Check that a singular record has been provided before proceeding
        if c == 1:
            # Find taxa ID as single entry in the list
            tID = int(record['IdList'][0])
            # Use ID lookup function to generate as a NCBITaxon object
            mtaxon = self.id_lookup(nnme,tID)
        # Catch cases where no record is found
        elif c == 0:
            LOGGER.error(f'Taxa {nnme} cannot be found')
            return
        else:
            # Preallocate container for the ranks
            t_ranks = []
            # First check if multiple taxa have the same rank
            for i in range(c):
                # Find taxa ID as single entry in the list
                tID = int(record['IdList'][i])
                # Use ID lookup function to generate a tempoary taxon
                temp_taxon = self.id_lookup(nnme,tID)
                # Add rank to list
                t_ranks.append(list(temp_taxon.taxa_hier.keys())[-1])

            # Record whether ranks match expected rank
            mtch = [True if x == f_key else False for x in t_ranks]

            # If there's only one match, we have determined the correct entry
            if sum(mtch) == 1:
                # Find relevant index
                ind = ([i for i, x in enumerate(mtch) if x])[0]
                # Find taxa ID as single entry in the list
                tID = int(record['IdList'][ind])
                # Use ID lookup function to generate as a NCBITaxon object
                mtaxon = self.id_lookup(nnme,tID)
            # Then check whether multiple taxonomic levels have been provided
            elif len(taxah) == 1:
                # If not raise an error
                LOGGER.error(f'Taxa {nnme} cannot be found using only one '
                             f'taxonomic level, more should be provided')
                return
            else:
                # Find second from last dictonary key
                f_key = list(taxah.keys())[-2]
                # Then find corresponding entry as a search term
                s_term = taxah[f_key]

                # Then find details of this parent record
                try:
                    handle = Entrez.esearch(db="taxonomy", term=s_term,
                                            api_key=user_key)
                except urllib.error.HTTPError:
                    NCBIError('Connection error to remote server')

                p_record = Entrez.read(handle)
                handle.close()

                # Store count of the number of records found
                pc = int(p_record['Count'])

                # Check that single parent taxa exists in records
                if pc == 0:
                    LOGGER.error(f'Provided parent taxa for {nnme} not found')
                    return
                elif pc > 1:
                    LOGGER.error(f'More than one possible parent taxa for {nnme} found')
                    return
                else:
                    # Find parent taxa ID as single entry in the list
                    pID = int(p_record['IdList'][0])
                    # Then use ID lookup function to find generate as a NCBITaxon object
                    ptaxon = self.id_lookup("parent",pID)

                # Save parent taxa rank and name
                p_key = list(ptaxon.taxa_hier.keys())[-1]
                p_val = ptaxon.taxa_hier[p_key]

                # Use list comprehension to make list of potential taxa
                potents = [self.id_lookup(f"c{i}",int(record['IdList'][i]))
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
                    tID = int(record['IdList'][child.index(True)])
                    # Use ID lookup function to find generate as a NCBITaxon object
                    mtaxon = self.id_lookup(nnme,tID)

        # Find last dictonary key
        f_key = list(mtaxon.taxa_hier.keys())[-1]

        # Check if taxonomic rank supplied is used
        if mtaxon.diverg == None:
            # Check that this is also the last dictonary key that was supplied
            if f_key != list(taxah.keys())[-1]:
                # If not raise an error
                LOGGER.error(f'{list(taxah.values())[-1]} is a {f_key}'
                             f' not a {list(taxah.keys())[-1]}')
                return
            # Then check whether orginally supplied name is still used
            elif taxah[f_key] != mtaxon.taxa_hier[f_key]:
                # If not print a warning
                LOGGER.warning(f'{taxah[f_key]} not accepted usage should be '
                               f'{mtaxon.taxa_hier[f_key]} instead')
                # And record the fact that usage is superseeded
                mtaxon.superseed = True
        else:
            # Check that the divergence rank matches the initially supplied rank
            if mtaxon.diverg != list(taxah.keys())[-1]:
                # If not raise an error
                LOGGER.error(f'{list(taxah.values())[-1]} is a {mtaxon.diverg}'
                             f' not a {list(taxah.keys())[-1]}')
                return

        return mtaxon

# OKAY SO WHAT DO I NEED THIS TOP LEVEL FUNCTION TO DO?
# FIRST IT NEEDS TO READ IN WORKSHEET DATA
# IT SHOULD STORE IT IN NCBI FORMAT (AND VALIDATE)
# THEN IT SHOULD CONVERT THIS INFO TO GBIF FORM
# THIS INFO IS STORED FOR OUTPUT
class GenBankTaxa:
    # REWRITE THIS HEAVILY
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

    def __init__(self):
        """Sets the initial properties of a GenBankTaxa object
        """

        # WHAT INFO FROM NCBI SHOULD BE STORED???
        # name: str ???
        # taxa_hier: dict NO NEED TO STORE FULL hierarchy
        # genbank_id: YES
        # synonyms: list[str] THESE ARE ALSO ONLY OF TRANSIENT INTEREST
        # diverg: str = dataclasses.field(init=False) # SHOULD GENERATE A WARNING
        # superseed: # SHOULD GENERATE A WARNING

        self.taxon_index = []
        self.taxon_names = set()
        self.ncbi_t = dict()
        self.hierarchy = set()
        self.n_errors = None
        self.taxon_names_used = set()

        # At the moment we only have one validator defined
        self.validator = RemoteNCBIValidator()

    # def load():
    # FUNCTION TO LOAD IN A TAXA WORKSHEET SHOULD BE ADDED IN HERE

    def validate_and_add_taxon(self, gb_taxon_input):
        # REWRITE THIS TO MATCH ALTERED FUNCTION
        """ User information is provided that names a taxon and (optionally) gives
        a NCBI taxonomy ID. This information is then used to find the closest GBIF
        backbone compatible entry in the NCBI database. This information along with
        the NCBI ID for the orginal entry is then saved, along with a list of possible
        synoyms. All this information is then used to find the closest taxa match
        in the GBIF database. This information is then used to update the relevant
        GenBankTaxa instance.

        The gb_taxon_input has the form:

        ['worksheet_name',
         ['name', 'taxa_hier'],
          'genbank_id']

        If there is no NCBI ID, the structure is:
        ['worksheet_name',
         ['name', 'taxa_hier'],
          None]

        Args:
            taxon_input: Taxon information in standard form as above

        Returns:
            IS THE BELOW CORRECT FOR MY CASE???
            Updates the taxon_names and taxon_index attributes of the class instance.
        """

        m_name, taxon_hier, genbank_id = gb_taxon_input

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
        if not(genbank_id is None or
               (isinstance(genbank_id, float) and genbank_id.is_integer()) or
               isinstance(genbank_id, int)) :
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

        # Now that inputs are sanitised, continue with checking...
        # First gather the info needed to index the entry
        if genbank_id != None:
            ncbi_info = [list(taxon_hier.values())[-1], list(taxon_hier.keys())[-1],
                         int(genbank_id)]
        else:
            ncbi_info = [list(taxon_hier.values())[-1], list(taxon_hier.keys())[-1],
                         None]

        # Set pre-processing as initially false
        p_proc = False

        # Check if pair has already been processed
        if tuple(ncbi_info) in self.ncbi_t:
            hr_taxon = self.ncbi_t[tuple(ncbi_info)]
            p_proc = True
        else:
            # If not go straight ahead and search for the taxon
            hr_taxon = self.validator.taxa_search(m_name, taxon_hier)
            # Catch case where errors are returned rather than a taxon
            if hr_taxon == None:
                LOGGER.error(f'Search based on taxon hierarchy failed')
                return

            # Then check if a genbank ID number has been provided
            if genbank_id != None:
                id_taxon = self.validator.id_lookup(m_name, int(genbank_id))
                # Check whether this matches what was found using hierarchy
                if id_taxon != hr_taxon:
                    # Check if taxonomy hierarchy superseeded
                    if hr_taxon.superseed == True:
                        LOGGER.warning(f'Taxonomic classification superseeded for '
                        f'{m_name}, using new taxonomic classification')
                    elif id_taxon.superseed == True:
                        LOGGER.warning(f'NCBI taxa ID superseeded for {m_name}'
                        f', using new taxa ID')
                    else:
                        LOGGER.error(f"The NCBI ID supplied for {m_name} does "
                        f"not match hierarchy: expected {hr_taxon.genbank_id}"
                        f" got {genbank_id}")
                        return
            else:
                # Warn user that superseeded taxonomy won't be used
                if hr_taxon.superseed:
                    LOGGER.warning(f'Taxonomic classification superseeded for '
                    f'{m_name}, using new taxonomic classification')

            # Store the NCBI taxon keyed by NCBI information (needs tuple)
            self.ncbi_t[tuple(ncbi_info)] = hr_taxon

        # Check if this has succeded without warnings or errors
        if hr_taxon.superseed == False and hr_taxon.diverg == None:
            # Straight forward when there isn't a genbank id, or previously processed
            if genbank_id == None or p_proc == True:
                # If so inform the user of this
                LOGGER.info(f'Taxon ({m_name}) found in NCBI database')
            # Otherwise need to check for superseeded ID's
            elif id_taxon.superseed == False:
                LOGGER.info(f'Taxon ({m_name}) found in NCBI database')

        return

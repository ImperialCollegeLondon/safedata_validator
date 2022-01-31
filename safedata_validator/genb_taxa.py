from typing import Union, Optional
import dataclasses
from Bio import Entrez

import requests
from enforce_typing import enforce_types

from safedata_validator.logger import (CH, FORMATTER, LOGGER,
                                       loggerinfo_push_pop)

# CHANGE THIS TO A PURPOSE DEFINED EMAIL AT SOMEPOINT
# MAYBE ADD TO RESOURCE FILE, AS A CHECK THAT USER HAS PROVIDED ONE
Entrez.email = "jc2017@ic.ac.uk"

# Extended version of backbone ranks to capture superkingdoms
# WILL NEED THOUGHT ABOUT THE MAPPING
BACKBONE_RANKS_EX = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                    'family', 'genus', 'species', 'subspecies']

"""
This module describes classes and methods used to both assess the validity of
GenBank taxonic data, and to check that the taxonomic data is consistent with
that stored in GBIF. When a conflict arises the user will be informed and asked
to reduce the level of taxonomic detail they provide to the point where GenBank
and GBIF are in agreement.

The NCBITaxon dataclass is used to store data about a taxon entry in a dataset. It
is initialised with user data and then the RemoteNCBIValidator class can be used
to update a Taxon object with the result of NCBI (GenBank) validation. Further
information (e.g. synonymus names) is also extract from the remote database, in
order to help with the next validation step.

FURTHER DETAILS OF STUFF DEFINED.

MORE SPIEL (ANYTHING THAT ITS IMPORTANT THAT THE END USER NOTES).
"""

class NCBIError(Exception):
    """Exception class for remote NCBI errors

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message="NCBI taxa ID not found"):
        self.message = message
        super().__init__(self.message)

# TODO - Unit testing, work out how I set up unit tests. Probably best begun early
# Worth asking David how to do this if I can't work it out myself

# TODO - Validate against GBIF
# So check if taxa provided exists if GBIF, if not check up hierachy until one that does is found
# Then tell user that they have to contract their taxonomic specification to this levels

# TODO - Link NCBI and GBIF steps
# Some kind of higher level function that links NCBI functions with GBIF functions
# into a coherent whole.

# TODO - Correctly read in xslx data
# Lot of this has already been written by David, main thing is that I will have
# to decide on the formatting, like what data are we taking in.

# TODO - Think about data output
# Does this make a GBIF compatible csv (or xslx) file?
# Or does it just validate that the provided data is GBIF compatible?

# POTENTIAL - Take steps to increase the speed
# Could involve using an api to speed up entrez
# Or alternatively by downloading the relevant part of NCBI's database and
# running the validation locally

# TODO - Modify the resource file to ask the user to provide an email address
# This should only be done if the user actually wants to use this module as it
# isn't need elsewhere (as far as I know)

# QUESTIONS FOR DAVID
# SHOULD A YOU ARE NOT CONNCTED TO THE INTERNET ERROR BE SETUP?

@enforce_types
@dataclasses.dataclass
class NCBITaxon:
    """Holds taxonomic information from a user on a microbial taxa. This can be
    populated using NCBI GenBank validation.

    There are 3 class properties that can be used to create an instance:
        * name
        * taxa_hier: Dictionary of valid taxonomic hierachy
        * genbank_id: GenBank ID for full taxa (i.e. including non-backbone ranks)
    The remaining properties are populated by processing functions not when an
    instance is created.
        * diverg: does the GenBank ID and stored taxa info diverge, and if so how?
        * synonyms: list of synonymus names for the taxa
        * superseed: is supplied taxon name/ID still the accepted usage
    """

    # Init properties
    name: str
    taxa_hier: dict
    genbank_id: Optional[Union[int, float]] = None
    diverg: str = dataclasses.field(init=False)
    synonyms: list[str] = dataclasses.field(init=False)
    superseed: bool = dataclasses.field(init=False)

    def __post_init__(self):
        """Sets the defaults for the post-init properties and checks inputs
        """

        if self.genbank_id is not None:
            if isinstance(self.genbank_id, float) and not self.genbank_id.is_integer():
                raise ValueError('GenBank Id is not an integer')
            self.genbank_id = int(self.genbank_id)

        self.diverg = None
        self.synonyms = []
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
            A MicTaxon object
        """

        if not isinstance(genbank_id, int):
            raise TypeError()

        if not genbank_id > 0:
            raise ValueError()

        # Set status of taxa validity as initally false
        valid = False

        # Use efetch to find taxonomy details based on the index
        try:
            handle = Entrez.efetch(db="taxonomy",
                                   id=f"{genbank_id}",
                                   retmode="xml")
        except urllib.error.HTTPError:
            NCBIError('Connection error to remote server')

        # Extract as a single dictonary and close
        tax_dic = Entrez.read(handle)[0]
        handle.close()

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
        red_taxa = {f"{actual_bb_rnks[i]}":linx[vinds[i]]["ScientificName"]
                    for i in range(0,rnk-m_rnk)}

        # Then if taxa is valid then add taxa as final entry
        if valid == True:
            red_taxa[f"{BACKBONE_RANKS_EX[rnk]}"] = tax_dic["ScientificName"]

        # Create and populate microbial taxon
        mtaxon = MicTaxon(name=nnme,genbank_id=genbank_id,taxa_hier=red_taxa)

        # Check for non-backbone rank cases
        if valid == False:
            # Store (non-standard) taxon rank to explain divergence
            t = tax_dic["Rank"]
            mtaxon.diverg=f"{t}"

        # Check that OtherNames exist in the dictonary
        if 'OtherNames' in tax_dic.keys():
            # If so find and add synonyms to potentially search GBIF for
            mtaxon.synonyms = (tax_dic["OtherNames"])["Synonym"]

        # Check if AkaTaxIds exists in taxonomic information
        if 'AkaTaxIds' in tax_dic.keys():
            # Warn user that they've given a superseeded taxa ID
            LOGGER.warning(f"NCBI ID {(tax_dic['AkaTaxIds'])[0]} has been "
                            f"superseeded by ID {tax_dic['TaxId']}")
            # Record that a superseeded GenBank ID has been provided
            mtaxon.superseed = True

        return mtaxon

    # New function to read in taxa information
    def taxa_search(self, nnme: str, taxa: dict):
        """Method that takes in taxonomic information, and finds the corresponding
        genbank ID. This genbank ID is then used to generate a MicTaxon object,
        which is returned. This function also makes use of parent taxa information
        to distinguish between ambigious taxa names.

        Params:
            taxa: A dictonary containing taxonomic information

        Returns:
            A MicTaxon object
        """

        # Find last dictonary key
        f_key = list(taxa.keys())[-1]

        # Then find corresponding entry as a search term
        s_term = taxa[f_key]

        # Search the online database
        try:
            handle = Entrez.esearch(db="taxonomy", term=s_term)
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
            # Use ID lookup function to find generate as a MicTaxon object
            mtaxon = self.id_lookup(nnme,tID)
        # Catch cases where no record is found
        elif c == 0:
            LOGGER.error(f'Taxa {nnme} cannot be found')
            return
        else:
            # Check whether multiple taxonomic levels have been provided
            if len(taxa) == 1:
                # If not raise an error
                LOGGER.error(f'Taxa {nnme} cannot be found using only one '
                             f'taxonomic level, more should be provided')
                return
            else:
                # Find second from last dictonary key
                f_key = list(taxa.keys())[-2]
                # Then find corresponding entry as a search term
                s_term = taxa[f_key]

                # Then find details of this parent record
                try:
                    handle = Entrez.esearch(db="taxonomy", term=s_term)
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
                    # Then use ID lookup function to find generate as a MicTaxon object
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
                    # Use ID lookup function to find generate as a MicTaxon object
                    mtaxon = self.id_lookup(nnme,tID)

        # Check if taxonomic rank supplied is used
        if mtaxon.diverg == None:
            # Find last dictonary key
            f_key = list(taxa.keys())[-1]
            # Then check whether orginally supplied name is still used
            if taxa[f_key] != mtaxon.taxa_hier[f_key]:
                # If not print a warning
                LOGGER.warning(f'{taxa[f_key]} not accepted usage should be '
                               f'{mtaxon.taxa_hier[f_key]} instead')
                # And record the fact that usage is superseeded
                mtaxon.superseed = True

        return mtaxon


# THIS IS SOMEWHAT LIKE THE STRUCTURE THAT DAVID SUGGESTED FOR UNIT TESTING
# Make validator
val = RemoteNCBIValidator()
# E coli (562)
d1 = {'genus': 'Escherichia', 'species': 'Escherichia coli'}
# Enterobacteria family (543)
d2 = {'order': 'Enterobacterales', 'family': 'Enterobacteria'}
# E coli strain (1444049)
d3 = {'species': 'Escherichia coli', 'strain': 'Escherichia coli 1-110-08_S1_C1'}
# Streptophytina subphylum (131221)
d4 = {'phylum': 'Streptophyta', 'subphylum': 'Streptophytina'}
# Opisthokonta clade (33154)
d5 = {'superkingdom': 'Eukaryota', 'clade': 'Opisthokonta'}
# Vulpes vulpes (9627)
d6 = {'genus': 'Vulpes', 'species': 'Vulpes vulpes'}
# Morus (NA)
d7 = {'genus': 'Morus'}
# Moraceae morus (3497)
d8 = {'family': 'Moraceae', 'genus': 'Morus'}
# Sulidae morus (37577)
d9 = {'family': 'Sulidae', 'genus': 'Morus'}
# Chordata morus (NA)
d10 = {'phylum': 'Chordata', 'genus': 'Morus'}
# Eukaryota morus(NA)
d11 = {'superkingdom': 'Eukaryota', 'genus': 'Morus'}
# Microcopris hidakai (2602157)
d12 = {'genus': 'Microcopris', 'species': 'Microcopris hidakai'}
# Nonsense garbage (NA)
d13 = {'genus': 'Nonsense', 'species': 'Nonsense garbage'}
# Tenacibaculum maritimum (107401)
d14 = {'genus': 'Tenacibaculum', 'species': 'Tenacibaculum maritimum'}
# Cytophaga marina (1000)
d15 = {'genus': 'Cytophaga', 'species': 'Cytophaga marina'}
# Maybe find synonym of this one to test as well
# Then test output of taxa search
test = val.taxa_search("Nickname",d15)
# Search for ID
# test = val.id_lookup("Nickname",1000)
print(test)

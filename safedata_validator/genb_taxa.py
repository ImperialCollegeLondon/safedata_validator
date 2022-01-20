from typing import Union, Optional
import dataclasses
import requests
from enforce_typing import enforce_types
from Bio import Entrez
import xml.dom.minidom

# CHANGE THIS TO A PURPOSE DEFINED EMAIL AT SOMEPOINT
Entrez.email = "jc2017@ic.ac.uk"

# Extended version of backbone ranks to capture superkingdoms
BACKBONE_RANKS_EX = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                    'family', 'genus', 'species', 'subspecies']

"""
This module describes classes and methods used to both assess the validity of
GenBank taxonic data, and to check that the taxonomic data is consistent with
that stored in GBIF. When a conflict arises the user will be informed and asked
to reduce the level of taxonomic detail they provide to the point where GenBank
and GBIF are in agreement.

DETAILS OF STUFF DEFINED.

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

# TODO - So we want a GenBank taxonomy to be provided
# From name alone we should be able to find an ID
# But probably worth asking for a taxon ID in case of name conflicts
# First stage is basically does this exist in GenBank at all

# TODO - Synonym checking
# GenBank lists hetrotypic synonyms this can be used for synonym checking
# Problem is what if the synonyms preferred with GenBank are not those preferred by GBIF?
# Can we save and store all synonyms and test them all in that case?

# TODO - Validate against GBIF
# So check if taxa provided exists if GBIF, if not check up hierachy until one that does is found
# Then tell user that they have to contract their taxonomic specification to this levels

# QUESTIONS FOR DAVID
# WHERE SHOULD WARNINGS BE SENT TO? HALF SORTED THIS, BUT STILL NEED TO WORK OUT THE LOGGER
# WHAT SHOULD WE DO ABOUT SUPERKINGDOM, DOMAIN, KINGDOM ISSUE?
# DO WE WANT A LocalNCBIValidator?
# HOW DO I ACTUALLY SET UP TESTING?
# HOW TO HANDLE SPECIES NAMES, WHICH ARE OFTEN GIVEN AS E.G. VULPES VULPES
# CAN/SHOULD WE SET UP AN EMAIL?
# SHOULD A YOU ARE NOT CONNCTED TO THE INTERNET ERROR BE SETUP?

@enforce_types
@dataclasses.dataclass
class MicTaxon:
    """Holds taxonomic information from a user on a microbial taxa. This can be
    populated using NCBI GenBank validation.

    There are 3 class properties that can be used to create an instance:
        * name
        * taxa_hier: Dictionary of valid taxonomic hierachy
        * genbank_id: GenBank ID for full taxa (i.e. including non-backbone ranks)
    The remaining properties are populated by processing functions not when
    an instance is created.
        * diverg: does the GenBank ID and stored taxa info diverge, and if so how?
    """

    # Init properties
    name: str
    taxa_hier: dict
    genbank_id: Optional[Union[int, float]] = None
    diverg: str = dataclasses.field(init=False)

    def __post_init__(self):
        """Sets the defaults for the post-init properties and checks inputs
        """

        if self.genbank_id is not None:
            if isinstance(self.genbank_id, float) and not self.genbank_id.is_integer():
                raise ValueError('GenBank Id is not an integer')
            self.genbank_id = int(self.genbank_id)

        self.diverg = None

@enforce_types
class RemoteNCBIValidator:
    # ADD MORE COMMENTS HERE AS AND WHEN I FIGURE THEM OUT
    """This provides a validate method for a MicrobeTaxon using various online
    NCBI APIs. This doesn't need an __init__ method and just contains methods.
    """
    # Functionality to find taxa information from genbank ID
    def id_lookup(self, nnme: str, genbank_id: int):
        """Method to return full taxonomic information from a GenBank ID.

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
        handle = Entrez.efetch(db="taxonomy",
                               id=f"{genbank_id}",
                               retmode="xml")
        # urllib.error.HTTPError: HTTP Error 500: Internal Server Error
        # MAKE SURE THAT THIS ERROR CAN BE CAUGHT, BUT DON'T KNOW HOW TO DO THIS WITHOUT THE ERROR JUST FIRING

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
                    raise NCBIError("""NCBI taxa ID cannot be mapped onto
                    backbone ranks""")
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

        return mtaxon

    # HOW ARE SYNONYMS HANDLED?
    # New function to read in taxa information
    def taxa_search(self, nnme: str, taxa: dict):
        """Method to find GenBank ID from taxonomic information.

        Params:
            taxa: A dictonary containing taxonomic information

        Returns:
            A MicTaxon object
        """

        # Find last dictonary key
        f_key = list(taxa.keys())[-1]

        # Then find corresponding entry as a search term
        s_term = taxa[f_key]

        # NETWORK ERRORS!
        # "Raises an IOError exception if there’s a network error"
        handle = Entrez.esearch(db="taxonomy", term=s_term)
        record = Entrez.read(handle)
        print(record)
        handle.close()

        # Store count of the number of records found
        c = int(record['Count'])

        # Check that a singular record has been provided before proceeding
        if c == 1:
            # Find taxa ID as single entry in the list
            tID = int(record['IdList'][0])
            # Use ID lookup function to find generate as a MicTaxon object
            mtaxon = id_lookup(nnme,tID)
        # Catch cases where no record is found
        elif c == 0:
            # NEED TO LOG AN ERROR HERE I RECKON
            print("Not found error")
        else:
            # Check whether multiple taxonomic levels have been provided
            if len(taxa) == 1:
                # If not raise an error
                # LOG ERROR!!!!!!!!!!!
                print("Ambiguity error")
            else:
                # Unsure exactly what to do here
                b = 1000
                # Definetly need to search the higher taxonomic level
                # Verify it exists and that I haven't been fed garbage
                # Store details of this parent taxa
                # Then load in parent taxa details for all potential taxa
                # And then select the one with the right details

        return mtaxon


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
# vulpes vulpes (9627)
d6 = {'genus': 'Vulpes', 'species': 'Vulpes vulpes'}
# Morus (NA)
d7 = {'genus': 'Morus'}
# Moraceae morus (3497)
d8 = {'family': 'Moraceae', 'genus': 'Morus'}
# Sulidae morus (37577)
d9 = {'family': 'Sulidae', 'genus': 'Morus'}
# Microcopris hidakai (2602157)
d10 = {'genus': 'Microcopris', 'species': 'Microcopris hidakai'}
# Nonsense garbage (NA)
d11 = {'genus': 'Nonsense', 'species': 'Nonsense garbage'}
# Look up an ID
# test = val.id_lookup("Nickname",562)
# Then test output of taxa search
test = val.taxa_search("Nickname",d8)

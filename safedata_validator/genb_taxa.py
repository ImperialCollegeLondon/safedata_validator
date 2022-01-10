from typing import Union, Optional
import dataclasses
import requests
from enforce_typing import enforce_types
from Bio import Entrez

# CHANGE THIS TO A PURPOSE DEFINED EMAIL AT SOMEPOINT
Entrez.email = "jc2017@ic.ac.uk"

# Extended version of backbone ranks to capture superkingdoms
BACKBONE_RANKS_EX = ['superkingdom', 'kingdom', 'phylum', 'order', 'class',
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

# TODO - Citation
# Work out how to appropraitely cite and credit anyones work that I make use of

# QUESTIONS FOR DAVID
# WHERE SHOULD WARNINGS BE SENT TO? HALF SORTED THIS, BUT STILL NEED TO WORK OUT THE LOGGER
# WHAT SHOULD WE DO ABOUT SUPERKINGDOM, DOMAIN, KINGDOM ISSUE?
# DO WE WANT A LocalNCBIValidator?
# HOW DO I ACTUALLY SET UP TESTING?
# HOW TO HANDLE SPECIES NAMES, WHICH ARE OFTEN GIVEN AS E.G. VULPES VULPES
# CAN/SHOULD WE SET UP AN EMAIL?

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

        # Make appropriate url
        up_url = f"""http://api.unipept.ugent.be/api/v1/taxonomy.json?input[]=
                 {genbank_id}&extra=true&names=true"""

        # Retrive data using the url
        taxon_row = requests.get(up_url)

        # Filter out problems with accessing the remote server
        if taxon_row.status_code != 200:
            raise NCBIError('Connection error to remote server')

        # Extract the response (as a list of dictionaries)
        response = taxon_row.json()

        # Check if list is empty
        if not response:
            raise NCBIError()

        # Only extract first dictonary in the list
        tax_dic = response[0]

        # Check that the taxon rank provided is a backbone rank
        if tax_dic["taxon_rank"] in BACKBONE_RANKS_EX:
            # In this case use provided rank
            rnk = BACKBONE_RANKS_EX.index(tax_dic["taxon_rank"])

        # Filter out ID's without ranks (e.g. strains)
        else:
            # Set as not a valid taxa
            vld_tax = False
            # Find lowest index
            r_ID = len(BACKBONE_RANKS_EX) - 1

            # While loop that runs until valid taxa is found
            while vld_tax == False:
                # Check if taxa id is found
                if tax_dic[f"{BACKBONE_RANKS_EX[r_ID]}_id"] != None:
                    # Close loop and store rank number
                    vld_tax = True
                    rnk = r_ID
                # Raise error once backbone ranks have been exhausted
                elif r_ID < 1:
                    raise NCBIError("""NCBI taxa ID cannot be mapped onto
                    backbone ranks""")
                else:
                    r_ID -= 1

        # Create dictonary of reduced taxa info using a list
        red_taxa = {f"{BACKBONE_RANKS_EX[i]}":tax_dic[f"{BACKBONE_RANKS_EX[i]}_name"]
                    for i in range(0,rnk+1)}

        # Create and populate microbial taxon
        mtaxon = MicTaxon(name=nnme,genbank_id=genbank_id,taxa_hier=red_taxa)

        # Check for non-backbone rank cases
        if tax_dic["taxon_rank"] not in BACKBONE_RANKS_EX:
            # Check whether taxa is non standard or just non backbone
            if tax_dic["taxon_rank"] == "no rank":
                # No rank, i.e. non-standard like strain or clade
                mtaxon.diverg='no rank'
            else:
                # Otherwise store (non-standard) taxon rank to explain divergence
                t = tax_dic["taxon_rank"]
                mtaxon.diverg=f"{t}"

        return mtaxon

    # HOW ARE SYNONYMS HANDLED?
    # New function to read in taxa information
    def taxa_search(self, taxa: dict):
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

        # Search taxonomy database
        handle = Entrez.esearch(db="taxonomy", term=s_term)
        record = Entrez.read(handle)
        handle.close()

        # Store count of the number of records found
        c = int(record['Count'])
        print(record)

        # Check that a singular record has been provided before proceeding
        if c == 1:
            # Find taxa ID as single entry in the list
            tID = int(record['IdList'][0])
        else:
            # NEED TO LOG AN ERROR HERE I RECKON
            print("NOT GREAT")
            # NEED TO HANDLE MULTIPLE ENTRY CASE HERE

        # NEED TO HANDLE NETWORK ERRORS ETC
        # OR EVEN HOW TO READ ALL THIS
        # YEAH NEED TO WORK OUT HOW THIS FUNCTION EVEN WORKS
        handle = Entrez.efetch(db="taxonomy", id="562", rettype="uilist", retmode="text")
        print(handle.readline().strip())
        handle.close()

        # DO I READ THE DATA IN A DIFFERENT WAY OR JUST USE THE GENBANK ID?

        # Fill in missing taxonomic information

        # Check against backbone
        # Generate instance
        return


# Make validator
val = RemoteNCBIValidator()
# E coli (562)
d1 = {'genus': 'Escherichia', 'species': 'Escherichia coli'}
# E coli strain (1444049)
d2 = {'species': 'Escherichia coli', 'strain': 'Escherichia coli 1-110-08_S1_C1'}
# Streptophytina subphylum (131221)
d3 = {'phylum': 'Streptophyta', 'subphylum': 'Streptophytina'}
# Opisthokonta clade (33154)
# vulpes vulpes (9627)
d4 = {'genus': 'Vulpes', 'species': 'Vulpes vulpes'}
# Moraceae morus (3497)
d5 = {'family': 'Moraceae', 'genus': 'Morus'}
# Sulidae morus (37577)
d6 = {'family': 'Sulidae', 'genus': 'Morus'}
# Look up an ID
# test = val.id_lookup("E coli",1444049)
# Then test output of taxa search
test = val.taxa_search(d2)
# Print out whatever gets returns
print(test)

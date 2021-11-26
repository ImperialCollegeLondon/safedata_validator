from typing import Union, Optional
import dataclasses
import requests
from enforce_typing import enforce_types

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

# Meta questions:
# What do we want to store? GBIF compatable taxonomy? With or without GenBank ID?
# How do we address newly defined species?
# Do we want a LocalNCBIValidator?

# TO DO - So we want a GenBank taxonomy to be provided
# From name alone we should be able to find an ID
# But probably worth asking for a taxon ID in case of name conflicts
# First stage is basically does this exist in GenBank at all

# TO DO - Synonym checking
# GenBank lists hetrotypic synonyms this can be used for synonym checking
# Problem is what if the synonyms preferred with GenBank are not those preferred by GBIF?
# Can we save and store all synonyms and test them all in that case?

# TO DO - Validate against GBIF
# So check if taxa provided exists if GBIF, if not check up hierachy until one that does is found
# Then tell user that they have to contract their taxonomic specification to this levels

# TO DO - Citation
# Work out how to appropraitely cite and credit anyones work that I make use of

@enforce_types
class RemoteNCBIValidator:
    # ADD MORE COMMENTS HERE AS AND WHEN I FIGURE THEM OUT
    """This provides a validate method for a MicrobeTaxon using various online
    NCBI APIs. This doesn't need an __init__ method and just contains methods.
    """
    # Functionality to find taxa information from genbank ID
    def id_lookup(self, genbank_id: int):
        """Method to return full taxonomic information from a GenBank ID. THERE'S PROBABLY MORE TO SAY THOUGH

        Params:
            genbank_id: An integer

        Returns:
            DECIDE WHAT I WANT IT TO RETURN
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

        # ALSO WORK OUT WHAT INFO I SHOULD BE PROVIDING TO THE LOGGER
        # Check that the taxon rank provided is a backbone rank
        if tax_dic["taxon_rank"] in BACKBONE_RANKS_EX:
            # In this case use provided rank
            rnk = tax_dic["taxon_rank"]

        # Filter out ID's without ranks (e.g. strains)
        else:
            # Check whether taxa is non standard or just non backbone
            if tax_dic["taxon_rank"] == "no rank":
                # Check for strains
                print("Either a strain or a clade")
            else:
                t = tax_dic["taxon_rank"]
                print(f"Rank {t} not a backbone rank")

            # Set as not a valid taxa
            vld_tax = False
            # Find lowest index
            r_ID = len(BACKBONE_RANKS_EX) - 1

            # While loop that runs until valid taxa is found
            while vld_tax == False:
                # Check if taxa id is found
                if tax_dic[f"{BACKBONE_RANKS_EX[r_ID]}_id"] != None:
                    # Close loop and store rank
                    vld_tax = True
                    rnk = BACKBONE_RANKS_EX[r_ID]
                # Raise error once backbone ranks have been exhausted
                elif r_ID < 1:
                    raise NCBIError("""NCBI taxa ID cannot be mapped onto
                    backbone ranks""")
                else:
                    r_ID -= 1

        # Loop over all ranks down to the taxon_rank
        for i in range(1+BACKBONE_RANKS_EX.index(rnk)):
            print(i)
            print(BACKBONE_RANKS_EX[i])



        return tax_dic
        # NEED TO WORK OUT WHAT TO DO BELOW HERE
        # Basically want to produce an object that contains all relevant taxonomic information
        # E.g. kingdom to lowest known level
        # From this function want to populate a taxa with ID number and full list of taxa inforamation (backbone only)

        # # Create and populate taxon
        # taxon = Taxon(name=response['canonicalName'],
        #               rank=response['rank'].lower(),
        #               gbif_id=response['nubKey'])
        #
        # # First, set the parent key - in the GBIF API, this is always provided
        # # and always points to the true parent taxon, unlike parent_key in
        # # the simple local database. Note that parentKey is not included in the
        # # response for Kingdom level taxa, hence get().
        # taxon.parent_id = response.get('parentKey')
        # taxon.taxon_status = response['taxonomicStatus'].lower()
        # taxon.lookup_status = 'found'

        # # Add the taxonomic hierarchy from the accepted usage - these are tuples
        # # to be used to extend a set for the taxonomic hierarchy
        # taxon.hierarchy = [(rk, response[ky])
        #                    for rk, ky in [(r, r + 'Key') for r in BACKBONE_RANKS[:-1]]
        #                    if ky in response]
        #
        # # Now populate the canon details, which requires another look up if the
        # # user details are for a synonym etc.
        # if taxon.taxon_status in ('accepted', 'doubtful'):
        #     taxon.is_canon = True
        # else:
        #     taxon.is_canon = False
        #     # acceptedKey is not provided in the response for canon taxa.
        #     taxon.canon_usage = self.id_lookup(response['acceptedKey'])
        #
        # return taxon

# Make validator
val = RemoteNCBIValidator()
# Look up an ID
# E coli (562)
# E coli strain (1444049)
# Streptophytina subphylum (131221)
# Opisthokonta clade (33154)
test = val.id_lookup(131221)
# Print out whatever gets returns
print(test)

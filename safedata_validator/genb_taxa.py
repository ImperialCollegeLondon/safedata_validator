from typing import Union, Optional
import dataclasses
import requests
from enforce_typing import enforce_types

# DELETE THE BELOW WHEN I AM DONE, AND INSTEAD IMPORT FROM TAXA
BACKBONE_RANKS = ['kingdom', 'phylum', 'order', 'class', 'family',
                  'genus', 'species', 'subspecies']
# from safedata_validator.taxa import BACKBONE_RANKS

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
# Do we want a LocalNCBIValidator

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

@enforce_types
class RemoteNCBIValidator:
    # ADD MORE COMMENTS HERE AS AND WHEN I FIGURE THEM OUT
    """This provides a validate method for a MicrobeTaxon using various online
    NCBI APIs. This doesn't need an __init__ method and just contains methods.
    """
    # Functionality to find taxa information from genbank ID
    def id_lookup(self, genbank_id: int):
        """Method to return taxonomic information from a GenBank ID. THERE'S PROBABLY MORE TO SAY THOUGH

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

        # Extract the response
        response = taxon_row.json()

        # Check if list is empty
        if not response:
            raise NCBIError()

        return response
        # NEED TO WORK OUT WHAT TO DO BELOW HERE

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
test = val.id_lookup(3333)
# Print out whatever gets returns
print(test)

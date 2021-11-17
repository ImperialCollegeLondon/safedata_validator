# Bunch of imports and similar

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

# IMMEDIATE TASKS
# 1) WORK OUT HOW TO INTERFACE WITH GENBANK DATABASE
# 2) WRITE FUNCTION THAT CHECKS IF DATA EXISTS IN GENBANK

# Meta questions:
# What do we want to store? GBIF compatable taxonomy? With or without GenBank ID?
# How do we address newly defined species?

# So we want a GenBank taxonomy to be provided
# From name alone we should be able to find an ID
# But probably worth asking for a taxon ID in case of name conflicts
# First stage is basically does this exist in GenBank at all

# TO DO - Synonym checking
# GenBank lists hetrotypic synonyms this can be used for synonym checking
# Problem is what if the synonyms preferred with GenBank are not those preferred by GBIF?
# Can we save and store all synonyms and test them all in that case?

# TO DO - Validate against GBIF
# So check if taxa provided exists if GBIF, if not check up hierachy until one that does is found
# Then tell user that they have to contract their taxonomic specification to this level

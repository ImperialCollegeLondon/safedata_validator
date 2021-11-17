from typing import Union, Optional
import dataclasses

# SHOULD I DO THIS OR DO I NEED MORE BACKBONE RANKS???
from .taxa import BACKBONE_RANKS

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

# MAKE A NEW CLASS FOR THE MICROBIAL TAXA, AS THEY NEED DIFFERENT FUNCTIONS + DIFFERENT DATA
# THINK I SHOULD BE DRAWING DATA FROM https://www.ncbi.nlm.nih.gov/taxonomy
@enforce_types
@dataclasses.dataclass
class MicrobeTaxon:
    """Holds taxonomic information from a user, which can be populated using GBIF validation.

    There are 4 class properties that can be used to create an instance:
        * name
        * rank
        * genbank_id
        * ignore_gbnk: a GenBank id match to reject
    The remaining properties are populated by processing functions not when
    an instance is created.
        * is_backbone: the taxon is at a taxonomic level included in the GBIF backbone
        * is_canon: the taxon is considered canon in GenBank
        # THESE VALUES PROBABLY CHANGE
        * lookup_status: the outcome of the lookup with one of the following values:
          found, no_match, validation_fail, unknown_id, id_mismatch
        # DON'T KNOW IF THIS WOULD BE RETURNED IN OUR CASE
        * taxon_status: the taxonomic status of the taxon with one of the following values:
          accepted, doubtful, synonym etc. etc.
        # AGAIN UNUSURE IF THIS IS SOMETHING WE WANT IN THE NEW CASE
        * parent_id: a GBIF id for the accepted parent taxon.
        # WILL NEED TO FIND OUT IF ANY EXTRA INFORMATION IS ACTUALLY PROVIDED
        * note: a string of any extra information provided by the search
        # UNSURE ON THIS ONE
        * hierarchy: a list of 2-tuples of rank and GBIF ID for the taxonomic hierarchy
    """

    # Init properties
    name: str
    rank: str
    genbank_id: Optional[Union[int, float]] = None
    ignore_gbnk: Optional[Union[int, float]] = None
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

        if self.genbank_id is not None:
            if isinstance(self.genbank_id, float) and not self.genbank_id.is_integer():
                raise ValueError('GenBank Id is not an integer')
            self.genbank_id = int(self.genbank_id)

        if self.ignore_gbnk is not None:
            if isinstance(self.ignore_gbnk, float) and not self.ignore_gbnk.is_integer():
                raise ValueError('Ignore GenBank is not an integer')
            self.ignore_gbnk = int(self.ignore_gbnk)

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

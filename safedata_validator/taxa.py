"""This module describes classes and methods used to compile taxonomic data from datasets
and to validate taxonomy against the GBIF backbone database and/or the NCBI taxonomy
database.

The two parallel Taxon dataclasses (GBIFTaxon and NCBITaxon) are used to store data
about a taxon entry in a dataset. They are initialised with user data and then the
relevant taxon Validator classes (GBIF or NCBI) can be used to update a Taxon object
with the result of validation against a local version of the taxon databases
(`GBIFValidator` and `NCBIValidator`).

Parallel 'Taxa' worksheets (GBIFTaxa and NCBITaxa) are defined, which are used to load
and collate the set of taxonomic entries from a dataset. These are then collected in a
higher level Taxa object, which additionally records the names used in the Data
worksheets. This allows us to check that all defined names are used, all used names are
defined, and that no names are defined in both Taxa worksheets (if both sheets are
provided).

Note that we explicitly exclude form and variety from the set of GBIF backbone taxonomic
levels because they cannot be matched into the backbone hierarchy without extra API
calls.

When validating against the NCBI database supplied taxa of any rank (i.e. strain or
clade) which can be successfully validated will be recorded. However, associated higher
taxa will only be recorded if their ranks are either a GBIF backbone rank or
superkingdom.
"""  # noqa D415

import dataclasses
import sqlite3
from collections import Counter
from io import StringIO
from itertools import compress, groupby
from typing import Optional, Union

from dominate import tags
from dominate.util import raw
from openpyxl import worksheet

from safedata_validator.logger import (
    COUNTER_HANDLER,
    FORMATTER,
    LOGGER,
    loggerinfo_push_pop,
)
from safedata_validator.resources import Resources
from safedata_validator.validators import (
    GetDataFrame,
    HasDuplicates,
    IsLower,
    blank_value,
)

BACKBONE_RANKS = [
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
]

# Extended version of backbone ranks to capture superkingdoms
BACKBONE_RANKS_EX = ["superkingdom"] + BACKBONE_RANKS

# These are the extra ranks defined by NCBI, they may have to be updated in future
EXTRA_NCBI_RANKS = [
    "biotype",
    "clade",
    "cohort",
    "forma",
    "forma specialis",
    "genotype",
    "infraclass",
    "infraorder",
    "isolate",
    "morph",
    "no rank",
    "parvorder",
    "pathogroup",
    "section",
    "series",
    "serogroup",
    "serotype",
    "species group",
    "species subgroup",
    "strain",
    "subclass",
    "subcohort",
    "subfamily",
    "subgenus",
    "subkingdom",
    "suborder",
    "subphylum",
    "subsection",
    "subtribe",
    "superclass",
    "superfamily",
    "superorder",
    "superphylum",
    "tribe",
    "varietas",
]


class GBIFError(Exception):
    """Exception class for GBIF errors.

    Attributes:
        message: explanation of the error
    """

    def __init__(self, message="GBIF ID not found"):
        self.message = message
        super().__init__(self.message)


class NCBIError(Exception):
    """Exception class for NCBI errors.

    Attributes:
        message: explanation of the error
    """

    def __init__(self, message="No entry found for ID, probably bad ID"):
        self.message = message
        super().__init__(self.message)


@dataclasses.dataclass
class GBIFTaxon:
    """Represent and validate a GBIF taxon.

    Initialised using user taxonomic information for single taxon, which can be then be
    validated against the GBIF database. Attributes are populated when an instance is
    passed to GBIFValidator.

    Args:
        name: A taxonomic name
        rank: A taxonomic rank
        gbif_id: A specific GBIF ID

    Attributes:
        is_backbone: the taxon is at a taxonomic level included in the GBIF backbone
        is_canon: the taxon is considered canon in GBIF
        lookup_status: the outcome of the lookup with one of the following values:
            found, no_match, validation_fail, unknown_id, id_mismatch
        taxon_status: the taxonomic status of the taxon with one of the following
            values: accepted, doubtful, synonym etc. etc.
        parent_id: a GBIF id for the accepted parent taxon.
        canon_usage: a GBIFTaxon instance holding the canonical usage for the taxon
        note: a string of any extra information provided by the search
        hierarchy: a list of 2-tuples of rank and GBIF ID for the taxonomic hierarchy
    """

    # Init properties
    name: str
    rank: str
    gbif_id: Optional[int] = None
    is_backbone: bool = dataclasses.field(init=False)
    is_canon: bool = dataclasses.field(init=False)
    # https://stackoverflow.com/questions/33533148
    canon_usage: Optional["GBIFTaxon"] = dataclasses.field(init=False)
    parent_id: Optional[int] = dataclasses.field(init=False)
    taxon_status: Optional[str] = dataclasses.field(init=False)
    lookup_status: str = dataclasses.field(init=False)
    hierarchy: list = dataclasses.field(init=False)

    def __post_init__(self) -> None:
        """Validates inputs and sets defaults for the post-init properties."""

        if not isinstance(self.name, str):
            raise TypeError("Provided taxon name not a string")

        if not isinstance(self.rank, str):
            raise TypeError("Provided rank not in string form")

        if self.gbif_id is not None:
            if isinstance(self.gbif_id, float) and not isinstance(self.gbif_id, int):
                raise ValueError("GBIF ID is a non-integer float")
            elif not isinstance(self.gbif_id, int):  # Catch non int or float case
                raise TypeError("GBIF ID is neither an int or a float")
            self.gbif_id = int(self.gbif_id)

        self.rank = self.rank.lower()
        self.is_backbone = self.rank in BACKBONE_RANKS
        self.is_canon = False
        self.canon_usage = None
        self.parent_id = None
        self.taxon_status = None
        self.lookup_status = "unvalidated"
        self.hierarchy = []

    def __repr__(self) -> str:
        """Provides a simple representation of the class."""
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
    def found(self) -> bool:
        """Confirms that a taxon is a backbone taxon found in GBIF."""

        # Shorthand property
        return self.is_backbone and self.lookup_status == "found"

    # TODO - could make _eq_ method and allow GBIFTaxon to be used directly
    #        instead of using taxon tuples as data and parent index.


@dataclasses.dataclass
class NCBITaxon:
    """Represent and validate an NCBI taxon.

    Initialised using user taxonomic information for single taxon, which can be then be
    validated against the NCBI database. Attributes are populated when an instance is
    passed to NCBIValidator.

    Args:
        name: A taxonomic name
        rank: A taxonomic rank
        taxa_hier: A dictionary of valid taxonomic hierarchy (with ID's)
        ncbi_id: NCBI ID for full taxa (i.e. including non-backbone ranks)

    Attributes:
        superseed: is supplied taxon name/ID still the accepted usage
        orig: Records rank of the original taxa if it was not found
    """

    # Init properties
    name: str
    rank: str
    taxa_hier: dict
    ncbi_id: Optional[int] = None
    superseed: bool = dataclasses.field(init=False)
    orig: str = dataclasses.field(init=False)

    def __post_init__(self):
        """Sets the defaults for the post-init properties and checks inputs."""

        if not isinstance(self.name, str):
            raise TypeError("Provided taxon name not a string")

        if not isinstance(self.rank, str):
            raise TypeError("Provided rank not in string form")

        if not isinstance(self.taxa_hier, dict):
            raise TypeError("Taxonomic hierarchy not provided in dictionary form")

        if self.ncbi_id is not None:
            if isinstance(self.ncbi_id, float) and not self.ncbi_id.is_integer():
                raise TypeError("NCBI ID is a non-integer float")
            elif not isinstance(self.ncbi_id, int):  # Catch non int or float case
                raise TypeError("NCBI ID is neither an int or a float")
            self.ncbi_id = int(self.ncbi_id)

        if self.taxa_hier is not None:
            if len(self.taxa_hier) == 0:
                raise ValueError("Taxa hierarchy dictionary empty")
            elif all(isinstance(x, str) for x in self.taxa_hier.keys()) is False:
                raise ValueError("Not all taxa dictionary keys are strings")
            elif all(isinstance(x, tuple) for x in self.taxa_hier.values()) is False:
                raise ValueError("Not all taxa dictionary values are tuples")
            elif (
                all(
                    list(map(type, x)) == [str, int, int]
                    or list(map(type, x)) == [str, int, type(None)]
                    for x in self.taxa_hier.values()
                )
                is False
            ):
                raise ValueError("Taxa tuples not all in [string integer integer] form")

        if self.rank.lower() != list(self.taxa_hier.keys())[-1]:
            raise ValueError(
                f"Provided rank ({self.rank.lower()}) does not match lowest rank in "
                f"supplied hierarchy ({list(self.taxa_hier.keys())[-1]})"
            )

        if self.name != (list(self.taxa_hier.values())[-1])[0]:
            raise ValueError(
                f"Provided taxon name ({self.name}) does not match"
                f" name of the lowest entry in supplied hierarchy ("
                f"{(list(self.taxa_hier.values())[-1])[0]})"
            )

        if (
            self.ncbi_id is not None
            and self.ncbi_id != (list(self.taxa_hier.values())[-1])[1]
        ):
            raise ValueError(
                f"Provided NCBI ID ({self.ncbi_id}) does not match"
                f" ID of the lowest entry in supplied hierarchy ("
                f"{(list(self.taxa_hier.values())[-1])[1]})"
            )

        self.rank = self.rank.lower()
        self.superseed = False
        self.orig = None

    def __repr__(self):
        """Provides a simple representation of the class."""

        if self.orig is not None:
            return f"{self.name} (resolved as {self.rank} rather than {self.orig})"
        elif self.superseed:
            return f"{self.name} (superseded taxon details provided)"
        else:
            return f"{self.name}"


def gen_invalid_NCBITaxon() -> NCBITaxon:
    """Generates a standardised invalid NCBITaxon instance."""

    return NCBITaxon("INVALID", "invalid_rank", {"invalid_rank": ("INVALID", -9, None)})


class GBIFValidator:
    """Validate taxon data against the GBIF database.

    This class connects to a local copy of the GBIF database and provides methods to
    validate GBIFTaxon instances and look up GBIF ID values.

    Args:
        resources: A Resources instance linking to the local GBIF database
    """

    def __init__(self, resources: Resources) -> None:
        conn = sqlite3.connect(resources.gbif_database)
        conn.row_factory = sqlite3.Row
        self.gbif_conn = conn

    def __del__(self) -> None:
        """Delete a LocalGBIFValidator instance.

        This method ensures that the database connection is closed correctly.
        """
        self.gbif_conn.close()

    def search(self, taxon: GBIFTaxon) -> GBIFTaxon:
        """Validate a GBIFTaxon instance.

        The method looks for the taxon in the GBIF database using name and rank and
        an optional GBIF ID for disambiguation. The input GBIFTaxon is updated in place
        and so there is no value returned.

        Args:
            taxon: A GBIFTaxon instance
        """

        if not taxon.is_backbone:
            raise ValueError("Cannot validate non-backbone taxa")

        if taxon.gbif_id is not None:
            # get the record associated with the provided ID
            try:
                id_taxon = self.id_lookup(taxon.gbif_id)
            except GBIFError as err:
                taxon.lookup_status = f"GBIF ID problem: {err.message}"
                return taxon

            # Check that name and rank are congruent with id
            if (id_taxon.name != taxon.name) or (id_taxon.rank != taxon.rank):
                taxon.lookup_status = "ID does not match name and rank"
                return taxon

            return id_taxon

        else:
            # get the set of records associated with the taxon and rank

            sql = (
                f"select * from backbone where canonical_name ='{taxon.name}' "
                f"and rank= '{taxon.rank.upper()}';"
            )

            taxon_rows = self.gbif_conn.execute(sql).fetchall()
            selected_row = None

            if len(taxon_rows) == 0:
                # No matching rows
                taxon.lookup_status = "No match found"
                return taxon
            elif len(taxon_rows) == 1:
                # one matching row - extract it from the list
                selected_row = taxon_rows[0]
            else:
                # More than one row - try to mimic the preferred hits reported
                # by the GBIF API to select a single hit by looking at the counts
                # of the different statuses.

                # First, get the taxon statuses
                tx_status = [tx["status"].lower() for tx in taxon_rows]
                tx_counts = Counter(tx_status)

                if "accepted" in tx_counts.keys():
                    if tx_counts["accepted"] == 1:
                        # Single accepted hits are first preference, and if there are
                        # multiple accepted hits then parent resolution needed.
                        selected_row = taxon_rows[tx_status.index("accepted")]
                elif "doubtful" in tx_counts.keys():
                    if tx_counts["doubtful"] == 1:
                        # Single doubtful hits get next preference - not quite sure
                        # about this! - and if there are multiple accepted hits then
                        # resolution needed.
                        selected_row = taxon_rows[tx_status.index("doubtful")]
                else:
                    # Rows now contain only synonyms (of varying kinds) and
                    # misapplied. Both of these types have accepted usage
                    # values, so look for a unique accepted usage, trapping the
                    # edge case of kingdoms, which have no parent_key.
                    tx_acc = {
                        tx["parent_key"]
                        for tx in taxon_rows
                        if tx["parent_key"] is not None
                    }

                    if len(tx_acc) == 1:
                        # A single accepted usage - pick the first row to index
                        selected_row = taxon_rows[0]

            if selected_row is None:
                # No single row has been accepted as the best, so return no
                # match and a note, as the API interface does.
                taxon.lookup_status = f"Multiple equal matches for {taxon.name}"
                return taxon

            # Should now have a single row for the preferred hit, which can be
            # extracted from the database
            return self.id_lookup(selected_row["id"])

    def id_lookup(self, gbif_id: int) -> GBIFTaxon:
        """Get a GBIFTaxon by GBIF ID.

        This method returns a GBIFTaxon directly from a GBIF ID. It will raise
        a GBIFError if the provided ID cannot be found.

        Args:
            gbif_id: A GBIF ID number.

        Returns:
            A populated GBIFTaxon instance for the ID.
        """

        if not isinstance(gbif_id, int):
            raise ValueError("Non-integer GBIF code")

        if not gbif_id >= 0:
            # 0 is kingdom placeholder for incertae sedis
            raise ValueError("Negative GBIF code")

        # get the record associated with the provided ID
        sql = f"select * from backbone where id = {gbif_id}"
        taxon_row = self.gbif_conn.execute(sql).fetchone()

        # check there is a result and that it is congruent with any
        # provided taxon or rank information
        if taxon_row is None:
            raise GBIFError()

        # Create and populate taxon
        taxon = GBIFTaxon(
            name=taxon_row["canonical_name"],
            rank=taxon_row["rank"].lower(),
            gbif_id=taxon_row["id"],
        )
        taxon.lookup_status = "found"
        taxon.taxon_status = taxon_row["status"].lower()
        # Trap empty string rather None in parent ID - 2016 GBIF doesn't use \N to
        # represent None and root kingdoms end up with parent_key = ''
        taxon.parent_id = (
            None if taxon_row["parent_key"] == "" else taxon_row["parent_key"]
        )

        # Detect deleted taxa - these contain a deletion date and (somewhat oddly)
        # have had hierarchy above phylum removed, so parent taxon points at the phylum
        if taxon_row["status"].lower() == "deleted":
            taxon.taxon_status = "deleted"
            taxon.lookup_status = "Deleted taxon"
            return taxon

        # Add the taxonomic hierarchy, using a mapping of backbone ranks (except
        # subspecies) to backbone table fields. This needs to omit missing keys and
        # more nested taxon levels: so for example a genus will have 'species_key' but
        # it will be None (or possibly an empty string in older backbone versions that
        # use that rather than explicit \\N in conversion)
        taxon.hierarchy = [
            (rk, taxon_row[ky])
            for rk, ky in [(r, r + "_key") for r in BACKBONE_RANKS[:-1]]
            if ky in taxon_row.keys()
            and taxon_row[ky] is not None
            and not taxon_row[ky] == ""
        ]

        # parent key in the local database has the odd property that the parent
        # tax_gbif['parent_key'] does dual duty: points up to parent for canon
        # taxa and 'up' to canon for non-canon taxa, so need to look through both
        # to get the canon and parent populated.
        if taxon.taxon_status in ["accepted", "doubtful"]:
            taxon.is_canon = True
        elif taxon.parent_id is None:
            LOGGER.warning("Non-canon taxa does not have valid parent id")
            taxon.is_canon = False
        else:
            taxon.is_canon = False
            taxon.canon_usage = self.id_lookup(taxon.parent_id)
            taxon.parent_id = taxon.canon_usage.parent_id

        return taxon


class NCBIValidator:
    """Validate taxon data against the NCBI database.

    This class connects to a local copy of the NCBI database and provides methods to
    validate NCBITaxon instances and look up NCBI ID values.

    Args:
        resources: A Resources instance linking to the local NCBI database
    """

    def __init__(self, resources: Resources) -> None:
        conn = sqlite3.connect(resources.ncbi_database)
        conn.row_factory = sqlite3.Row
        self.ncbi_conn = conn

    def __del__(self) -> None:
        """Delete a LocalNCBIValidator instance.

        This method ensures that the database connection is closed correctly.
        """

        self.ncbi_conn.close()

    def id_lookup(self, nnme: str, ncbi_id: int) -> NCBITaxon:
        """Get an NCBITaxon by NCBI ID.

        This method returns a populated NCBITaxon instance, given an NCBI ID and
        nickname.  It will raise a NCBIError if the provided ID cannot be found.

        Args:
            nnme: A nickname to identify the taxon
            ncbi_id: Unique identifier for the taxon

        Returns:
            A populated NCBITaxon instance
        """

        if not isinstance(ncbi_id, int):
            raise TypeError("Non-integer NCBI taxonomy ID")

        if not isinstance(nnme, str):
            raise TypeError("Non-string nickname")

        if not ncbi_id > 0:
            raise ValueError("Negative NCBI taxonomy ID")

        superseed = False

        # Look for record associated with the provided ID
        sql = f"select * from nodes where tax_id = {ncbi_id}"
        taxon_row = self.ncbi_conn.execute(sql).fetchone()

        # If nothing found check if this ID has been merged
        if taxon_row is None:
            sql = f"select * from merge where old_tax_id = {ncbi_id}"
            taxon_row = self.ncbi_conn.execute(sql).fetchone()
            # If it's not found then give bad ID error
            if taxon_row is None:
                raise NCBIError()
            else:
                sql = f"select * from nodes where tax_id = {taxon_row['new_tax_id']}"
                taxon_row = self.ncbi_conn.execute(sql).fetchone()
                superseed = True
                # Warn user that they've given a superseded taxa ID
                LOGGER.warning(
                    f"NCBI ID {ncbi_id} has been superseded by ID "
                    f"{taxon_row['tax_id']}"
                )

        # Extract relevant info from the taxon row
        t_rank = taxon_row["rank"]
        good_id = taxon_row["tax_id"]
        # Then use to find and store name
        sql = f"select * from names where tax_id = {good_id} and "
        sql += "name_class = 'scientific name'"
        name_row = self.ncbi_conn.execute(sql).fetchone()
        t_name = name_row["name_txt"]

        # Then setup loop to find the whole lineage
        lin_fnd = False
        linx = []

        while lin_fnd is False:
            tmp_dic = {}
            # Find node and name of the parent taxon
            sql = f"select * from nodes where tax_id = {taxon_row['parent_tax_id']}"
            taxon_row = self.ncbi_conn.execute(sql).fetchone()
            sql = f"select * from names where tax_id = {taxon_row['tax_id']} and "
            sql += "name_class = 'scientific name'"
            name_row = self.ncbi_conn.execute(sql).fetchone()
            # Store all relevant info
            tmp_dic["TaxID"] = taxon_row["tax_id"]
            tmp_dic["ScientificName"] = name_row["name_txt"]
            tmp_dic["Rank"] = taxon_row["rank"]
            # And add all of it to the Lineage
            linx.append(tmp_dic)
            # End this when the parent taxon is root (ID=1)
            if taxon_row["parent_tax_id"] == 1:
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
            while vld_tax is False:
                # Check if taxa id is found
                if any(
                    [rnks[i] == f"{BACKBONE_RANKS_EX[r_ID]}" for i in range(tx_len)]
                ):
                    # Close loop and store rank number
                    vld_tax = True
                    # Add 1 to the rank as only including lineage in this case
                    rnk = r_ID + 1
                    # Warn user that non-backbone rank has been supplied
                    LOGGER.warning(f"{nnme} of non-backbone rank: {t_rank}")
                # Raise error once backbone ranks have been exhausted
                elif r_ID < 1:
                    LOGGER.error(
                        f"Taxon hierarchy for {nnme} contains no backbone ranks"
                    )
                    return gen_invalid_NCBITaxon()
                else:
                    r_ID -= 1

        # Make list of backbone ranks we are looking for
        actual_bb_rnks = BACKBONE_RANKS_EX[0:rnk]

        # Number of missing ranks initialised to zero
        m_rnk = 0

        # Check that all desired backbone ranks appear in the lineage
        if all(item in rnks for item in actual_bb_rnks) is False:
            # Find all missing ranks
            miss = list(set(actual_bb_rnks).difference(rnks))
            # Count missing ranks
            m_rnk = len(miss)
            # Remove missing ranks from our list of desired ranks
            for i in range(0, m_rnk):
                actual_bb_rnks.remove(miss[i])

        # Find valid indices (e.g. those corresponding to backbone ranks)
        vinds = [idx for idx, element in enumerate(rnks) if element in actual_bb_rnks]

        # Create dictionary of valid taxa lineage using a list
        if len(actual_bb_rnks) != 0:
            red_taxa: dict[str, tuple] = {
                f"{actual_bb_rnks[0]}": (
                    str(linx[vinds[0]]["ScientificName"]),
                    int(linx[vinds[0]]["TaxID"]),
                    None,
                )
            }

            # Recursively add all the hierarchy data in
            for i in range(1, rnk - m_rnk):
                red_taxa[f"{actual_bb_rnks[i]}"] = (
                    str(linx[vinds[i]]["ScientificName"]),
                    int(linx[vinds[i]]["TaxID"]),
                    int(linx[vinds[i - 1]]["TaxID"]),
                )

            # Then add taxa information as a final entry
            red_taxa[f"{t_rank}"] = (
                str(t_name),
                int(good_id),
                int(linx[vinds[-1]]["TaxID"]),
            )
        else:
            # Just make taxa with individual rank if no valid lineage provided
            red_taxa = {f"{t_rank}": ((t_name), int(good_id), None)}

        # Create and populate microbial taxon
        mtaxon = NCBITaxon(
            name=t_name, rank=t_rank, ncbi_id=int(good_id), taxa_hier=red_taxa
        )

        # Record whether a superseded NCBI ID has been provided
        mtaxon.superseed = superseed

        return mtaxon

    def taxa_search(self, nnme: str, taxah: dict) -> NCBITaxon:
        """Find an NCBI taxon given a taxon hierarchy.

        Method that takes in taxonomic information, and finds the corresponding
        NCBI ID. This NCBI ID is then used to generate a NCBITaxon object,
        which is returned. This function also makes use of parent taxa information
        to distinguish between ambiguous taxa names.

        Args:
            nnme: A nickname to identify the taxon
            taxah: A dictionary containing taxonomic information

        Returns:
            A populated NCBITaxon instance
        """

        if isinstance(taxah, dict) is False:
            raise TypeError("Taxa hierarchy should be a dictionary")
        elif all(isinstance(x, str) for x in taxah.keys()) is False:
            raise ValueError("Not all taxa dictionary keys are strings")
        elif all(isinstance(x, str) for x in taxah.values()) is False:
            raise ValueError("Not all taxa dictionary values are strings")

        # Warn user if no backbone ranks have been provided
        if bool(set(taxah.keys()) & set(BACKBONE_RANKS_EX)) is False:
            LOGGER.warning(f"No backbone ranks provided in {nnme}'s taxa hierarchy")

        # Find last dictionary key
        f_key = list(taxah.keys())[-1]

        # Then find corresponding entry as a search term
        s_term = taxah[f_key]

        # Search for name in the local names database
        sql = "select * from names where name_txt = '%s'" % s_term
        taxon_rows = self.ncbi_conn.execute(sql).fetchall()

        # Store count of the number of rows found
        c = len(taxon_rows)

        # Check that a singular record has been provided before proceeding
        if c == 1:
            # Find taxa ID as single entry in the list
            taxon_row = self.ncbi_conn.execute(sql).fetchone()
            tID = taxon_row["tax_id"]
            # Use ID lookup function to generate as a NCBITaxon object
            mtaxon = self.id_lookup(nnme, tID)

        # Catch cases where no record is found
        elif c == 0:
            # Check if there actually is any higher taxonomy provided
            if len(taxah.keys()) == 1:
                LOGGER.error(
                    f"Taxa {nnme} cannot be found and its higher "
                    f"taxonomic hierarchy is absent"
                )
                return gen_invalid_NCBITaxon()

            # If there is then set up a loop over it
            fnshd = False
            cnt = 1
            while fnshd is False:
                # Increment counter
                cnt += 1
                # Use to find higher taxonomic level to use in search
                f_key = list(taxah.keys())[-cnt]
                new_s_term = taxah[f_key]

                # Search for name in the local names database
                sql = "select * from names where name_txt = '%s'" % new_s_term
                taxon_row = self.ncbi_conn.execute(sql).fetchall()

                # Process the response
                c = len(taxon_row)

                # Check how many records have been found
                if c == 1:
                    # Find taxa ID as single entry in the list
                    taxon_row = self.ncbi_conn.execute(sql).fetchone()
                    tID = taxon_row["tax_id"]
                    # Use ID lookup function to generate as a NCBITaxon object
                    mtaxon = self.id_lookup(nnme, tID)
                    # Store originally supplied rank
                    mtaxon.orig = list(taxah.keys())[-1]
                    # Warn the user that a higher taxonomic rank is being used
                    LOGGER.warning(
                        f"{s_term} not registered with NCBI, but "
                        f"higher level taxon {new_s_term} is"
                    )
                    fnshd = True
                elif c > 1:  # Not going to handle ambiguities in this case
                    LOGGER.error(
                        f"Taxa {nnme} cannot be found and its higher "
                        f"taxonomic hierarchy is ambiguous"
                    )
                    return gen_invalid_NCBITaxon()
                # Catch when all the provided hierarchy has been exhausted
                elif cnt == len(taxah.keys()):
                    fnshd = True

                # If valid higher taxon not found print error
                if "mtaxon" not in locals():
                    LOGGER.error(
                        f"Taxa {nnme} cannot be found and neither can "
                        f"its higher taxonomic hierarchy"
                    )
                    return gen_invalid_NCBITaxon()

        # Case where only one rank has been provided
        elif len(taxah) == 1:
            # Preallocate container for the ranks
            t_ranks = []
            # First check if multiple taxa have the same rank
            for i in range(c):
                # Find taxa ID as single entry in the list
                tID = taxon_rows[i]["tax_id"]
                # Use ID lookup function to generate a temporary taxon
                temp_taxon = self.id_lookup(nnme, tID)
                # Add rank to list
                t_ranks.append(list(temp_taxon.taxa_hier.keys())[-1])

            # Record whether ranks match expected rank
            if f_key == "kingdom":
                mtch = [
                    True if x == f_key or x == "superkingdom" else False
                    for x in t_ranks
                ]
            else:
                mtch = [True if x == f_key else False for x in t_ranks]

            # If there's only one match, we have determined the correct entry
            if sum(mtch) == 1:
                # Find relevant index
                ind = ([i for i, x in enumerate(mtch) if x])[0]
                # Find taxa ID as single entry in the list
                tID = taxon_rows[ind]["tax_id"]
                # Use ID lookup function to generate as a NCBITaxon object
                mtaxon = self.id_lookup(nnme, tID)
            # Then check whether multiple taxonomic levels have been provided
            else:
                # If not raise an error
                LOGGER.error(
                    f"Taxa {nnme} cannot be found using only one "
                    f"taxonomic level, more should be provided"
                )
                return gen_invalid_NCBITaxon()
        # Higher ranks provided
        else:
            # Find second from last dictionary key
            f_key = list(taxah.keys())[-2]
            # Then find corresponding entry as a search term
            s_term = taxah[f_key]

            # Search for name in the local names database
            sql = "select * from names where name_txt = '%s'" % s_term
            p_taxon_rows = self.ncbi_conn.execute(sql).fetchall()

            # Store count of the number of records found
            pc = len(p_taxon_rows)

            # Check that single parent taxa exists in records
            if pc == 0:
                LOGGER.error(f"Provided parent taxa for {nnme} not found")
                return gen_invalid_NCBITaxon()
            elif pc > 1:
                LOGGER.error(f"More than one possible parent taxa for {nnme} found")
                return gen_invalid_NCBITaxon()
            else:
                # Find parent taxa ID as single entry in the list
                p_taxon_row = self.ncbi_conn.execute(sql).fetchone()
                pID = p_taxon_row["tax_id"]
                # Then use ID lookup function to find generate as a NCBITaxon object
                ptaxon = self.id_lookup("parent", pID)

            # Save parent taxa rank and name
            p_key = list(ptaxon.taxa_hier.keys())[-1]
            p_val = ptaxon.taxa_hier[p_key]

            # Use list comprehension to make list of potential taxa
            potents = [
                self.id_lookup(f"c{i}", taxon_rows[i]["tax_id"]) for i in range(0, c)
            ]

            # Check if relevant rank exists in the child taxa
            child = [p_key in potents[i].taxa_hier for i in range(0, c)]

            # Then look for matching rank name pairs
            for i in range(0, c):
                # Only check cases where rank exists in the child taxa
                if child[i] is True:
                    child[i] = potents[i].taxa_hier[f"{p_key}"] == p_val

            # Check for errors relating to finding too many or few child taxa
            if sum(child) == 0:
                LOGGER.error(f"Parent taxa not actually a valid parent of {nnme}")
                return gen_invalid_NCBITaxon()
            elif sum(child) > 1:
                LOGGER.error(
                    f"Parent taxa for {nnme} refers to multiple " f"possible child taxa"
                )
                return gen_invalid_NCBITaxon()
            else:
                # Find index corresponding to correct child taxa
                tID = int(taxon_rows[child.index(True)]["tax_id"])
                # Use ID lookup function to find generate as a NCBITaxon object
                mtaxon = self.id_lookup(nnme, tID)

        # Check whether mtaxon has actually been populated by the id_lookup
        if mtaxon.name == "INVALID":
            return mtaxon

        # Find last dictionary key
        f_key = list(mtaxon.taxa_hier.keys())[-1]

        # Check if taxonomic rank supplied is used
        if (
            mtaxon.orig is None
            and f_key != list(taxah.keys())[-1]
            and f_key != "superkingdom"
        ):
            # If not raise an error
            LOGGER.error(
                f"{list(taxah.values())[-1]} is a {f_key}"
                f" not a {list(taxah.keys())[-1]}"
            )
            return gen_invalid_NCBITaxon()
        elif list(taxah.keys())[-1] == "kingdom" and f_key == "superkingdom":
            # If not print a warning
            LOGGER.warning(
                f"NCBI records {(mtaxon.taxa_hier[f_key])[0]} as "
                f"a superkingdom rather than a kingdom"
            )
        # Then check whether originally supplied name is still used
        elif mtaxon.orig is None and taxah[f_key] != (mtaxon.taxa_hier[f_key])[0]:
            # If not print a warning
            LOGGER.warning(
                f"{taxah[f_key]} not accepted usage should be "
                f"{(mtaxon.taxa_hier[f_key])[0]} instead"
            )
            # And record the fact that usage is superseded
            mtaxon.superseed = True

        return mtaxon


class GBIFTaxa:
    """Manage a set of GBIF taxon data and associated GBIFTaxon instances.

    A class to hold a list of taxon names and a validated taxonomic
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
            gbif_status (str)]

        Where a taxon is not accepted or doubtful on GBIF, two entries are
        inserted for the taxon, one under the canon name and one under the
        provided name. They will share the same worksheet name and so can
        be paired back up for description generation. The worksheet name
        for parent taxa and deeper taxonomic hierarchy is set to None.

    The index_higher_taxa method can be used to extend the taxon_index to
    include all of the higher taxa linking the validated taxa.

    The index can then be used:

    a) to generate the taxonomic coverage section of the dataset description, and
    b) to populate a database table to index the taxonomic coverage of datasets.

    Args:
        resources: A Resources instance.

    Attributes:
        taxon_index: A list containing taxon index lists
        taxon_names: A set of worksheet names
        parents: A dictionary linking tuples of taxon parent information to
            GBIFTaxon instances
        hierarchy: A set of lists containing the complete taxonomic hierarchy for taxa
            in the GBIFTaxa instance.
        n_errors: A count of processing errors when loading and validating taxa
        taxon_names_used: A set used to track which taxon names have been used in data
            worksheets
    """

    def __init__(self, resources: Resources) -> None:
        self.taxon_index: list[list] = []
        self.taxon_names: set[str] = set()
        self.parents: dict[tuple, GBIFTaxon] = dict()
        self.hierarchy: set[list] = set()
        self.n_errors: int = 0
        self.taxon_names_used: set[str] = set()

        # Get the validator instance
        self.validator = GBIFValidator(resources)

    @loggerinfo_push_pop("Loading GBIFTaxa worksheet")
    def load(self, worksheet: worksheet) -> None:
        """Populate a GBIFTaxa instance from an Excel worksheet.

        This method loads a set of taxa from the rows of a SAFE formatted GBIFTaxa
        worksheet and populates the taxonomic hierarchy for those rows. The GBIFTaxa
        instance is updated.

        Args:
            worksheet: An openpyxl worksheet instance using the GBIFTaxa formatting
        """

        start_errors = COUNTER_HANDLER.counters["ERROR"]

        # Get the data read in.
        LOGGER.info("Reading taxa data")
        FORMATTER.push()
        dframe = GetDataFrame(worksheet)

        if not dframe.data_columns:
            LOGGER.error("No data or only headers in GBIFTaxa worksheet")
            FORMATTER.pop()
            return

        # Dupe headers likely cause serious issues, so stop
        if "duplicated" in dframe.bad_headers:
            LOGGER.error("Cannot parse taxa with duplicated headers")
            FORMATTER.pop()
            return

        # Get the headers
        headers = IsLower(dframe.headers).values

        # Field cleaning
        core_fields = {"name", "taxon name", "taxon type"}
        missing_core = core_fields.difference(headers)

        if missing_core:
            # core names are not found so can't continue
            LOGGER.error("Missing core fields: ", extra={"join": missing_core})
            FORMATTER.pop()
            return

        # TODO - Test this new behaviour
        # Fields used to describe taxa (not including comments)
        tx_fields = [
            "name",
            "taxon name",
            "taxon type",
            "taxon id",
            "ignore id",
            "parent name",
            "parent type",
            "parent id",
        ]
        all_fields = set(tx_fields + ["comments"])

        # Now check that there are no unexpected (i.e. likely misspelled) fields
        unexpected_headers = set(headers).difference(all_fields)

        if unexpected_headers:
            # An unexpected header (which might well be misspelled) was found
            LOGGER.error(
                "Unexpected (or misspelled) headers found:",
                extra={"join": unexpected_headers},
            )
            FORMATTER.pop()
            return

        # Any duplication in names
        dupl_taxon_names = HasDuplicates([dframe.data_columns[headers.index("name")]])

        if dupl_taxon_names:
            LOGGER.error(
                "Duplicated names found: ", extra={"join": dupl_taxon_names.duplicated}
            )

        # get dictionaries of the taxa
        taxa = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]
        FORMATTER.pop()

        # check number of taxa found
        if len(taxa) == 0:
            LOGGER.info("No taxon rows found")
            return

        # Standardise to the expected fields, filling in None for any
        # completely missing fields (parent fields could be missing).
        taxa = [{fld: tx.get(fld) for fld in tx_fields} for tx in taxa]

        # Standardize the taxon representation into lists of taxon and parent data
        # Note that parent tuples cannot have an ignore id.
        #     [name,
        #       [taxon name, taxon type, taxon id, ignore id],
        #       [parent name, parent type, parent id]]

        for idx, row in enumerate(taxa):
            # Standardise blank values to None
            row = {ky: None if blank_value(vl) else vl for ky, vl in row.items()}
            taxon_info = [
                row["taxon name"],
                row["taxon type"],
                row["taxon id"],
                row["ignore id"],
            ]
            parent_info: Optional[list] = [
                row["parent name"],
                row["parent type"],
                row["parent id"],
            ]

            # If there is no parent information, replace the parent tuple with None
            if parent_info == [None, None, None]:
                parent_info = None

            self.taxon_names.update([row["name"]])
            LOGGER.info(f"Validating row {idx + 1}: {row['name']}")
            FORMATTER.push()
            self.validate_and_add_taxon((row["name"], taxon_info, parent_info))
            FORMATTER.pop()

        # Add the higher taxa
        self.index_higher_taxa()

        # summary of processing
        self.n_errors = COUNTER_HANDLER.counters["ERROR"] - start_errors
        if self.n_errors is None:
            LOGGER.critical("GBIFTaxa error logging has broken!")
        elif self.n_errors > 0:
            LOGGER.info("GBIFTaxa contains {} errors".format(self.n_errors))
        else:
            LOGGER.info("{} taxa loaded correctly".format(len(self.taxon_names)))

        FORMATTER.pop()

    # TODO - would be nice to use the decorator, but more complex than
    #        I had anticipated: https://stackoverflow.com/questions/11731136/
    #        Could do this via e.g. @loggerinfo_push_pop(f'Validating
    #        {self._row_description}') but this implementation ties
    #        validate_and_add_taxon() to needing that property populated

    def validate_and_add_taxon(self, taxon_input: tuple) -> None:
        """Add a GBIF formatted taxon row to the GBIFTaxa instance.

        This method takes user information on a taxon, and optionally a parent taxon,
        validates it and updates the GBIFTaxa instance to include the new details.

        This is typically used to process rows found in a dataset with a GBIFTaxa
        formatted table, can also be used to populate a GBIFTaxa instance
        programmatically.

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
        """

        m_name, taxon_info, parent_info = taxon_input

        # Sanitise worksheet names for taxa - only keep unpadded strings.
        if m_name is None or not isinstance(m_name, str) or m_name.isspace():
            LOGGER.error("Worksheet name missing, whitespace only or not text")
        elif m_name != m_name.strip():
            LOGGER.error(f"Worksheet name has whitespace padding: {repr(m_name)}")
            m_name = m_name.strip()
            self.taxon_names.add(m_name)
        else:
            self.taxon_names.add(m_name)

        # Check the parent details
        p_fail = False
        if parent_info is not None:
            # Name and rank must be unpadded strings - can still check cleaned padded
            # strings
            for idx, idx_name in ((0, "Parent name"), (1, "Parent rank")):
                val = parent_info[idx]

                if val is None or not isinstance(val, str):
                    LOGGER.error(f"{idx_name} missing or not text")
                    p_fail = True
                elif val != val.strip():
                    LOGGER.error(f"{idx_name} has whitespace padding: {repr(val)}")
                    parent_info[idx] = val.strip()

            # ID can be None or an integer (openpyxl loads all values as float)
            if not (
                parent_info[2] is None
                or (isinstance(parent_info[2], float) and parent_info[2].is_integer())
                or isinstance(parent_info[2], int)
            ):
                LOGGER.error("Parent GBIF ID contains value that is not an integer")
                p_fail = True

        # Check the main taxon details
        mfail = False

        # Name and rank must be unpadded strings - can still check cleaned padded
        # strings
        for idx, idx_name in ((0, "Taxon name"), (1, "Taxon rank")):
            val = taxon_info[idx]

            if val is None or not isinstance(val, str) or val.isspace():
                LOGGER.error(f"{idx_name} missing, whitespace only or not text")
                mfail = True
            elif val != val.strip():
                LOGGER.error(f"{idx_name} has whitespace padding: {repr(val)}")
                taxon_info[idx] = val.strip()

        # GBIF ID and Ignore ID can be None or an integer (openpyxl loads all values as
        # float)
        for idx, idx_name in ((2, "GBIF ID"), (3, "Ignore ID")):
            val = taxon_info[idx]

            if not (
                val is None
                or (isinstance(val, float) and val.is_integer())
                or isinstance(val, int)
            ):
                LOGGER.error(f"{idx_name} contains value that is not an integer: {val}")
                mfail = True

        if p_fail:
            LOGGER.error("Parent taxon details not properly formatted, cannot validate")

        if mfail:
            LOGGER.error("Taxon details not properly formatted, cannot validate")

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
            p_taxon = GBIFTaxon(
                name=parent_info[0], rank=parent_info[1], gbif_id=parent_info[2]
            )

            # Look for a match
            if p_taxon.is_backbone:
                p_taxon = self.validator.search(p_taxon)

                # Update the hierarchy and index with the search results
                self.hierarchy.update(
                    [rw for rw in p_taxon.hierarchy if rw[1] is not None]
                )

                # For non-canon taxa, add the canon hierarchy and the non-canon usage as
                # this is not included in the hierachy (starts with the canon usage)
                if (
                    p_taxon.is_backbone
                    and p_taxon.found
                    and not p_taxon.is_canon
                    and p_taxon.canon_usage
                ):
                    self.hierarchy.update(
                        [
                            rw
                            for rw in p_taxon.canon_usage.hierarchy
                            if rw[1] is not None
                        ]
                        + [(p_taxon.rank, p_taxon.gbif_id)]
                    )

            # Store the parent taxon keyed by parent information (needs tuple)
            self.parents[tuple(parent_info)] = p_taxon

        # Report on the parent information
        if p_taxon is not None:
            if not p_taxon.is_backbone:
                LOGGER.error(f"Parent taxon ({p_taxon.name}) is not of a backbone rank")

            elif not p_taxon.found:
                LOGGER.error(f"Parent taxon ({p_taxon.name}) {p_taxon.lookup_status}")

            elif not p_taxon.is_canon and p_taxon.canon_usage:
                LOGGER.warning(
                    f"Parent taxon ({p_taxon.name}) considered a {p_taxon.taxon_status}"
                    f" of {p_taxon.canon_usage.name} in GBIF backbone"
                )
            else:
                LOGGER.info(f"Parent taxon ({p_taxon.name}) accepted")

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
        m_taxon = GBIFTaxon(
            name=taxon_info[0], rank=taxon_info[1], gbif_id=taxon_info[2]
        )
        ignore_gbif = taxon_info[3]

        if ignore_gbif is not None:
            # Handle ignored matches first

            # The taxon must be a backbone taxon - can't ignore impossible matches
            if not m_taxon.is_backbone:
                LOGGER.error(
                    "Ignore ID can only be used with GBIF backbone taxon ranks"
                )
            else:
                # It should also be found and the ignore ID should match to the actual
                # usage or canon usage.
                m_taxon = self.validator.search(m_taxon)

                if not m_taxon.found:
                    LOGGER.error("Taxon with Ignore ID not found in GBIF backbone")
                elif m_taxon.is_canon and (m_taxon.gbif_id != ignore_gbif):
                    LOGGER.error(
                        f"Ignore ID does not match the canon GBIF usage ("
                        f"{m_taxon.gbif_id})"
                    )
                elif (
                    not m_taxon.is_canon
                    and m_taxon.canon_usage
                    and (m_taxon.canon_usage.gbif_id != ignore_gbif)
                ):
                    LOGGER.error(
                        f"Taxon is non-canon and Ignore ID does not match the canon "
                        f"GBIF usage ({m_taxon.canon_usage.gbif_id})"
                    )
                else:
                    LOGGER.info("Canon GBIF usage ignored")

            # It must also have a valid parent.
            if p_taxon is None:
                LOGGER.error("Taxa with Ignore ID must provide parent information.")
            elif not p_taxon.found:
                LOGGER.error("Taxon with Ignore ID has invalid parent information.")
            else:
                LOGGER.info(
                    "Taxon with ignored canon usage has valid parent information."
                )
                # Update index - no taxon hierarchy except for parent - link to the
                # canon usage of parent
                if p_taxon.is_canon:
                    pid = p_taxon.gbif_id
                elif p_taxon.canon_usage is not None:
                    pid = p_taxon.canon_usage.gbif_id

                self.taxon_index.append(
                    [m_name, -1, pid, m_taxon.name, m_taxon.rank, "user"]
                )

        elif not m_taxon.is_backbone:
            # Now handle non-backbone cases - just needs a valid parent.
            if p_taxon is None:
                LOGGER.error(
                    f"Taxon of type {m_taxon.rank} must provide parent information."
                )
            elif not p_taxon.found:
                # Non backbone with bad parent information
                LOGGER.error(
                    f"Taxon of type {m_taxon.rank} has invalid parent information."
                )
            else:
                # Non backbone with with good parent info
                LOGGER.info(
                    f"Taxon of type {m_taxon.rank} has valid parent information"
                )
                # Update index - no taxon hierarchy except for parent - link to the
                # canon usage of parent
                if p_taxon.is_canon:
                    pid = p_taxon.gbif_id
                elif p_taxon.canon_usage is not None:
                    pid = p_taxon.canon_usage.gbif_id

                self.taxon_index.append(
                    [m_name, -1, pid, m_taxon.name, m_taxon.rank, "user"]
                )

        else:
            # Otherwise try and validate backbone taxon
            m_taxon = self.validator.search(m_taxon)

            if m_taxon.found:
                # Add the index entry and update hierarchy
                self.taxon_index.append(
                    [
                        m_name,
                        m_taxon.gbif_id,
                        m_taxon.parent_id,
                        m_taxon.name,
                        m_taxon.rank,
                        m_taxon.taxon_status,
                    ]
                )

                # default to named taxon hierarchy
                hier_to_use = m_taxon.hierarchy

                # Now check for non-canon usage, adding the canon taxa and switching the
                # taxon hierarchy to step up above the canon usage:
                #    (non-canon -> canon -> rest of hierarchy)
                if m_taxon.canon_usage is not None:
                    # Add the canon index entry under the worksheet name and update to
                    # using the canon hierarchy
                    self.taxon_index.append(
                        [
                            m_name,
                            m_taxon.canon_usage.gbif_id,
                            m_taxon.canon_usage.parent_id,
                            m_taxon.canon_usage.name,
                            m_taxon.canon_usage.rank,
                            m_taxon.canon_usage.taxon_status,
                        ]
                    )

                    hier_to_use = m_taxon.canon_usage.hierarchy

                # Now update the hierarchy
                self.hierarchy.update([rw for rw in hier_to_use if rw[1] is not None])

                # Reporting
                if p_taxon is None:
                    if m_taxon.is_canon:
                        LOGGER.info(
                            f"Taxon found in GBIF backbone ({m_taxon.taxon_status})"
                        )
                    elif m_taxon.canon_usage is not None:
                        LOGGER.warning(
                            f"Taxon considered a {m_taxon.taxon_status} "
                            f"of {m_taxon.canon_usage.name} in GBIF backbone"
                        )

                elif p_taxon is not None:
                    if p_taxon.found:
                        # Good backbone with good parent - are they compatible? Check if
                        # all entries in the parent hierarchy appear in the taxon
                        # hierarchy
                        if not set(p_taxon.hierarchy).issubset(m_taxon.hierarchy):
                            LOGGER.error(
                                f"Taxon in GBIF backbone ({m_taxon.taxon_status}) with "
                                f"incompatible parent information"
                            )
                        else:
                            LOGGER.info(
                                f"Taxon in GBIF backbone ({m_taxon.taxon_status}) with "
                                f"compatible parent information"
                            )

                    else:
                        # Good backbone with bad parent
                        LOGGER.error(
                            f"Taxon in GBIF backbone ({m_taxon.taxon_status}) but with "
                            f"invalid parent information."
                        )

            elif not m_taxon.found:
                if p_taxon is None:
                    # Taxon is a backbone type but is not found in GBIF and has no
                    # parent info
                    if m_taxon.lookup_status == "No match found":
                        LOGGER.error("Taxon name and rank combination not found")
                    else:
                        LOGGER.error(f"GBIF issue: {m_taxon.lookup_status}")

                elif not p_taxon.found:
                    # Taxon is a backbone type but not found and parent not found either
                    LOGGER.error(
                        "Taxon not found in GBIF and has invalid parent information."
                    )
                else:
                    # Taxon is a backbone type but not found but does have valid parent
                    # info
                    LOGGER.info(
                        "Taxon not found in GBIF but has valid parent information"
                    )

                    # Add to index  - parent already in hierarchy so nothing to add
                    self.taxon_index.append(
                        [
                            m_name,
                            -1,
                            p_taxon.gbif_id,
                            m_taxon.name,
                            m_taxon.rank,
                            "user",
                        ]
                    )

    @loggerinfo_push_pop("Indexing taxonomic hierarchy")
    def index_higher_taxa(self) -> None:
        """Extend the taxon index to include higher taxa.

        This method uses the taxon hierarchy entries to add higher taxa to the taxon
        index for the instance. It does not duplicate taxa already in the index.
        """

        known = [tx[1] for tx in self.taxon_index if tx[1] != -1]
        to_add = [tx for tx in self.hierarchy if tx[1] not in known]
        to_add.sort(key=lambda val: BACKBONE_RANKS.index(val[0]))

        # Look up the taxonomic hierarchy
        for tx_lev, tx_id in to_add:
            higher_taxon = self.validator.id_lookup(tx_id)
            self.taxon_index.append(
                [
                    None,
                    higher_taxon.gbif_id,
                    higher_taxon.parent_id,
                    higher_taxon.name,
                    higher_taxon.rank,
                    higher_taxon.taxon_status,
                ]
            )
            LOGGER.info(f"Added {tx_lev} {higher_taxon}")

    @property
    def is_empty(self) -> bool:
        """Check if a GBIFTaxa instance contains any taxa."""
        return len(self.taxon_names) == 0


class NCBITaxa:
    """Manage a set of NCBI taxon data and associated NCBITaxon instances.

    A class to hold a list of taxon names and a validated taxonomic
    index for those taxa and their taxonomic hierarchy. The validate_taxon
    method checks that the provided taxon hierarchy and (optional) NCBI ID can be
    validated against the NCBI database (and both refer to the same taxon).

    i)  the taxon_names attribute of the dataset, which is just a set of
        names used as a validation list for taxon names used in data worksheets.
    ii) the taxon_index attribute of the dataset, which contains a set
        of lists structured as:

            [worksheet_name (str),
            ncbi_id (int),
            ncbi_parent_id (int),
            canonical_name (str),
            taxonomic_rank (str),
            ncbi_status (str)]

        Where a taxon is not accepted or doubtful on GBIF, two entries are
        inserted for the taxon, one under the canon name and one under the
        provided name. They will share the same worksheet name and so can
        be paired back up for description generation. The worksheet name
        for parent taxa and deeper taxonomic hierarchy is set to None.

    The index_higher_taxa method can be used to extend the taxon_index to
    include all of the higher taxa linking the validated taxa.

    The index can then be used:

    a) to generate the taxonomic coverage section of the dataset description, and
    b) to populate a database table to index the taxonomic coverage of datasets.

    Args:
        resources: A Resources instance.

    Attributes:
        taxon_index: A list containing taxon index lists
        taxon_names: A set of worksheet names
        parents: A dictionary linking tuples of taxon parent information to
            GBIFTaxon instances
        hierarchy: A set of lists containing the complete taxonomic hierarchy for taxa
            in the GBIFTaxa instance.
        n_errors: A count of processing errors when loading and validating taxa
        taxon_names_used: A set used to track which taxon names have been used in data
            worksheets
    """

    def __init__(self, resources: Resources) -> None:
        self.taxon_index: list[list] = []
        self.taxon_names: set[str] = set()
        self.hierarchy: set[tuple] = set()
        self.n_errors: int = 0

        # Get the validator instance
        self.validator = NCBIValidator(resources)

    @loggerinfo_push_pop("Loading NCBITaxa worksheet")
    def load(self, worksheet: worksheet) -> None:
        """Populate an NCBITaxa instance from an Excel worksheet.

        This method loads a set of taxa from the rows of a SAFE formatted NCBITaxa
        worksheet and populates the taxonomic hierarchy for those rows. The GBIFTaxa
        instance is updated.

        Args:
            worksheet: An openpyxl worksheet instance using the GBIFTaxa formatting
        """

        start_errors = COUNTER_HANDLER.counters["ERROR"]

        # Get the data read in.
        LOGGER.info("Reading NCBI taxa data")
        FORMATTER.push()
        dframe = GetDataFrame(worksheet)

        if not dframe.data_columns:
            LOGGER.error("No data or only headers in Taxa worksheet")
            FORMATTER.pop()
            return

        # Dupe headers likely cause serious issues, so stop
        if "duplicated" in dframe.bad_headers:
            LOGGER.error("Cannot parse taxa with duplicated headers")
            FORMATTER.pop()
            return

        # Get the headers
        headers = IsLower(dframe.headers).values

        # Field cleaning
        core_fields = {"name", "ncbi id"}
        missing_core = core_fields.difference(headers)

        if missing_core:
            # core names are not found so can't continue
            LOGGER.error("Missing core fields: ", extra={"join": missing_core})
            FORMATTER.pop()
            return

        # Possible fields in this case are the two core fields + comments + all
        # taxonomic ranks defined in NCBI
        all_fields = set(
            list(core_fields) + BACKBONE_RANKS_EX + EXTRA_NCBI_RANKS + ["comments"]
        )

        # Now check that there are no unexpected (i.e. likely misspelled) fields
        unexpected_headers = set(headers).difference(all_fields)

        if unexpected_headers:
            # An unexpected header (which might well be misspelled) was found
            LOGGER.error(
                "Unexpected (or misspelled) headers found:",
                extra={"join": unexpected_headers},
            )
            FORMATTER.pop()
            return

        # Check that at least two backbone taxa have been provided
        fnd_rnks = set(BACKBONE_RANKS_EX).intersection(headers)

        if len(fnd_rnks) < 2:
            # can't continue if less than two backbone ranks are provided
            LOGGER.error("Less than two backbone taxonomic ranks are provided")
            FORMATTER.pop()
            return

        # Find core field indices and use to isolate non core (i.e. taxonomic) headers
        core_inds = [headers.index(item) for item in core_fields]
        non_core = [element for i, element in enumerate(headers) if i not in core_inds]

        # Check to see if comments is provided
        if "comments" in non_core:
            # Check that it is the last header
            if non_core[-1] != "comments":
                LOGGER.error(
                    "If 'Comments' is provided as a field it must be the last column"
                )
                FORMATTER.pop()
                return
            else:
                # If it's the last header go ahead and delete it
                del non_core[-1]

        # Check that backbone ranks are in the correct order
        l1 = [v for v in BACKBONE_RANKS_EX if v in fnd_rnks]
        l2 = [v for v in non_core if v in fnd_rnks]

        if l1 != l2:
            LOGGER.error("Backbone taxonomic ranks not provided in the correct order")

        # Check that if subspecies has been provided species is provided
        if "subspecies" in headers:
            if "species" not in headers:
                LOGGER.error("If subspecies is provided so must species")
                FORMATTER.pop()
                return

        # Same check
        if "species" in headers:
            if "genus" not in headers:
                LOGGER.error("If species is provided so must genus")
                FORMATTER.pop()
                return

        # Any duplication in names
        dupl_taxon_names = HasDuplicates(dframe.data_columns[headers.index("name")])

        if dupl_taxon_names:
            LOGGER.error(
                "Duplicated names found: ", extra={"join": dupl_taxon_names.duplicated}
            )

        # get dictionaries of the taxa
        taxa = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]
        FORMATTER.pop()

        # check number of taxa found
        if len(taxa) == 0:
            LOGGER.info("No taxon rows found")
            return

        # Change data format such that it is appropriate to be used to search
        # the NCBI database, i.e. convert from k__ notation, generate species
        # binomials and subspecies trinomials, this is then all collected in a
        # dictionary.

        for idx, row in enumerate(taxa):
            # Strip k__ notation so that the text is appropriate for the search
            for rnk in non_core:
                row[rnk], match = taxa_strip(row[rnk], rnk)
                if match is False:
                    print(
                        f"Implied rank of {row[rnk]} in row {idx + 1} does not"
                        f" match rank it is assigned"
                    )

            # Standardise blank values to None
            row = {ky: None if blank_value(vl) else vl for ky, vl in row.items()}
            # Replace any NA values with None
            row = {ky: None if vl == "NA" else vl for ky, vl in row.items()}

            # Start with empty dictionary for taxonomic hierarchy
            taxa_hier = {}

            # Loop over all ranks to populate the dictionary
            for rnk in non_core:
                if row[rnk] is not None:
                    if rnk == "species":
                        taxa_hier[rnk] = construct_bi_or_tri(
                            row["genus"], row[rnk], False
                        )
                    elif rnk == "subspecies":
                        taxa_hier[rnk] = construct_bi_or_tri(
                            row["species"], row[rnk], True
                        )
                    else:
                        taxa_hier[rnk] = row[rnk]

            self.taxon_names.update([row["name"]])
            LOGGER.info(f"Validating row {idx + 1}: {row['name']}")
            FORMATTER.push()
            self.validate_and_add_taxon((row["name"], taxa_hier, row["ncbi id"]))
            FORMATTER.pop()

        # Add the higher taxa
        self.index_higher_taxa()

        # summary of processing
        self.n_errors = COUNTER_HANDLER.counters["ERROR"] - start_errors
        if self.n_errors is None:
            LOGGER.critical("NCBITaxa error logging has broken!")
        elif self.n_errors > 0:
            LOGGER.info("NCBITaxa contains {} errors".format(self.n_errors))
        else:
            LOGGER.info("{} taxa loaded correctly".format(len(self.taxon_names)))

        FORMATTER.pop()

    def validate_and_add_taxon(self, ncbi_taxon_input: tuple) -> None:
        """Add a GBIF formatted taxon row to the GBIFTaxa instance.

        This method takes user information on a taxon, and optionally an NCBI taxonomy
        ID, validates it against the NCBI database and updates the NCBITaxa instance to
        include the new details.

        This is typically used to process rows found in a dataset with an NCBITaxa
        formatted table, can also be used to populate a NCBITaxa instance
        programmatically.

        The taxon_input has the form:

            ['m_name', 'taxon_hier', 'ncbi_id']

        If there is no NCBI ID, the structure is:

            ['m_name', 'taxon_hier', None]

        Args:
            ncbi_taxon_input: NCBITaxon information in standard form as above
        """

        m_name, taxon_hier, ncbi_id = ncbi_taxon_input

        # Sanitise worksheet names for taxa - only keep unpadded strings.
        if (
            m_name is None
            or not isinstance(m_name, str)
            or m_name.isspace()
            or not m_name
        ):
            LOGGER.error("Worksheet name missing, whitespace only or not text")
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
        if not (
            ncbi_id is None
            or (isinstance(ncbi_id, float) and ncbi_id.is_integer())
            or isinstance(ncbi_id, int)
        ):
            LOGGER.error("NCBI ID contains value that is not an integer")
            i_fail = True

        # Check the main taxon details
        h_fail = False

        # Check that a dictionary with at least one entry has been provided
        if not isinstance(taxon_hier, dict) or len(taxon_hier) == 0:
            LOGGER.error("Taxa hierarchy should be a (not empty) dictionary")
            h_fail = True
        # Otherwise check for padding of dictionary keys and values
        else:
            # Make a translation table
            translate = {}
            # Loop over all dictionary keys
            for idx in taxon_hier.keys():
                if not isinstance(idx, str):
                    LOGGER.error(f"Non-string dictionary key used: {repr(idx)}")
                    h_fail = True
                elif idx == "" or idx.isspace():
                    LOGGER.error("Empty dictionary key used")
                    h_fail = True
                elif idx != idx.strip():
                    LOGGER.error(f"Dictionary key has whitespace padding: {repr(idx)}")
                    # Save keys to swap to new translation table
                    translate[idx] = idx.strip()

                # Extract corresponding dictionary value
                val = taxon_hier[idx]
                # Then perform similar checks on dictionary value
                if not isinstance(val, str):
                    LOGGER.error(f"Non-string dictionary value used: {repr(val)}")
                    h_fail = True
                elif val == "" or val.isspace():
                    LOGGER.error("Empty dictionary value used")
                    h_fail = True
                elif val != val.strip():
                    LOGGER.error(
                        f"Dictionary value has whitespace padding: {repr(val)}"
                    )
                    taxon_hier[idx] = val.strip()

            # Use translation table to replace white-spaced dictionary keys
            for old, new in translate.items():
                taxon_hier[new] = taxon_hier.pop(old)

        # Now check that the taxa hierarchy is correctly ordered
        if i_fail:
            LOGGER.error("Improper NCBI ID provided, cannot be validated")

        if h_fail:
            LOGGER.error("Taxon details not properly formatted, cannot validate")

        if h_fail or i_fail:
            return

        # Now check that dictionary containing taxa hierarchy is properly ordered
        if len(taxon_hier) > 1:  # Only matters if it contains multiple entries
            # Find all keys in taxa hierarchy
            t_ord = list(taxon_hier.keys())
            # Find corresponding keys from backbone
            b_ord = [x for x in BACKBONE_RANKS_EX if x in t_ord]
            # Remove non-backbone keys from taxa order
            t_ord = [x for x in t_ord if x in b_ord]
            # Then catch cases where orders don't match
            if b_ord != t_ord:
                LOGGER.error("Taxon hierarchy not in correct order")
                return

        # Go straight ahead and search for the taxon
        hr_taxon = self.validator.taxa_search(m_name, taxon_hier)
        # Catch case where errors are returned rather than a taxon
        if hr_taxon.name == "INVALID":
            LOGGER.error("Search based on taxon hierarchy failed")
            return

        # Then check if a genbank ID number has been provided
        if ncbi_id is not None:
            id_taxon = self.validator.id_lookup(m_name, int(ncbi_id))
            # Check if taxonomy hierarchy superseded
            if hr_taxon.superseed is True:
                LOGGER.warning(
                    f"Taxonomic classification superseded for "
                    f"{m_name}, using new taxonomic classification"
                )
            elif id_taxon.superseed is True:
                LOGGER.warning(
                    f"NCBI taxa ID superseded for {m_name}, using new taxa ID"
                )
            elif id_taxon != hr_taxon:
                LOGGER.error(
                    f"The NCBI ID supplied for {m_name} does "
                    f"not match hierarchy: expected {hr_taxon.ncbi_id}"
                    f" got {ncbi_id}"
                )
                return
        else:
            # Warn user that superseded taxonomy won't be used
            if hr_taxon.superseed:
                LOGGER.warning(
                    f"Taxonomic classification superseded for "
                    f"{m_name}, using new taxonomic classification"
                )

        # Check that the hierarchy found matches
        self.compare_hier(m_name, hr_taxon, taxon_hier)

        # Find parent ID
        if len(hr_taxon.taxa_hier) > 1:
            # Find taxon one level up
            f_key = list(hr_taxon.taxa_hier.keys())[-2]
            parent_id = (hr_taxon.taxa_hier[f_key])[1]
        else:
            # Set to None if hierarchy is empty (i.e. top level taxa)
            parent_id = None

        # Catch cases where original taxon has not been found
        if hr_taxon.orig is not None:
            self.taxon_index.append(
                [
                    m_name,
                    -1,
                    hr_taxon.ncbi_id,
                    list(taxon_hier.values())[-1],
                    hr_taxon.orig,
                    "user",
                ]
            )
        else:
            # Then check if taxon is superseded
            if hr_taxon.superseed is True or (
                ncbi_id is not None and id_taxon.superseed is True
            ):
                if ncbi_id is None or id_taxon.superseed is False:
                    superseed_id = hr_taxon.ncbi_id
                    # Find supplied name using lowest found rank
                    f_key = list(hr_taxon.taxa_hier.keys())[-1]
                    superseed_name = taxon_hier[f_key]
                else:
                    # Supplied ID is superseded
                    superseed_id = ncbi_id
                    # Check if supplied name superseded
                    if hr_taxon.superseed is True:
                        f_key = list(hr_taxon.taxa_hier.keys())[-1]
                        superseed_name = taxon_hier[f_key]
                    else:
                        superseed_name = hr_taxon.name

                # Add superseded taxon to the index
                self.taxon_index.append(
                    [
                        m_name,
                        superseed_id,
                        parent_id,
                        superseed_name,
                        hr_taxon.rank,
                        "merged",
                    ]
                )

            # Then (also) add non-superseded taxon info to the index
            self.taxon_index.append(
                [
                    m_name,
                    hr_taxon.ncbi_id,
                    parent_id,
                    hr_taxon.name,
                    hr_taxon.rank,
                    "accepted",
                ]
            )

        self.hierarchy.update(list(hr_taxon.taxa_hier.items()))

        # Check if this has succeeded without warnings or errors
        if hr_taxon.superseed is False and hr_taxon.orig is None:
            # Straight forward when there isn't a genbank id, or previously processed
            if ncbi_id is None:
                # If so inform the user of this
                LOGGER.info(f"Taxon ({m_name}) found in NCBI database")
            # Otherwise need to check for superseded ID's
            elif id_taxon.superseed is False:
                LOGGER.info(f"Taxon ({m_name}) found in NCBI database")
        elif hr_taxon.superseed is False:
            LOGGER.info(f"Higher taxon for ({m_name}) resolved in NCBI")

    @loggerinfo_push_pop("Indexing taxonomic hierarchy")
    def index_higher_taxa(self) -> None:
        """Extend the taxon index to include higher taxa.

        This method uses the taxon hierarchy entries to add higher taxa to the taxon
        index for the instance. It does not duplicate taxa already in the index.
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
            self.taxon_index.append([None, tx_id, p_id, tx_nme, tx_lev, "accepted"])
            LOGGER.info(f"Added {tx_lev} {tx_nme}")

    def compare_hier(self, m_name: str, mtaxon: NCBITaxon, taxon_hier: dict) -> None:
        """Validate provided NCBI taxon hierarchy.

        This method compares the retrieved hierarchy of a taxon with the hierarchy that
        was initially supplied. This function only checks that provided information
        matches, missing levels or entries are not checked for.
        """
        # Not worth checking hierarchy in superseded case
        if mtaxon.superseed is True:
            return

        # Find all common ranks between the hierarchies
        rnks = list(set(taxon_hier.keys()) & set(mtaxon.taxa_hier.keys()))

        # Find if all taxa names match
        t_match = [taxon_hier[r] == (mtaxon.taxa_hier[r])[0] for r in rnks]
        if all(t_match) is False:
            # Find indices of non-matching values
            inds = [i for i in range(len(t_match)) if t_match[i] is False]
            # Then loop over these indices giving appropriate warnings
            for ind in inds:
                LOGGER.warning(
                    f"Hierarchy mismatch for {m_name} its {rnks[ind]}"
                    f" should be {(mtaxon.taxa_hier[rnks[ind]])[0]} not "
                    f"{taxon_hier[rnks[ind]]}"
                )

    @property
    def is_empty(self) -> bool:
        """Check if an NCBITaxaWe do instance contains any taxa."""
        return len(self.taxon_names) == 0


class Taxa:
    """Manage combined NCBITaxa and GBIFTaxa instances.

    This class wraps parallel instances of GBIFTaxa and NCBITaxa and provides shared
    properties across the two instances.

    Args:
        resources: A Resources instance


    We are interested in checking that no worksheet
    names are reused when both Taxa sheets are provided, that every worksheet
    name is used somewhere in the Data worksheets, and that every taxon name
    used across the Data worksheets is defined in a Taxa worksheet.

    This overarching class stores instances of the two lower level classes
    (GBIFTaxa, NCBITaxa). It can also store (as `taxon_names_used`) the set
    of all names used across the Data worksheets. The property `is_empty` can
    be used to check whether both of the lower level classes are empty, and
    the property `taxon_names` can be used to find the set of all taxon names
    defined in either GBIFTaxa or NCBITaxa. Finally, the property `repeat_names`
    can be used to find if any names are used in both GBIFTaxa and NCBITaxa
    worksheets.
    """

    def __init__(self, resources: Resources):
        self.gbif_taxa = GBIFTaxa(resources)
        self.ncbi_taxa = NCBITaxa(resources)
        self.taxon_names_used: set[str] = set()

    @property
    def is_empty(self) -> bool:
        """Reports if neither GBIF nor NCBI taxa any taxa loaded."""
        return self.gbif_taxa.is_empty and self.ncbi_taxa.is_empty

    @property
    def taxon_names(self) -> set[str]:
        """Provides loaded taxon names from both NCBI and GBIF taxa."""
        return self.gbif_taxa.taxon_names.union(self.ncbi_taxa.taxon_names)

    @property
    def repeat_names(self) -> set[str]:
        """Reports taxon names duplicated between NCBI and GBIF taxa."""
        return self.gbif_taxa.taxon_names.intersection(self.ncbi_taxa.taxon_names)


def taxon_index_to_text(
    taxa: list[dict], html: bool = False, indent_width: int = 4, auth: str = "GBIF"
) -> Union[str, tags.div]:
    """Render a GBIFTaxa instance as text or html.

    This function takes a taxon index and renders the contents into either a text or
    html representation of the taxonomic hierarchy used in the dataset. Taxonomic ranks
    are indented to render a nested hierarchy. The `auth` argument is used to set
    whether the taxa are validated using GBIF or NCBI and this only affects the
    formatting of the names in the representation.

    Args:
        taxa: A list of taxon dictionaries containing the taxa for a dataset.
        html: Render as html or text.
        indent_width: The indentation width to use for successive taxonomic ranks.
        auth: The taxonomic authority that the taxa are taken from.

    Returns:
        Either a HTML or text representation of the taxa tree.
    """

    def _indent(n: int, use_html: bool = html):
        if use_html:
            return raw("&ensp;-&ensp;" * n)
        else:
            return " " * indent_width * (n - 1)

    def _format_name(tx: dict, use_html: bool = html, auth: str = "GBIF"):
        if auth == "GBIF":
            # format the canonical name
            if tx["taxon_rank"] in ["genus", "species", "subspecies"]:
                if use_html:
                    return tags.i(tx["taxon_name"])
                else:
                    return f"_{tx['taxon_name']}_"
            elif tx["taxon_rank"] in ["morphospecies", "functional group"]:
                return f"[{tx['worksheet_name']}, {tx['taxon_rank']}]"
            else:
                return tx["taxon_name"]

        elif auth == "NCBI":
            # format the canonical name
            if tx["taxon_status"] == "user":
                if tx["taxon_rank"] in BACKBONE_RANKS_EX:
                    return f"[{tx['taxon_name']}]"
                else:
                    return (
                        f"[{tx['taxon_name']}]  (non-backbone rank: {tx['taxon_rank']})"
                    )
            else:
                if tx["taxon_rank"] in ["genus", "species", "subspecies"]:
                    if use_html:
                        return tags.i(tx["taxon_name"])
                    else:
                        return f"_{tx['taxon_name']}_"
                elif tx["taxon_rank"] not in BACKBONE_RANKS_EX:
                    return f"{tx['taxon_name']} (non-backbone rank: {tx['taxon_rank']})"
                else:
                    return tx["taxon_name"]
        else:
            raise ValueError(f"Unknown auth value: {auth}")

    # Container type depends on whether or not html output is required
    if html:
        # Container to hold the output
        html_out = tags.div()
    else:
        html_out = StringIO()

    # group by parent taxon, substituting 0 for None
    # secondary order is then alphabetic based on taxon name
    taxa.sort(key=lambda x: (x["parent_id"] or 0, x["taxon_name"]))

    # Preallocate container to store identity of surplus taxa
    surp_tx_ids = []
    # Define keys that would match in unwanted repeated entries
    match_keys = [
        "taxon_id",
        "parent_id",
        "taxon_name",
        "taxon_rank",
        "taxon_status",
    ]

    # Loop over taxa to filter for repeated entries
    for idx, taxon in enumerate(taxa):
        # Identify elements in taxa where all 5 of the desired keys match
        matches = list(
            map(
                lambda x: x == 5,
                [sum([taxon[k] == item[k] for k in match_keys]) for item in taxa],
            )
        )
        if sum(matches) > 1:
            # Generate reduced list of matching taxa
            taxa_mtch = list(compress(taxa, matches))
            ws_names = [item["worksheet_name"] for item in taxa_mtch]
            # Find first non-None worksheet names
            first_nm = next(name for name in ws_names if name is not None)
            # If it doesn't match worksheet name of taxon, add index to be deleted
            if first_nm != taxon["worksheet_name"]:
                surp_tx_ids.append(idx)

    # Delete taxa that are superfluous by index
    for index in sorted(surp_tx_ids, reverse=True):
        del taxa[index]

    # group taxa by their parent id
    grouped = {k: list(v) for k, v in groupby(taxa, lambda x: x["parent_id"])}

    # start the stack with root taxa, which will have None as a parent (kingdoms for
    # GBIF and superkingdoms for NCBI)
    stack = [({"current": grouped[None][0]}, {"next": grouped[None][1:]})]

    while stack:
        # Handle the current top of the stack: format the canonical name
        current = stack[-1][0]["current"]
        canon_name = _format_name(current)

        # Look for a non-None entry in next that shares the same worksheet name
        next_ws_names = [
            tx["worksheet_name"]
            for tx in stack[-1][1]["next"]
            if tx["worksheet_name"] is not None
        ]

        if current["worksheet_name"] in next_ws_names:
            # pop out the matching entry and find which is 'accepted'
            name_pair = stack[-1][1]["next"].pop(
                next_ws_names.index(current["worksheet_name"])
            )
            if current["taxon_status"] == "accepted":
                as_name = _format_name(name_pair)
                as_status = name_pair["taxon_status"]
            else:
                as_name = canon_name
                as_status = current["taxon_status"]
                canon_name = _format_name(name_pair)

            if html:
                html_txt = [
                    _indent(len(stack)),
                    canon_name,
                    " (as ",
                    as_status,
                    ": ",
                    as_name,
                    ")",
                    tags.br(),
                ]
            else:
                txt = (
                    f"{_indent(len(stack))} {canon_name} (as {as_status}: {as_name})\n"
                )
        else:
            if html:
                html_txt = [_indent(len(stack)), canon_name, tags.br()]
            else:
                txt = f"{_indent(len(stack))} {canon_name}\n"

        if html:
            html_out += html_txt
        else:
            html_out.write(txt)

        # Is this taxon a parent for other taxa - if so add that taxon to the top of
        # the stack, otherwise start looking for a next taxon to push onto the stack.
        # If there is none at the top, pop and look down.
        parent_id = current["taxon_id"]
        if parent_id in grouped:
            stack.append(
                ({"current": grouped[parent_id][0]}, {"next": grouped[parent_id][1:]})
            )
        else:
            while stack:
                push = stack.pop()
                if push[1]["next"]:
                    stack.append(
                        ({"current": push[1]["next"][0]}, {"next": push[1]["next"][1:]})
                    )
                    break

    if html:
        return html_out
    else:
        return html_out.getvalue()


def taxa_strip(name: str, rank: str) -> tuple[str, bool]:
    """Strip NCBI style rank prefixes from taxon names.

    This function removes NCBI `k__` type notation from taxa names. It returns a tuple
    containing the stripped name and a boolean indicating if the provided rank is
    consistent with the rank implied by its prefix.

    Args:
        name: An NCBI taxon name with `k__` style rank prefix
        rank: The expected taxonomic rank for the name.
    """
    if name is None:
        return (None, True)
    elif "__" in name:
        # Strip name down
        ind = name.rfind("_")
        s_name = name[ind + 1 :]
        # Check if ranks match
        match = name[0].lower() == rank[0].lower()
        return (s_name, match)
    else:
        return (name, True)


def construct_bi_or_tri(higher_nm: str, lower_nm: str, tri: bool) -> Optional[str]:
    """Generate a species binomial or a subspecies trinomial.

    The NCBI database sometimes includes extra tags in binomials, such as 'candidatus'.
    This function cleans up those names to remove extra tags. It returns the cleaned
    name or None in the event of a parsing error.

    Args:
        higher_nm: The NCBI genus/species name
        lower_nm: The NCBI species/subspecies name
        tri: Are we looking at a subspecies trinomial
    """

    # Determine whether a species binomial or a subspecies trinomial is considered
    if tri:
        n = 3
        H = "Species"
        L = "Subspecies"
        type = "trinomial"
    else:
        n = 2
        H = "Genus"
        L = "Species"
        type = "binomial"

    # First check if lower level name is a single name
    if len(lower_nm.split()) == 1 and len(higher_nm.split()) == n - 1:
        return higher_nm.strip() + " " + lower_nm.strip()
    # Look for Candidatus formats
    elif lower_nm.lower().startswith("candidatus") or higher_nm.lower().startswith(
        "candidatus"
    ):
        if lower_nm.lower().startswith("candidatus"):
            # Construct full name with first word of lower name removed
            nm = higher_nm.strip()
            for i in range(1, len(lower_nm.split())):
                nm = nm + " " + lower_nm.split()[i]
            return nm
        else:
            return higher_nm.strip() + " " + lower_nm.strip()
    # Catch too short species name case
    elif tri and len(higher_nm.split()) == 1:
        LOGGER.error(f"Species name ({higher_nm}) too short")
        return None
    # Then check that lower name is more words than the higher name
    elif len(lower_nm.split()) > len(higher_nm.split()):
        if lower_nm.lower().startswith(higher_nm.lower()):
            return lower_nm
        else:
            LOGGER.error(
                f"{L} name ({lower_nm}) appears to be {type} but "
                f"does not contain {H.lower()} name ({higher_nm})"
            )
            return None
    else:
        LOGGER.error(f"{H} name ({higher_nm}) appears to be too long")
        return None

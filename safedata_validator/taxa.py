"""This module describes classes used to compile taxonomic data from datasets.

This taxonomy can be validated against the GBIF backbone database. Alternatively,
taxonomic information from sequencing can be taken on a trust basis, in which case
checks are performed to catch badly formatted taxonomy data but the provided taxonomy is
accepted as is without being checked against a taxonomic authority.

The GBIFTaxon dataclass is used to store data about a taxon entry in a dataset. They are
initialised with user data and then the GBIFValidator class can be used to update a
Taxon object with the result of validation against a local version of the GBIF taxon
database.

Parallel 'Taxa' worksheets (GBIFTaxa and SeqTaxa) are defined, which are used
to load and collate the set of taxonomic entries from a dataset. These are then
collected in a higher level Taxa object, which additionally records the names used in
the Data worksheets. This allows us to check that all defined names are used, all used
names are defined, and that no names are defined in both Taxa worksheets (if both sheets
are provided).

Note that we explicitly exclude form and variety from the set of GBIF backbone taxonomic
levels because they cannot be matched into the backbone hierarchy without extra API
calls.
"""

import dataclasses
import re
import sqlite3
from collections import Counter
from io import StringIO
from itertools import compress, groupby
from typing import Optional

from dominate import tags
from dominate.util import raw
from openpyxl import worksheet

from safedata_validator.logger import (
    FORMATTER,
    LOGGER,
    get_handler,
    loggerinfo_push_pop,
)
from safedata_validator.resources import Resources
from safedata_validator.validators import (
    GetDataFrame,
    HasDuplicates,
    IsLower,
    blank_value,
)

ALL_BACKBONE_RANKS = [
    "domain",
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
]

# GBIF backbone doesn't include the two highest-level ranks ("domain", "superkingdom")
GBIF_BACKBONE_RANKS = ALL_BACKBONE_RANKS[2:]

# For sequenced taxonomy the top level rank can be any of the top three ranks ("domain",
# "superkingdom", "kingdom")
SEQ_TOP_RANKS = ALL_BACKBONE_RANKS[0:3]

# Sequenced taxonomy can then also include anything else as a backbone rank apart from
# subspecies
SEQ_ADDITIONAL_RANKS = ALL_BACKBONE_RANKS[3:-1]

# NBCI name regex
NCBI_prefix_re = re.compile("^[a-z]__")


class GBIFError(Exception):
    """Exception class for GBIF errors.

    Attributes:
        message: explanation of the error
    """

    def __init__(self, message="GBIF ID not found"):
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
    gbif_id: int | None = None
    is_backbone: bool = dataclasses.field(init=False)
    is_canon: bool = dataclasses.field(init=False)
    # https://stackoverflow.com/questions/33533148
    canon_usage: Optional["GBIFTaxon"] = dataclasses.field(init=False)
    parent_id: int | None = dataclasses.field(init=False)
    taxon_status: str | None = dataclasses.field(init=False)
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
        self.is_backbone = self.rank in GBIF_BACKBONE_RANKS
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
            taxon_rows = self.gbif_conn.execute(
                f"select * from backbone where canonical_name ='{taxon.name}' "
                f"and rank= '{taxon.rank.upper()}';"
            ).fetchall()
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
        taxon_row = self.gbif_conn.execute(
            f"select * from backbone where id = {gbif_id}"
        ).fetchone()

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
        taxon.parent_id = taxon_row["parent_key"]

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
            for rk, ky in [(r, r + "_key") for r in GBIF_BACKBONE_RANKS[:-1]]
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

        This method loads a set of taxa from the rows of a `safedata` formatted GBIFTaxa
        worksheet and populates the taxonomic hierarchy for those rows. The GBIFTaxa
        instance is updated.

        Args:
            worksheet: An openpyxl worksheet instance using the GBIFTaxa formatting
        """
        handler = get_handler()
        start_errors = handler.counters["ERROR"]

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

        # Fields used to describe taxa
        tx_fields = {
            "name",
            "taxon name",
            "taxon type",
            "taxon id",
            "ignore id",
            "parent name",
            "parent type",
            "parent id",
        }

        # Now check for extra fields and report them to the user
        extra_fields = set(headers).difference(tx_fields)
        if extra_fields:
            LOGGER.info("Additional fields provided: ", extra={"join": extra_fields})

        # Any duplication in names
        dupl_taxon_names = HasDuplicates([dframe.data_columns[headers.index("name")]])

        if dupl_taxon_names:
            LOGGER.error(
                "Duplicated names found: ", extra={"join": dupl_taxon_names.duplicated}
            )

        # get dictionaries of the taxa
        taxa = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]

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
            parent_info: list | None = [
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
        self.n_errors = handler.counters["ERROR"] - start_errors
        if self.n_errors is None:
            LOGGER.critical("GBIFTaxa error logging has broken!")
        elif self.n_errors > 0:
            LOGGER.info(f"GBIFTaxa contains {self.n_errors} errors")
        else:
            LOGGER.info(f"{len(self.taxon_names)} taxa loaded correctly")

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
            LOGGER.error(f"Worksheet name has whitespace padding: {m_name!r}")
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
                    LOGGER.error(f"{idx_name} has whitespace padding: {val!r}")
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
                LOGGER.error(f"{idx_name} has whitespace padding: {val!r}")
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
                self.taxon_index.append(
                    [
                        None,
                        p_taxon.gbif_id,
                        p_taxon.parent_id,
                        p_taxon.name,
                        p_taxon.rank,
                        p_taxon.taxon_status,
                    ]
                )

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
                    )
                    self.taxon_index.append(
                        [
                            None,
                            p_taxon.canon_usage.gbif_id,
                            p_taxon.canon_usage.parent_id,
                            p_taxon.canon_usage.name,
                            p_taxon.canon_usage.rank,
                            p_taxon.canon_usage.taxon_status,
                        ]
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
        # else:
        #         LOGGER.info('No parent taxon provided')

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
                # Update index - no taxon hierarchy except for parent
                self.taxon_index.append(
                    [m_name, -1, p_taxon.gbif_id, m_taxon.name, m_taxon.rank, "user"]
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
                # Update index - no taxon hierarchy except for parent
                self.taxon_index.append(
                    [m_name, -1, p_taxon.gbif_id, m_taxon.name, m_taxon.rank, "user"]
                )

        else:
            # Otherwise try and validate backbone taxon
            m_taxon = self.validator.search(m_taxon)

            if m_taxon.found and p_taxon is None:
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

                self.hierarchy.update(
                    [rw for rw in m_taxon.hierarchy if rw[1] is not None]
                )

                # Good backbone with no parent, provide info on taxon status
                if m_taxon.is_canon:
                    LOGGER.info(
                        f"Taxon found in GBIF backbone ({m_taxon.taxon_status})"
                    )
                elif m_taxon.canon_usage:
                    LOGGER.warning(
                        f"Taxon considered a {m_taxon.taxon_status} "
                        f"of {m_taxon.canon_usage.name} in GBIF backbone"
                    )

                    # Add the canon index entry and update hierarchy
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
                    self.hierarchy.update(
                        [
                            rw
                            for rw in m_taxon.canon_usage.hierarchy
                            if rw[1] is not None
                        ]
                    )

            elif m_taxon.found and p_taxon is not None:
                if p_taxon.found:
                    # Good backbone with good parent - are they compatible? Check if all
                    # entries in the parent hierarchy appear in the taxon hierarchy
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

                # Add to index and hierarchy
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
                self.hierarchy.update(
                    [rw for rw in m_taxon.hierarchy if rw[1] is not None]
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
        to_add.sort(key=lambda val: GBIF_BACKBONE_RANKS.index(val[0]))

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


class SeqTaxa:
    """Manage a set of taxon data derived from a sequencing workflow.

    This class to manage the generation of a taxon index from taxon tables generated
    through bioinformatics pipelines. It is a high-trust taxon table implementation that
    accepts a typically machine-generated taxon table and simply compiles a taxon
    hierarchy from the table.

    i)  the taxon_names attribute of the dataset, which is just a set of
        names used as a validation list for taxon names used in data worksheets.
    ii) the taxon_index attribute of the dataset, which contains a set
        of lists structured as:

            [worksheet_name (str),
            taxon_id (int),
            parent_id (int),
            canonical_name (str),
            taxonomic_rank (str),
            status (str)]

    Each taxon is assigned an arbitrary (negative) ID number. These are needed so that
    the taxon index follows the same format as for the GBIF validated case. These ID
    numbers are all negative so as to prevent any possible confusion with GBIF ID
    numbers, which refer to actual entries in the GBIF taxonomy database.

    The index can then be used:

    a) to generate the taxonomic coverage section of the dataset description, and
    b) to populate a database table to index the taxonomic coverage of datasets.

    Args:
        sheet_name: The name of the sheet that the specific SeqTaxa instance corresponds
            to
        database_name: The name of the database that the sequencing taxonomy has been
            resolved using
        database_version: The specific database version used to resolve the taxonomy
        database_link: Link (optional) to where the database can be found

    Attributes:
        taxon_index: A list containing taxon index lists
        taxon_names: A set of worksheet names
        hierarchy: A set of lists containing the complete taxonomic hierarchy for taxa
            in the SeqTaxa instance.
    """

    def __init__(
        self,
        sheet_name: str,
        database_name: str,
        database_version: str,
        database_link: str | None,
    ) -> None:
        self.taxon_index: list[tuple] = []
        self.taxon_names: set[str] = set()
        self.n_errors: int = 0
        self.sheet_name = sheet_name
        self.database_name = database_name
        self.database_version = database_version
        self.database_link = database_link

    @loggerinfo_push_pop("Loading sequenced taxonomy worksheet")
    def load(self, worksheet: worksheet) -> None:
        """Populate an SeqTaxa instance from an Excel worksheet.

        This method loads a set of taxa from the rows of a `safedata` formatted SeqTaxa
        worksheet and populates the taxonomic hierarchy for those rows.

        Args:
            worksheet: An openpyxl worksheet instance using the SeqTaxa formatting
        """
        handler = get_handler()
        start_errors = handler.counters["ERROR"]

        # Get the data read in, handling header issues like whitespace padding
        LOGGER.info(f"Reading bioinformatics taxon data from {self.sheet_name}")
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

        # Only the name field is indispensable
        if "name" not in headers:
            LOGGER.error("Sequencing taxa sheet is missing the name fields")
            FORMATTER.pop()
            return

        # Check that at least one top-level rank is provided, and that both domain and
        # superkingdom aren't provided
        top_ranks = set(SEQ_TOP_RANKS).intersection(headers)
        if len(top_ranks) == 0:
            LOGGER.error("At least one top-level taxonomic rank must be provided!")
            FORMATTER.pop()
            return
        elif "domain" in top_ranks and "superkingdom" in top_ranks:
            LOGGER.error(
                "Cannot provide both 'domain' and 'superkingdom' as taxonomic ranks!"
            )
            FORMATTER.pop()
            return

        if "domain" in top_ranks:
            highest_rank = "domain"
        elif "superkingdom" in top_ranks:
            highest_rank = "superkingdom"
        else:
            highest_rank = "kingdom"

        # It is acceptable to not provide any additional ranks beyond the top level one.
        # But if additional ranks are provided there can be no gaps between the lowest
        # provided rank and the top level ranks
        lower_ranks = set(SEQ_ADDITIONAL_RANKS).intersection(headers)
        if lower_ranks:
            lowest_rank = next(
                x for x in reversed(SEQ_ADDITIONAL_RANKS) if x in lower_ranks
            )
            lowest_rank_index = SEQ_ADDITIONAL_RANKS.index(lowest_rank) + 1
        else:
            lowest_rank_index = 0

        missing_ranks = set(SEQ_ADDITIONAL_RANKS[:lowest_rank_index]).difference(
            headers
        )
        if missing_ranks:
            LOGGER.error(
                "Need to provide all taxonomic ranks higher than current lowest "
                f"rank ({lowest_rank}) in SeqTaxa, missing ranks are as follows: ",
                extra={"join": missing_ranks},
            )
            FORMATTER.pop()
            return

        # List the ranks used in descending order. When only one top rank is provided,
        # then its just added. If two are provided the second one has to be kingdom and
        # first rank is filled by the other one.
        if len(top_ranks) == 1:
            ordered_ranks = [
                highest_rank,
                *SEQ_ADDITIONAL_RANKS[:lowest_rank_index],
            ]
        else:
            ordered_ranks = [
                highest_rank,
                "kingdom",
                *SEQ_ADDITIONAL_RANKS[:lowest_rank_index],
            ]

        # Now report extra fields (non-backbone ranks and other information)
        extra_fields = set(headers).difference(
            [*SEQ_TOP_RANKS, *SEQ_ADDITIONAL_RANKS, "name"]
        )
        if extra_fields:
            LOGGER.info("Additional fields provided: ", extra={"join": extra_fields})

        # Get dictionaries of the taxa
        taxa = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]
        FORMATTER.pop()

        # check number of taxa found
        if len(taxa) == 0:
            LOGGER.info("No taxon rows found")
            return

        # Store cleaned information as lists of taxon tuple - this is used to provide a
        # clean indexing system to build internally consistent parent child taxon ids
        # for the table.
        cleaned_taxa = {}

        # Clean and validate each taxon row
        for idx, row in enumerate(taxa):
            # Start validating the row
            LOGGER.info(f"Loading row {idx + 1}: {row['name']}")
            FORMATTER.push()

            # Get the worksheet row name
            worksheet_name = row["name"]

            if not isinstance(worksheet_name, str):
                LOGGER.error(f"Worksheet name is not a string: {worksheet_name!r}")
                worksheet_name = str(worksheet_name)
            else:
                worksheet_name_strip = worksheet_name.strip()
                if worksheet_name != worksheet_name_strip:
                    LOGGER.error(
                        f"Worksheet name has whitespace padding: {worksheet_name!r}"
                    )
                    worksheet_name = worksheet_name_strip

            self.taxon_names.add(worksheet_name)

            # Standardise blank and NA values to None
            row = {
                ky: None if blank_value(vl) or vl == "NA" else vl
                for ky, vl in row.items()
            }

            # Loop over rank fields to populate a cleaned taxon hierarchy:
            # - Tackle in taxonomic order by iterating over required ranks
            # - Drop empty entries
            # - Validate non-empty entries as unpadded strings
            # - Strip any NCBI k__ notation to match entries in names.names_txt db
            #   field. Runs from root, so cleans genus, species, subspecies.

            taxon_rank_tuple: list[tuple[str, str]] = []

            for rnk in ordered_ranks:
                # Get the name value associated with the rank
                value = row[rnk]

                # Genus name is needed to construct the species binomial name. It is set
                # as unknown until such time as it turns out that a valid name has been
                # provided (which then overwrites it)
                if rnk == "genus":
                    last_genus = "<genus unknown>"

                # Don't copy empty entries
                if value is None:
                    if rnk == highest_rank:
                        LOGGER.error(
                            f"Highest taxonomic rank ({rnk}) must be populated!"
                        )
                        break
                    else:
                        continue

                # The value must be an unpadded and not empty string
                if not isinstance(value, str) or value.isspace():
                    LOGGER.error(
                        f"Rank {rnk} has non-string or empty string value: {value!r}"
                    )
                    continue

                # The value must not be padded but processing can continue
                value_stripped = value.strip()
                if value != value_stripped:
                    LOGGER.error(f"Rank {rnk} has whitespace padding: {value!r}")
                    value = value_stripped

                # Strip k__ notation to provide clean name_txt search input - dropping
                # levels no taxonomic information is associated with the annotation (s__
                # etc. entries)
                value = taxa_strip(value, rnk)
                if value is None:
                    if rnk == highest_rank:
                        LOGGER.error(
                            f"Highest taxonomic rank ({rnk}) must be populated!"
                        )
                        break
                    else:
                        continue

                # Also remove any additional tags in front of the name, e.g. candidatus
                value = remove_additional_tags(value)

                # Genus name is now known to be properly formatted so we overwrite the
                # "<genus unknown>" placeholder with it
                if rnk == "genus":
                    last_genus = value

                if rnk == "species":
                    # Log an error if the species name appears to be a binomial
                    if len(value.split()) > 1:
                        LOGGER.error(
                            "Provided species name appears to be a binomial (which "
                            f"isn't allowed): {value}"
                        )
                        break

                    value = f"{last_genus} {value}"

                taxon_rank_tuple.append((rnk, value))

            # Add cleaned taxon tuples to list and report
            if taxon_rank_tuple:
                cleaned_taxa[worksheet_name] = taxon_rank_tuple
                leaf = taxon_rank_tuple[-1]
                LOGGER.info(f"Loaded {leaf[0]}: {leaf[1]}")
            else:
                LOGGER.info(f"Failed to load taxon {worksheet_name}")

            FORMATTER.pop()

        # Build the taxon index

        # Assign an negative arbitrary ID number to each unique pair of taxon rank and
        # name across the dataset. Negative numbers are used so that they cannot be
        # confused with the real ID numbers used in the GBIF case
        all_ranks = set([rank_pair for tx in cleaned_taxa.values() for rank_pair in tx])
        all_ranks_index = {
            rank_pair: val
            for val, rank_pair in zip(range(-1, -len(all_ranks) - 1, -1), all_ranks)
        }

        unique_taxa: set[tuple[None, int, int | None, str, str, str]] = set()

        for ws_name, taxon_details in cleaned_taxa.items():
            # Add taxa from the root to the tip, maintaining the chain of internal
            # ID values, and using the None index to represent the root node.
            lower_index = None
            for taxon_pair in taxon_details:
                this_index = all_ranks_index[taxon_pair]
                unique_taxa.add(
                    (
                        None,
                        this_index,
                        lower_index,
                        taxon_pair[1],
                        taxon_pair[0],
                        "loaded",
                    )
                )
                lower_index = this_index

            self.taxon_index = list(unique_taxa)

        FORMATTER.pop()

        # summary of processing
        self.n_errors = handler.counters["ERROR"] - start_errors
        if self.n_errors is None:
            LOGGER.critical("SeqTaxa error logging has broken!")
        elif self.n_errors > 0:
            LOGGER.info(f"{self.sheet_name} contains {self.n_errors} errors")
        else:
            LOGGER.info(
                f"{len(self.taxon_names)} taxa loaded correctly from {self.sheet_name}"
            )

        FORMATTER.pop()

    @property
    def is_empty(self) -> bool:
        """Check if an SeqTaxa instance contains any taxa."""
        return len(self.taxon_names) == 0


class Taxa:
    """Manage combined taxon sheet instances.

    This class wraps taxon sheets and provides shared properties across the instances

    Args:
        resources: A Resources instance


    We are interested in checking that:
    * no worksheet names are reused when more than one taxon sheet are provided,
    * every worksheet name is used somewhere in the Data worksheets, and
    * every taxon name used across the Data worksheets is defined in a Taxa worksheet.

    This overarching class stores instances of the lower level taxon sheet handler
    classes. It can also store (as `taxon_names_used`) the set of all names used across
    Data worksheets. The property `is_empty` can be used to check whether lower level
    classes are empty, and the property `taxon_names` can be used to find the set of all
    taxon names defined across taxon handlers. Finally, the property `repeat_names` can
    be used to find if any names are used in both GBIFTaxa and SeqTaxa worksheets.
    """

    def __init__(self, resources: Resources):
        self.gbif_taxa = GBIFTaxa(resources)
        self.seq_taxa_sheets: list[SeqTaxa] = []
        self.taxon_names_used: set[str] = set()

    @property
    def is_empty(self) -> bool:
        """Reports if neither GBIF nor sequenced taxa are loaded."""
        return self.gbif_taxa.is_empty and len(self.seq_taxa_sheets) == 0

    @property
    def taxon_names(self) -> set[str]:
        """Provides loaded taxon names from all taxon handlers."""

        return set(
            [
                *self.gbif_taxa.taxon_names,
                *[name for sheet in self.seq_taxa_sheets for name in sheet.taxon_names],
            ]
        )

    @property
    def repeat_names(self) -> set[str]:
        """Reports taxon names duplicated between taxon handlers."""

        seen = set()
        duplicated = set()

        all_names = [
            *self.gbif_taxa.taxon_names,
            *[name for sheet in self.seq_taxa_sheets for name in sheet.taxon_names],
        ]

        for this_name in all_names:
            if this_name not in seen:
                seen.add(this_name)
            else:
                duplicated.add(this_name)

        return duplicated


def taxon_index_to_text(
    taxa: list[dict],
    html: bool = False,
    indent_width: int = 4,
    lowest_taxa: str | None = None,
) -> str | tags.div:
    """Render a taxon index as text or html.

    This function takes a taxon index and renders the contents into either a text or
    html representation of the taxonomic hierarchy used in the dataset. Taxonomic ranks
    are indented to render a nested hierarchy.

    Args:
        taxa: A list of taxon dictionaries containing the taxa for a dataset.
        html: Render as html or text.
        indent_width: The indentation width to use for successive taxonomic ranks.
        lowest_taxa: The lowest taxonomic rank that the index renders, if no rank is
            provided then the index is rendered for all ranks.

    Returns:
        Either a HTML or text representation of the taxa tree.
    """

    def _indent(n: int, use_html: bool = html):
        if use_html:
            return raw("&ensp;-&ensp;" * n)
        else:
            return " " * indent_width * (n - 1)

    def _format_name(tx: dict, use_html: bool = html):
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

    # Eliminate any taxa with ranks below the minimum
    if lowest_taxa:
        # Check that the lowest rank appears in the full set of taxa
        if lowest_taxa not in ALL_BACKBONE_RANKS:
            raise ValueError(
                f"Rank provided to render taxa tree down to {lowest_taxa} is not a "
                f"backbone rank! Should be one of: {ALL_BACKBONE_RANKS}"
            )

        # Generate the full list of ranks that should be rendered
        rendered_ranks = ALL_BACKBONE_RANKS[: ALL_BACKBONE_RANKS.index(lowest_taxa) + 1]

        # Then add any taxa that have ranks that aren't in the list of rendered ranks to
        # the superfluous taxa index
        for idx, taxon in enumerate(taxa):
            if taxon["taxon_rank"] not in rendered_ranks:
                surp_tx_ids.append(idx)

    # Delete taxa that are superfluous by index
    for index in sorted(set(surp_tx_ids), reverse=True):
        del taxa[index]

    # group taxa by their parent id
    grouped = {k: list(v) for k, v in groupby(taxa, lambda x: x["parent_id"])}

    # start the stack with root taxa, which will have None as a parent (kingdoms for
    # GBIF, kingdoms/superkingdoms/domains for sequenced taxa)
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


def taxa_strip(name: str, rank: str) -> str | None:
    """Strip NCBI style rank prefixes from taxon names.

    This function removes NCBI `k__` type notation from taxa names. It returns the
    stripped name and emits an error if the removed prefix is inconsistent with the
    provide rank. If a name consists _only_ of the `s__` style annotation, indicating a
    match to an unknown sequence at species level, then None is returned.

    Args:
        name: An NCBI taxon name with `k__` style rank prefix
        rank: The expected taxonomic rank for the name.

    Returns:
        A trimmed string or None in the case of anonymous `k__` style annotation.
    """

    prefix_match = NCBI_prefix_re.match(name)

    if prefix_match is not None:
        # Strip the name down and check rank consistency
        s_name = name[prefix_match.end() :]
        if name[0].lower() != rank[0].lower():
            LOGGER.error(f"Prefix of taxon {name} inconsistent with rank {rank}")

        return None if s_name == "" else s_name

    return name


def remove_additional_tags(taxon_name: str) -> str:
    """Remove additional tags from taxon names.

    Reference taxonomy databases sometimes attach additional tags in front of taxon
    names (e.g. 'candidatus' to indicate that a taxon has never been cultured). We do
    not want to clog up our taxon names with these tags, so this function removes them.
    It does this by checking if the provided taxon name starts with a known tag of this
    type and removing it if so.

    Args:
        taxon_name: The name provided for the taxon.

    Returns:
        The taxon name with additional tags stripped out.
    """

    known_tags = ["candidatus"]

    if len(taxon_name.split()) >= 1 and any(
        taxon_name.lower().startswith(tag) for tag in known_tags
    ):
        return " ".join(taxon_name.split()[1:])
    else:
        return taxon_name

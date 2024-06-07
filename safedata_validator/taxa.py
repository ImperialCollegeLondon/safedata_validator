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
import re
import sqlite3
from collections import Counter
from io import StringIO
from itertools import compress, groupby, pairwise
from typing import Any, Optional

from dominate import tags
from dominate.util import raw
from openpyxl import worksheet

from safedata_validator.logger import (
    FORMATTER,
    LOGGER,
    get_handler,
    log_and_raise,
    loggerinfo_push_pop,
)
from safedata_validator.resources import Resources
from safedata_validator.validators import (
    GetDataFrame,
    HasDuplicates,
    IsLower,
    blank_value,
)

GBIF_BACKBONE_RANKS = [
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
NCBI_BACKBONE_RANKS = ["superkingdom", *GBIF_BACKBONE_RANKS]

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


class NCBIError(Exception):
    """Exception class for NCBI errors.

    Attributes:
        message: explanation of the error
    """

    def __init__(
        self, message="No entry found for NCBI ID: merged, deleted or invalid"
    ):
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


@dataclasses.dataclass
class NCBITaxon:
    """Represent and validate an NCBI taxon.

    Initialised using user taxonomic information for single taxon, which can be then be
    validated against the NCBI database. Attributes are populated when an instance is
    passed to NCBIValidator.

    The format of the hierarchy is expected to be a list in the format generated by the
    [`_get_canon_hierarchy`][safedata_validator.taxa.NCBIValidator._get_canon_hierarchy]
    method.

    Args:
        name: A taxonomic name
        rank: A taxonomic rank
        taxa_hier: A taxonomic hierarchy list
        ncbi_id: NCBI ID for full taxa (i.e. including non-backbone ranks)
        parent_ncbi_id: The NCBI ID of the parent backbone taxon
        non_canon_name: An optional string giving a non canon name
        non_canon_name_class: The NCBI name class of the non canon name
    """

    # Init properties
    name: str
    rank: str
    taxa_hier: list[tuple[Any, ...]]
    ncbi_id: int
    parent_ncbi_id: int
    non_canon_name: str | None = None
    non_canon_name_class: str | None = None

    def __post_init__(self):
        """Sets the defaults for the post-init properties and checks inputs."""

        if not isinstance(self.name, str):
            raise TypeError("Provided taxon name not a string")

        if not isinstance(self.rank, str):
            raise TypeError("Provided rank not in string form")

        # Screen ncbi id values: int-like
        if not (
            isinstance(self.ncbi_id, int)
            or (isinstance(self.ncbi_id, float) and self.ncbi_id.is_integer())
        ):
            raise TypeError("NCBI ID is not an integer")

        self.ncbi_id = int(self.ncbi_id)

        # Screen parent id values: None or int-like
        if self.parent_ncbi_id is None:
            pass
        elif (
            isinstance(self.parent_ncbi_id, float) and self.parent_ncbi_id.is_integer()
        ):
            self.parent_ncbi_id = int(self.parent_ncbi_id)
        elif not isinstance(self.parent_ncbi_id, int):
            raise TypeError("Parent NCBI ID is not an integer or None")

        if not isinstance(self.taxa_hier, list):
            raise TypeError("Taxon hierarchy not a list")

        if len(self.taxa_hier) == 0:
            raise ValueError("Taxon hierarchy empty")

        if not all(isinstance(x, tuple) for x in self.taxa_hier):
            raise ValueError("Taxon hierarchy values not all tuples")

        # Check the tuple values: rank, name, ncbi_id and parent_ncbi_id/None,
        # where None shows the root of the taxonomy.
        tuple_contents = {tuple(map(type, x)) for x in self.taxa_hier}
        valid_tuples = {(str, str, int, int), (str, str, int, type(None))}
        if not tuple_contents.issubset(valid_tuples):
            raise ValueError("Taxon hierarchy tuples malformed")

        # rank to lowercase
        self.rank = self.rank.lower()

        # Extract elements that should match the hierarchy leaf
        leaf_rank, leaf_name, leaf_ncbi_id, parent_ncbi_id = self.taxa_hier[0]

        if self.rank != leaf_rank:
            raise ValueError(
                f"Provided rank ({self.rank}) does not match "
                f"first rank in hierarchy ({leaf_rank})"
            )

        if self.name != leaf_name:
            raise ValueError(
                f"Provided taxon name ({self.name}) does not match "
                f"first name in hierarchy ({leaf_name})"
            )

        if self.ncbi_id != leaf_ncbi_id:
            raise ValueError(
                f"Provided NCBI ID ({self.ncbi_id}) not does not match "
                f"first ID in hierarchy ({leaf_ncbi_id})"
            )

        if self.parent_ncbi_id != parent_ncbi_id:
            raise ValueError(
                f"Provided parent NCBI ID ({self.parent_ncbi_id}) not does not match "
                f"first parent ID in hierarchy ({parent_ncbi_id})"
            )

    def __repr__(self):
        """Provides a simple representation of the class."""

        if self.non_canon_name is not None:
            return (
                f"{self.name} (resolved as {self.rank} "
                f"rather than {self.non_canon_name})"
            )
        else:
            return f"{self.name}"


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
        self.ncbi_conn: sqlite3.Connection = conn
        """A connection to the NCBI sqlite3 database."""
        # Retrieve the taxon ranking and extra taxa from the DB, explicitly
        # sorting from root to leaf
        self.ncbi_ranks_root_to_leaf: list[str] = [
            rw["rank"]
            for rw in self.ncbi_conn.execute(
                "select rank from unique_ncbi_ranks order by rank_index;"
            )
        ]
        """A root to leaf list of NCBI taxonomic ranks, excluding clade and no rank."""

        self.superkingdoms: list[str] = [
            rw["name_txt"]
            for rw in self.ncbi_conn.execute(
                "select name_txt "
                "from names join nodes using (tax_id) "
                "where rank='superkingdom'"
            ).fetchall()
        ]
        """A list of all names assigned superkingdom rank in the NCBI database."""

        self.kingdoms: list[str] = [
            rw["name_txt"]
            for rw in self.ncbi_conn.execute(
                "select name_txt "
                "from names join nodes using (tax_id) "
                "where rank='kingdom'"
            ).fetchall()
        ]
        """A list of all names assigned kingdom rank in the NCBI database."""

    def __del__(self) -> None:
        """Delete a LocalNCBIValidator instance.

        This method ensures that the database connection is closed correctly.
        """

        self.ncbi_conn.close()

    def id_lookup(self, ncbi_id: int) -> NCBITaxon:
        """Get an NCBITaxon by NCBI ID.

        This method returns a populated NCBITaxon instance given an NCBI ID.

        Args:
            ncbi_id: Unique identifier for the taxon

        Returns:
            A populated NCBITaxon instance

        Raises:
            NCBIError: if the provided ID cannot be found.
        """

        if not isinstance(ncbi_id, int):
            raise TypeError("Non-integer NCBI taxonomy ID")

        if not ncbi_id > 0:
            raise ValueError("Negative NCBI taxonomy ID")

        # Get complete hierarchy dictionaries on the taxon ID, which traps bad IDs
        hierarchy = self._get_canon_hierarchy(ncbi_id)

        # Reduce the hierarchy to the ranks only also present in the backbone hierarchy
        bb_hierarchy = self._canon_to_backbone_hierarchy(hierarchy)

        # Create and populate taxon using details from the leaf taxon at position 0.
        # Note that this handles the edge case where there is only a single level in the
        # bb_hierarchy and that has just had the parent ID set to None to mark the root.
        taxon_rank, taxon_name, ncbi_id, parent_ncbi_id = bb_hierarchy[0]

        mtaxon = NCBITaxon(
            name=taxon_name,
            rank=taxon_rank,
            ncbi_id=ncbi_id,
            parent_ncbi_id=parent_ncbi_id,
            taxa_hier=bb_hierarchy,
        )

        return mtaxon

    def _get_canon_hierarchy(self, ncbi_id: int) -> list[tuple[Any, ...]]:
        """Extract the complete canonical NCBI hierarchy for a taxon ID.

        This function compiles complete taxon hierarchy as a list ordered from tip to
        root. Each list entry is a tuple providing the taxon rank, canonical scientific
        name, taxon id and parent taxon id, starting with the provided taxon id at the
        tip.

        The canonical usage is taken from the row associated with each NCBI taxon ID
        that is recorded as having the 'scientific name' name class. The SQL query below
        demonstrates that all IDs have a canonical usage: there are no rows found with
        taxon IDs that do not also appear with 'scientific name' status.

        ``` SQL
            SELECT count(*) FROM names
                WHERE name_class != "scientific name"
                AND tax_id NOT IN (
                    SELECT tax_id FROM names
                    WHERE name_class == "scientific name"
                );
        ```

        Args:
            ncbi_id: A valid NCBI taxon ID.

        Returns:
            An ordered dictionary of parent details.

        Raises:
            ValueError: This function is intended only to be used with valid current
                NCBI IDs and raises an error for invalid, merged or deleted NCBI IDs.
        """

        # SQL to retrieve canon name usage for a taxon id
        canon_sql = """
            select rank, name_txt, tax_id, parent_tax_id
                from names join nodes using (tax_id)
                where tax_id = {}
                and name_class = 'scientific name';"""

        # Get information on the taxon ID
        name_row = self.ncbi_conn.execute(canon_sql.format(ncbi_id)).fetchone()

        if name_row is None:
            raise NCBIError()

        # Now build the taxonomic hierarchy - this will yield a list of tuples of
        # taxa in taxonomic order from leaf to root.
        # TODO - ideally this would be typed as tuple[str, str, int, int] but
        #        tuple(sqlite3.Row) is typed as tuple[Any, ...]
        hierarchy = [tuple(name_row)]

        while name_row["parent_tax_id"] != 1:
            # Find node and name of the parent taxon
            name_row = self.ncbi_conn.execute(
                canon_sql.format(name_row["parent_tax_id"])
            ).fetchone()

            hierarchy.append(tuple(name_row))

        return hierarchy

    @staticmethod
    def _canon_to_backbone_hierarchy(
        canon_hier: list[tuple[Any, ...]],
    ) -> list[tuple[Any, ...]]:
        """Reduce a canonical NCBI hierarchy to the backbone ranks.

        This method reduces a canonical hierarchy to only the ranks present in the list
        of backbone taxa, which are used in the taxon index for the elements. This
        cleanly omits missing backbone elements, if e.g. genus has an order parent.
        """

        bb_elements = [hr for hr in canon_hier if hr[0] in NCBI_BACKBONE_RANKS]

        # If the first element of the canon hierarchies is not of a backbone rank,
        # re-insert that element to record the leaf details
        if canon_hier[0][0] not in NCBI_BACKBONE_RANKS:
            bb_elements.insert(0, canon_hier[0])

        # Update the parent taxa ids to follow the backbone ranks, inserting None as the
        # root parent ID. If there is only one bb_element (backbone or non-backbone
        # leaf) that should be connected directly to the root.
        if len(bb_elements) > 1:
            bb_hierarchy = []

            for child, parent in pairwise(bb_elements):
                bb_hierarchy.append(tuple(list(child[:3]) + list([parent[2]])))

            bb_hierarchy.append(tuple(list(parent[:3]) + list([None])))
        else:
            bb_hierarchy = [tuple(list(bb_elements[0][:3]) + list([None]))]

        return bb_hierarchy

    def _check_congruent_hierarchies(
        self,
        expected_hier: list[tuple[Any, ...]],
        provided_hier: list[tuple[str, str]],
        report: bool,
    ) -> bool:
        """Check if a provided hierarchy is congruent with an expected hierarchy.

        The NCBITaxa sheet can provide a taxonomic hierarchy below a leaf taxon. This
        method can be used to compare the elements of that provided hierarchy with the
        complete expected hierarchy for a taxon.

        The expected hierarchy should be provided as the canonical hierarchy (see
        [_get_canon_hierarchy][safedata_validator.taxa.NCBIValidator._get_canon_hierarchy]).

        The method logs warnings when non-canon names are used and errors if any
        elements are completely incongruent between the two hierarchies.

        Args:
            expected_hier: The expected hierarchy
            provided_hier: A provided hierarchy that should contain matching elements.
            report: Should the function log errors when mismatches are found.

        Returns:
            A boolean indicating if the provided details are congruent with the
            expected.
        """

        congruent = True
        expected_hier_lookup = {v[0]: v for v in expected_hier}

        # Is each provided taxon_hier entry compatible with the expected hierarchy?
        for rank, provided_name in provided_hier:
            # Is the rank present in the expected hierarchy
            if rank not in expected_hier_lookup:
                congruent = False
                if report:
                    LOGGER.error(
                        f"Taxonomy mismatch for {provided_name} at rank {rank}: "
                        f"rank not found in expected hierarchy"
                    )
                continue

            # Get the canon name and id and check if the provided name is canon
            _, canon_name, tax_id, _ = expected_hier_lookup[rank]
            if canon_name == provided_name:
                continue

            # If not, is there a matching alternative name at the same rank and id?
            row = self.ncbi_conn.execute(
                f"select name_class\n"
                f"from names join nodes using (tax_id)\n"
                f"where name_txt = '{provided_name}' and rank = '{rank}' and tax_id = "
                f"'{tax_id}';"
            ).fetchone()

            if row is None:
                congruent = False
                if report:
                    LOGGER.error(
                        f"Taxonomy mismatch for {provided_name} at rank {rank}: "
                        f"expecting {canon_name}"
                    )
            elif report:
                LOGGER.warning(
                    f"Non-canon name {provided_name} at rank {rank}: "
                    f"{row['name_class']} for {canon_name}"
                )

        return congruent

    def _get_name_and_rank_matches(self, name: str, rank: str) -> list[sqlite3.Row]:
        """Get NCBI taxon details matching a given name and rank."""

        row_sql = """
        select tax_id, name_txt, parent_tax_id, rank, name_class
            from names join nodes using (tax_id)
            where name_txt = '{}' and rank = '{}';
        """

        return self.ncbi_conn.execute(row_sql.format(name, rank)).fetchall()

    def taxa_search(self, nnme: str, taxon_hier: list[tuple[str, str]]) -> NCBITaxon:
        """Find an NCBI taxon given a taxon hierarchy.

        This method takes a provided NCBI taxonomic hierarchy, and attempts to find a
        congruent NCBI ID. This NCBI ID is then used to generate a NCBITaxon object,
        which is returned. Where possible, this method makes use of parent taxa
        information to distinguish between ambiguous taxa names.

        Args:
            nnme: A nickname to identify the taxon
            taxon_hier: A list of the (rank, name) tuples provided in the taxon
                worksheet, sorted in root to tip taxonomic order

        Returns:
            A populated NCBITaxon instance

        Raises:
            ValueError: incorrect argument values
            NCBIError: no match found
        """

        # TODO - validation of inputs at multiple points - maybe make less paranoid? All
        #        of this has been sanitised upstream in validate_and_add_taxa.
        if not isinstance(nnme, str):
            raise ValueError("Taxon nickname should be a string")

        if isinstance(taxon_hier, list) is False:
            raise ValueError("Taxon hierarchy should be a list")

        if len(taxon_hier) == 0:
            raise ValueError("Empty taxon hierarchy")

        for entry in taxon_hier:
            if not (
                isinstance(entry, tuple)
                and len(entry) == 2
                and all(isinstance(x, str) for x in entry)
            ):
                raise ValueError(
                    "Not all taxa hierachy entries are 2 tuples of strings"
                )

        # Now ensure the taxon hierarchy is in root to leaf order and then get the leaf
        # and look for that name and rank combination in the database.
        taxon_hier.sort(key=lambda x: self.ncbi_ranks_root_to_leaf.index(x[0]))
        leaf_rank, leaf_name = taxon_hier.pop()
        taxon_rows = self._get_name_and_rank_matches(leaf_name, leaf_rank)

        # The outcome of that query can be:
        # - No records found with the leaf name and rank - fail
        # - One unique row found with the name and rank - continue
        # - Multiple records with the name and rank, in which case is one and only one
        #   compatible with the remaining taxon hierarchy data, if there is any.

        n_records = len(taxon_rows)
        mtaxon = None

        if n_records == 0:
            # No records found
            LOGGER.error(
                f"Taxa {nnme} not found with name {leaf_name} and rank {leaf_rank}"
            )
            raise NCBIError("Taxon not found")

        if n_records == 1:
            # Single record found - get the canon hierarchy and check congruence
            match_taxon = dict(taxon_rows[0])
            canon_hier = self._get_canon_hierarchy(match_taxon["tax_id"])

            # Warn if the usage is non-canon
            not_canon = match_taxon["name_class"] != "scientific name"
            if not_canon:
                LOGGER.warning(
                    f"Non-canon usage: {leaf_name} is {match_taxon['name_class']} "
                    f"for {canon_hier[0][1]}"
                )

            # Now check that any taxonomy above the taxon is congruent
            congruent = self._check_congruent_hierarchies(
                provided_hier=taxon_hier,
                expected_hier=canon_hier,
                report=True,
            )

            if congruent:
                LOGGER.info(f"Match found for {nnme}")
            else:
                LOGGER.error(f"Match found for {nnme} with incongruent taxonomy")
                raise NCBIError("Taxon not found")

            # TODO - id_lookup simply repeats the _get_canon_hierarchy step and then
            #        packages an NCBI taxon. We could redefine id_lookup as
            #        build_ncbi_taxon(..., canon_hier, ...)

            mtaxon = self.id_lookup(taxon_rows[0]["tax_id"])

        elif not taxon_hier:
            # No further taxonomic data to resolve multiple hits
            LOGGER.error(
                f"Multiple matches for taxon {nnme}: no additional taxonomy provided."
            )
            raise NCBIError("Taxon not found")

        else:
            # Get the taxonomic hierarchy for each candidate
            candidate_hier = [
                self._get_canon_hierarchy(rw["tax_id"]) for rw in taxon_rows
            ]

            # Check which hierarachies are compatible with the taxon_hier data
            compatible_rows = []
            for cand_row, cand_hier in zip(taxon_rows, candidate_hier):
                # If the provided hierarchy is congruent with this expected hierachy
                # then add it to the list of compatible NCBI ID
                congruent = self._check_congruent_hierarchies(
                    expected_hier=cand_hier,
                    provided_hier=taxon_hier,
                    report=False,
                )

                if congruent:
                    compatible_rows.append((cand_row, cand_hier))

            # How many compatible taxa are found: none, one or more?

            n_compatible = len(compatible_rows)
            if n_compatible == 0:
                LOGGER.error(
                    f"Multiple matches for taxon {nnme}: "
                    "provided taxonomy not congruent with candidates"
                )
                raise NCBIError("Taxon not found")

            elif n_compatible > 1:
                LOGGER.error(
                    f"Multiple matches for taxon {nnme}: "
                    "provided taxonomy does not resolve candidates"
                )
                raise NCBIError("Taxon not found")

            else:
                # Report on the match
                match_taxon = dict(compatible_rows[0][0])
                canon_hier = compatible_rows[0][1]

                # Warn if the usage is non-canon
                not_canon = match_taxon["name_class"] != "scientific name"
                if not_canon:
                    LOGGER.warning(
                        f"Non-canon usage: {leaf_name} is {match_taxon['name_class']} "
                        f"for {canon_hier[0][1]}"
                    )

                LOGGER.info(f"Match found for {nnme}")
                mtaxon = self.id_lookup(match_taxon["tax_id"])

        if not_canon:
            mtaxon.non_canon_name = leaf_name
            mtaxon.non_canon_name_class = match_taxon["name_class"]

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
        self.taxon_index: list[tuple] = []
        self.taxon_names: set[str] = set()
        self.hierarchy: set[tuple] = set()
        self.n_errors: int = 0

        # Get the validator instance
        self.validator = NCBIValidator(resources)

    @loggerinfo_push_pop("Loading NCBITaxa worksheet")
    def load(self, worksheet: worksheet) -> None:
        """Populate an NCBITaxa instance from an Excel worksheet.

        This method loads a set of taxa from the rows of a `safedata` formatted NCBITaxa
        worksheet and populates the taxonomic hierarchy for those rows. The GBIFTaxa
        instance is updated.

        Args:
            worksheet: An openpyxl worksheet instance using the GBIFTaxa formatting
        """
        handler = get_handler()
        start_errors = handler.counters["ERROR"]

        # Get the data read in, handling header issues like whitespace padding
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

        # Only the name field is indispensible
        if "name" not in headers:
            LOGGER.error("NCBI taxa sheet is missing the name fields")
            FORMATTER.pop()
            return

        # The set of possible data fields is name, new and the NCBI ranks
        data_fields = ["name", "new", *self.validator.ncbi_ranks_root_to_leaf]

        # Now report extra fields
        extra_fields = set(headers).difference(data_fields)
        if extra_fields:
            LOGGER.info("Additional fields provided: ", extra={"join": extra_fields})

        # Extract a root to leaf sorted list of taxonomic rank fields and report.
        found_ranks = [
            rnk for rnk in self.validator.ncbi_ranks_root_to_leaf if rnk in headers
        ]

        if not found_ranks:
            LOGGER.error("NCBI taxa sheet contains no taxonomic rank fields.")
            FORMATTER.pop()
            return

        LOGGER.info(
            f"{len(found_ranks)} NCBI rank fields found: ", extra={"join": found_ranks}
        )

        # The NCBI database uses full binomial and trinomial names for species and
        # subspecies, but NCBI data commonly only provides parts of those names under
        # the field headings, so require the fields needed to construct bi/trinomial
        if "subspecies" in found_ranks and "species" not in found_ranks:
            LOGGER.error("A species field is required with subspecies.")
            FORMATTER.pop()
            return

        if "species" in found_ranks and "genus" not in found_ranks:
            LOGGER.error("A genus field is required with species.")
            FORMATTER.pop()
            return

        # Any duplication in worksheet taxon names
        dupl_taxon_names = HasDuplicates(dframe.data_columns[headers.index("name")])

        if dupl_taxon_names:
            LOGGER.error(
                "Duplicated names found: ", extra={"join": dupl_taxon_names.duplicated}
            )

        # Get dictionaries of the taxa
        taxa = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]
        FORMATTER.pop()

        # check number of taxa found
        if len(taxa) == 0:
            LOGGER.info("No taxon rows found")
            return

        # Clean and validate each taxon row
        for idx, row in enumerate(taxa):
            # Start validating the row
            LOGGER.info(f"Validating row {idx + 1}: {row['name']}")
            FORMATTER.push()

            # Validate row name
            m_name = row["name"]

            if not isinstance(m_name, str):
                LOGGER.error(f"Worksheet name is not a string: {m_name!r}")
                m_name = str(m_name)
            else:
                m_name_strip = m_name.strip()
                if m_name != m_name_strip:
                    LOGGER.error(f"Worksheet name has whitespace padding: {m_name!r}")
                    m_name = m_name_strip

            # Standardise blank and NA values to None
            row = {
                ky: None if blank_value(vl) or vl == "NA" else vl
                for ky, vl in row.items()
            }

            # Loop over rank fields to populate a cleaned taxon hierarchy:
            # - Tackle in taxonomic order by iterating over ordered found_ranks
            # - Drop empty entries
            # - Validate non-empty entries as unpadded strings
            # - Strip any NCBI k__ notation to match entries in names.names_txt db
            #   field. Runs from root, so cleans genus, species, subspecies.

            taxon_dict: dict[str, str] = {}
            validate = True

            for rnk in found_ranks:
                # Get the name value associated with the rank
                value = row[rnk]

                # Don't copy empty entries
                if value is None:
                    continue

                # The value must be an unpadded and not empty string
                if not isinstance(value, str) or value.isspace():
                    LOGGER.error(
                        f"Rank {rnk} has non-string or empty "
                        f"string value: {value!r}"
                    )
                    validate = False
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
                    continue

                # Finally promote kingdoms to superkingdoms if required
                if (
                    rnk == "kingdom"
                    and value not in self.validator.kingdoms
                    and value in self.validator.superkingdoms
                ):
                    taxon_dict["superkingdom"] = value
                    LOGGER.warning(
                        f"NCBI records {taxon_dict['superkingdom']} as "
                        f"a superkingdom rather than a kingdom"
                    )
                else:
                    taxon_dict[rnk] = value

            # Now convert species and subspecies ranks to binomial, trinomial where
            # possible. Abandon validation for taxa where construct_bi_or_tri fails.

            if "species" in taxon_dict:
                try:
                    binomial = construct_bi_or_tri(
                        higher_nm=taxon_dict["genus"],
                        lower_nm=taxon_dict["species"],
                        tri=False,
                    )
                    taxon_dict["species"] = binomial
                except ValueError:
                    validate = False

            if "subspecies" in taxon_dict:
                try:
                    trinomial = construct_bi_or_tri(
                        higher_nm=row["species"],
                        lower_nm=taxon_dict["subspecies"],
                        tri=True,
                    )
                    taxon_dict["subspecies"] = trinomial
                except ValueError:
                    validate = False

            # Check new settings - empty cell (new=false) or a yes/no string.
            if "new" in row:
                isnew = row["new"]
                if isnew is None:
                    new = False
                elif not isinstance(isnew, str) or isnew.lower() not in ["yes", "no"]:
                    LOGGER.error("Values in the 'new' field must be yes or no.")
                    validate = False
                else:
                    new = False if isnew.lower() == "no" else True
            else:
                new = False

            # Validate if possible
            if validate:
                self.validate_and_add_taxon(
                    m_name=m_name,
                    taxon_hier=[(k, v) for k, v in taxon_dict.items()],
                    new=new,
                )
            FORMATTER.pop()

        FORMATTER.pop()

        # Add the higher taxa
        if self.hierarchy:
            self.index_higher_taxa()

        # summary of processing
        self.n_errors = handler.counters["ERROR"] - start_errors
        if self.n_errors is None:
            LOGGER.critical("NCBITaxa error logging has broken!")
        elif self.n_errors > 0:
            LOGGER.info(f"NCBITaxa contains {self.n_errors} errors")
        else:
            LOGGER.info(f"{len(self.taxon_names)} taxa loaded correctly")

        FORMATTER.pop()

    def validate_and_add_taxon(
        self, m_name: str, taxon_hier: list[tuple[str, str]], new: bool = False
    ) -> None:
        """Add a NCBI formatted taxon row to the NCBITaxa instance.

        This method takes a worksheet name and taxon hierarchy, validates it against the
        NCBI database and updates the NCBITaxa instance to include the new details. This
        is typically used to process rows found in a dataset with an NCBITaxa formatted
        table, but can also be used to populate a NCBITaxa instance programmatically.

        Args:
            m_name: The worksheet name used to refer to the taxon.
            taxon_hier: A list of tuples of the provided pairs of NCBI taxon ranks and
                taxon names, sorted in order from root to leaf/
            new: Is the leaf taxon new and not expected to be found in the NCBI
                database.
        """

        # Validation steps - separate input value issues, which should log and continue,
        # from straightforward programming errors in the structure of the data being
        # passed.

        # Sanitise worksheet names for taxa - only keep unpadded strings.
        if (
            m_name is None
            or not isinstance(m_name, str)
            or m_name.isspace()
            or not m_name
        ):
            LOGGER.error("Worksheet name missing, whitespace only or not text")
            return

        # Clean whitespace padding and warn
        m_name_strip = m_name.strip()
        if m_name != m_name_strip:
            LOGGER.error(f"Worksheet name has whitespace padding: {m_name!r}")
            m_name = m_name_strip

        # TODO  - validation at multiple levels. Make less paranoid? Here, an empty list
        #         is a possiblity from the data but the upstream code should have
        #         sanitised the rest so a value error is raised.
        if isinstance(taxon_hier, list) is False:
            log_and_raise("Taxon hierarchy should be a list", ValueError)

        if len(taxon_hier) == 0:
            LOGGER.error("No taxonomy provided")
            return

        # Detect malformed taxon hierarchy lists
        for entry in taxon_hier:
            if not (
                isinstance(entry, tuple)
                and len(entry) == 2
                and all(isinstance(x, str) for x in entry)
                and all(not x.isspace() for x in entry)
            ):
                log_and_raise(
                    "Not all taxa hierachy entries are 2 tuples of non-empty strings",
                    ValueError,
                )

        # Detect padded taxon hierarchy entries
        for entry_idx, entry in enumerate(taxon_hier):
            stripped = (entry[0].strip(), entry[1].strip())
            if stripped != entry:
                LOGGER.error(
                    "Hierarchy contains whitespace: "
                    f"rank {entry[0]!r}, name {entry[1]!r}"
                )
                taxon_hier[entry_idx] = stripped

        # Handle new entries - remove the tip
        if new:
            new_leaf_rank, new_leaf_name = taxon_hier.pop()
            if not taxon_hier:
                LOGGER.error(
                    f"{new_leaf_name} marked as new taxon with no additional taxonomy"
                )
                return

            LOGGER.info(f"{new_leaf_name} marked as new taxon")

        # Search for the taxon
        try:
            hr_taxon = self.validator.taxa_search(m_name, taxon_hier)
        except NCBIError:
            LOGGER.error("Search based on taxon hierarchy failed")
            return

        # Populate taxon names, index and higher taxon hierarchy
        self.taxon_names.add(m_name)

        if new:
            # Add the new taxon as the named taxon, linked to the hr_taxon as a parent
            # with a -1 ncbi ID and "user" status
            self.taxon_index.append(
                (
                    m_name,
                    -1,
                    hr_taxon.ncbi_id,
                    new_leaf_name,
                    new_leaf_rank,
                    "user",
                )
            )
        else:
            # Otherwise add the hr_taxon
            self.taxon_index.append(
                (
                    m_name,
                    hr_taxon.ncbi_id,
                    hr_taxon.parent_ncbi_id,
                    hr_taxon.name,
                    hr_taxon.rank,
                    "accepted",
                )
            )

            if hr_taxon.non_canon_name is not None:
                self.taxon_index.append(
                    (
                        m_name,
                        hr_taxon.ncbi_id,
                        hr_taxon.parent_ncbi_id,
                        hr_taxon.non_canon_name,
                        hr_taxon.rank,
                        hr_taxon.non_canon_name_class,
                    )
                )

        self.hierarchy.update(hr_taxon.taxa_hier)

    @loggerinfo_push_pop("Indexing taxonomic hierarchy")
    def index_higher_taxa(self) -> None:
        """Extend the taxon index to include higher taxa.

        This method uses the taxon hierarchy entries to add higher taxa to the taxon
        index for the instance. It does not duplicate taxa already in the index.
        """

        # Use the taxon hierarchy entries to add higher taxa
        # - drop taxa with a NCBI ID already in the index

        # Only add if taxon id is not already included in known
        known = [tx[1] for tx in self.taxon_index]
        to_add = [tx for tx in self.hierarchy if tx[2] not in known]
        to_add.sort(key=lambda val: NCBI_BACKBONE_RANKS.index(val[0]))

        # Look up the taxonomic hierarchy
        for tx_lev, tx_nme, tx_id, p_id in to_add:
            # Add all this to the taxonomy
            self.taxon_index.append((None, tx_id, p_id, tx_nme, tx_lev, "accepted"))
            LOGGER.info(f"Added {tx_lev} {tx_nme}")

    @property
    def is_empty(self) -> bool:
        """Check if an NCBITaxa instance contains any taxa."""
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
) -> str | tags.div:
    """Render a taxon index as text or html.

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
                if tx["taxon_rank"] in NCBI_BACKBONE_RANKS:
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
                elif tx["taxon_rank"] not in NCBI_BACKBONE_RANKS:
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


def construct_bi_or_tri(higher_nm: str, lower_nm: str, tri: bool) -> str:
    """Generate a species binomial or a subspecies trinomial.

    NCBI tools often return the separate components of species and subspecies names: for
    example, genus _Escherichia_ and species _coli_. However the NCBI database
    associates taxon ID with complete binomial and trinomial names. This function parses
    the provided inputs to try and construct those names, whilst also handling extra
    tags in binomials, such as 'candidatus', that are included by the NCBI for some
    taxa.

    Args:
        higher_nm: The NCBI genus/species name
        lower_nm: The NCBI species/subspecies name
        tri: Are we looking at a subspecies trinomial

    Raises:
        ValueError: where the input cannot be safely parsed into a binomial or trinomial
            name.
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
        value = higher_nm.strip() + " " + lower_nm.strip()
    # Look for Candidatus formats
    elif lower_nm.lower().startswith("candidatus") or higher_nm.lower().startswith(
        "candidatus"
    ):
        if lower_nm.lower().startswith("candidatus"):
            # Construct full name with first word of lower name removed
            nm = higher_nm.strip()
            for i in range(1, len(lower_nm.split())):
                nm = nm + " " + lower_nm.split()[i]
            value = nm
        else:
            value = higher_nm.strip() + " " + lower_nm.strip()
    # Catch too short species name case
    elif tri and len(higher_nm.split()) == 1:
        msg = f"Species name ({higher_nm}) too short"
        LOGGER.error(msg)
        raise ValueError(msg)

    # Then check that lower name is more words than the higher name
    elif len(lower_nm.split()) > len(higher_nm.split()):
        if lower_nm.lower().startswith(higher_nm.lower()):
            value = lower_nm
        else:
            msg = (
                f"{L} name ({lower_nm}) appears to be {type} but "
                f"does not contain {H.lower()} name ({higher_nm})"
            )
            LOGGER.error(msg)
            raise ValueError(msg)

    else:
        msg = f"{H} name ({higher_nm}) appears to be too long"
        LOGGER.error(msg)
        raise ValueError(msg)

    return value

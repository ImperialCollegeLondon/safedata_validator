"""This module provides the Locations class, which used to validate a set of known and or
new location information, formatted in the safedata_validator style. An instance can
then also be used to track which of the validated set of locations is used in the rest
of the dataset.
"""  # noqa D415

from openpyxl.worksheet.worksheet import Worksheet
from shapely import wkt
from shapely.errors import WKTReadingError

from safedata_validator.extent import Extent
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
    IsLocName,
    IsLower,
    IsNotBlank,
    IsNumber,
    IsString,
)


class Locations:
    """An interface for Location metadata.

    A Locations instance is initialised using a Resources instance that provides data on
    known valid locations. The instance validates location names provided in the
    Locations table of a Dataset and then validates and updates the latitudinal and
    logitudinal extent of those locations. The instance can then be used to track the
    use of location names across data tables in the dataset.

    Args:
        resources: A Resources instance, used to provide information about
            known locations
        latitudinal_extent: An Extent instance tracking latititudinal extents.
        longitudinal_extent: An Extent instance tracking longitudinal extents.

    Attributes:
        locations: A list of
        locations_index:
        locations_used:
        valid_locations:
        location_aliases:
        known_loc_names:
        latitudinal_extent:
        longitudinal_extent:
    """

    def __init__(
        self,
        resources: Resources,
        latitudinal_extent: Extent | None = None,
        longitudinal_extent: Extent | None = None,
    ) -> None:
        self.n_errors = 0
        self.locations: set = set()
        self.location_index: list = []
        self.locations_used: set = set()

        self.valid_locations = resources.valid_locations
        self.location_aliases = resources.location_aliases
        self.known_loc_names = set(
            list(resources.valid_locations.keys())
            + list(resources.location_aliases.keys())
        )

        # Attach or create extents
        if latitudinal_extent is None:
            self.latitudinal_extent = Extent(
                "latitudinal extent",
                (float, int),
                hard_bounds=resources.extents.latitudinal_hard_extent,
                soft_bounds=resources.extents.latitudinal_soft_extent,
            )
        else:
            self.latitudinal_extent = latitudinal_extent

        if longitudinal_extent is None:
            self.longitudinal_extent = Extent(
                "latitudinal extent",
                (float, int),
                hard_bounds=resources.extents.longitudinal_hard_extent,
                soft_bounds=resources.extents.longitudinal_soft_extent,
            )
        else:
            self.longitudinal_extent = longitudinal_extent

    @loggerinfo_push_pop("Loading Locations worksheet")
    def load(self, worksheet: Worksheet):
        """Populate a Locations instance from an Excel Worksheet.

        Validates the contents of a locations table stored in an Excel Worksheet and
        then updates the geographic extent of the locations used.

        Args:
            worksheet: An openpyxl Worksheet instance containing the formatted set of
                locations used within a Dataset.
        """
        handler = get_handler()
        start_errors = handler.counters["ERROR"]

        # Load the locations data frame - which runs header checks
        dframe = GetDataFrame(worksheet)

        if not dframe.data_columns:
            LOGGER.error("No data or only headers in Locations worksheet")
            return

        # Dupe headers likely cause serious issues, so stop
        if "duplicated" in dframe.bad_headers:
            LOGGER.error("Cannot parse locations with duplicated headers")
            return

        # Reduce to lower case
        # TODO - not trapping dupes that are only case
        #        differences. Do that in GetDataFrame?
        headers = IsLower(dframe.headers).values

        # Check location names are available
        if "location name" not in headers:
            LOGGER.error("Location name column not found")
            return

        # Get dictionaries of values for each row
        locs = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]

        # Split up old and new if there are there any new ones?
        if "new" in headers:
            # Check the New column is just yes, no
            new_vals_with_blanks = IsLower([rw["new"] for rw in locs])
            new_vals = IsNotBlank(new_vals_with_blanks, keep_failed=False)
            if not new_vals:
                LOGGER.error("Missing values in 'new' field")

            # check only yes or no entries
            valid_new = {"yes", "no"}
            bad_new = set(new_vals) - valid_new
            if bad_new:
                LOGGER.error(
                    "Values other than yes and no in 'new' field: ",
                    extra={"join": bad_new},
                )

            # Parse locations that can be assigned to new or known.
            known_locs = [
                rw
                for rw in locs
                if isinstance(rw["new"], str) and rw["new"].lower() == "no"
            ]
            new_locs = [
                rw
                for rw in locs
                if isinstance(rw["new"], str) and rw["new"].lower() == "yes"
            ]
        else:
            new_locs = []
            known_locs = locs

        if known_locs:
            known_loc_names = [lc["location name"] for lc in known_locs]
            self.add_known_locations(known_loc_names)

        if new_locs:
            self.add_new_locations(new_locs)

        # summary of processing
        self.n_errors = handler.counters["ERROR"] - start_errors

        if self.n_errors > 0:
            LOGGER.info(f"Locations contains {self.n_errors} errors")
        else:
            LOGGER.info(f"{len(self.locations)} locations loaded correctly")

    @loggerinfo_push_pop("Checking new locations")
    def add_new_locations(self, locs: list[dict]):
        """Add new locations to a Locations instance.

        This method takes a list of dictionaries giving the details of new locations to
        be added to the instance. These are user-defined locations that are not included
        in the set of known locations loaded from the instance resources. These
        dictionaries  should contain keys `location name` and `type` and then at least
        one of:

        * `latitude` _and_ `longitude` as float values
        * `wkt` providing a WellKnownText geometry for the location.

        Either of these options _can_ be 'NA' to show that the location coordinates are
        not known, but they must be provided.

        Args:
            locs: The list of dictionaries of user-defined locations.
        """

        # Validation - TODO check locs is a list of dicts

        # - Do all the dicts have the same keys
        loc_keys = {tuple(k.keys()) for k in locs}
        if len(loc_keys) > 1:
            LOGGER.critical("Inconsistent keys in add_new_locations")
            return

        # - Do they provide location names...
        location_keys = loc_keys.pop()
        if "location name" not in location_keys:
            LOGGER.error("No location name entries in add_new_locations")
            return

        # - ... and are any of those names blank ...
        loc_names_with_blanks = [itm["location name"] for itm in locs]
        loc_names = IsNotBlank(loc_names_with_blanks, keep_failed=False)
        if not loc_names:
            LOGGER.error("Location names contains empty cells or whitespace text")

        # - ... or not strings. Only allow strings for new locations - known
        #   locations get a pass for integer site codes but not here.
        loc_names_as_str = IsString(loc_names)
        if not loc_names_as_str:
            LOGGER.error(
                "New location names include non-string values: ",
                extra={"join": loc_names_as_str.failed},
            )

        # Look for duplicated names in inputs - this includes all types and so
        # could give messy information but can't detect duplicates late on
        # cleaned names because we're dealing with sets by that point and there
        # can be no duplication in sets
        dupes = HasDuplicates(loc_names_as_str)
        if dupes:
            LOGGER.error(
                "New location names contain duplicated values: ",
                extra={"join": dupes.duplicated},
            )

        # Look for new names that duplicate known names
        duplicates_existing = [
            rw["location name"]
            for rw in locs
            if rw["location name"] in self.known_loc_names
        ]

        if duplicates_existing:
            LOGGER.error(
                "New location names duplicate known names and aliases: ",
                extra={"join": duplicates_existing},
            )

        # Type is required - used to indicate the kind of location in the absence
        # of any actual geodata
        if "type" not in location_keys:
            LOGGER.error("New locations do not provide the location type")
            for this_loc in locs:
                this_loc["type"] = "MISSING"
        else:
            # get lowercase types
            geo_types_with_blanks = set(IsLower([vl["type"] for vl in locs]))

            # Handle blanks
            geo_types = IsNotBlank(geo_types_with_blanks, keep_failed=False)
            if not geo_types:
                LOGGER.error(
                    "Types for new locations contains blank or whitespace entries."
                )

            # Handle unknown geo types
            valid_geotypes = {"point", "linestring", "polygon", "transect", "area"}
            bad_geo_types = set(geo_types) - valid_geotypes
            if bad_geo_types:
                LOGGER.error(
                    "New locations include unknown location types: ",
                    extra={"join": bad_geo_types},
                )

        # Geometry information
        # Record which geom columns are present
        if ("latitude" in location_keys) ^ ("longitude" in location_keys):
            LOGGER.error(
                "New locations should either latitude _and_ longitude or neither"
            )

        lonlat_provided = "latitude" in location_keys and "longitude" in location_keys
        wkt_provided = "wkt" in location_keys

        # TODO - supplying both is not an error, and probably shouldn't be. WKT
        # takes priority in the index when both are present, but not testing for
        # congruence at present

        # Look for geographic data (even if it is just the explicit statement
        # that none is available using NA)
        if not (lonlat_provided or wkt_provided):
            LOGGER.error(
                "New locations reported: you must provide Lat/Long or WKT,"
                "using NA explicitly when this data is missing."
            )
        else:
            # Check Lat Long and WKT

            if lonlat_provided:
                LOGGER.info("Validating lat / long data")
                FORMATTER.push()
                for axs, ext_attr in [
                    ("latitude", "latitudinal_extent"),
                    ("longitude", "longitudinal_extent"),
                ]:
                    # Allow NAs for unknown location points
                    axs_vals_with_blanks = [vl[axs] for vl in locs if vl[axs] != "NA"]
                    axs_vals = IsNotBlank(axs_vals_with_blanks, keep_failed=False)
                    if not axs_vals:
                        LOGGER.error(f"Blank {axs} values for new locations: use NA.")

                    # Check for data types _here_ to keep interpretable errors
                    # TODO - maybe simplify Extent objects to _assume_ types.
                    #        Probably no but there is some duplication of effort here
                    axs_vals_as_number = IsNumber(axs_vals, keep_failed=False)
                    if not axs_vals_as_number:
                        LOGGER.error(
                            f"Non-numeric {axs} values for new locations: ",
                            extra={"join": axs_vals_as_number.failed},
                        )

                    # Update extent instances
                    if axs_vals_as_number.values:
                        ext = getattr(self, ext_attr)
                        ext.update(axs_vals_as_number)

                FORMATTER.pop()

            if wkt_provided:
                LOGGER.info("Validating WKT data")
                FORMATTER.push()

                blank_wkt = []
                non_string_wkt = []
                bad_wkt = []
                bounds = []

                for this_new_loc in locs:
                    if this_new_loc["wkt"] is None:
                        blank_wkt.append(this_new_loc["location name"])
                    elif not isinstance(this_new_loc["wkt"], str):
                        non_string_wkt.append(this_new_loc["location name"])
                    elif this_new_loc["wkt"].isspace():
                        blank_wkt.append(this_new_loc["location name"])
                    elif this_new_loc["wkt"] == "NA":
                        pass

                    else:
                        # Run the potential WKT through the parser
                        try:
                            this_new_geom = wkt.loads(this_new_loc["wkt"])
                        except WKTReadingError:
                            bad_wkt.append(this_new_loc["location name"])
                        else:
                            # Is it a valid 2D geom
                            if not this_new_geom.is_valid or this_new_geom.has_z:
                                bad_wkt.append(this_new_loc["location name"])
                            # Store the extents to check for sensible coordinates
                            bounds.append(this_new_geom.bounds)

                if blank_wkt:
                    LOGGER.error("Blank WKT values for new locations: use NA.")

                if non_string_wkt:
                    LOGGER.error(
                        "WKT values for new location not a string: ",
                        extra={"join": non_string_wkt},
                    )

                if bad_wkt:
                    LOGGER.error(
                        "WKT information badly formatted, not geometrically valid or "
                        "3D: ",
                        extra={"join": bad_wkt},
                    )

                if bounds:
                    # Extract from bound tuples to lists of lats and longs
                    lat_bnds = [bnd[1] for bnd in bounds] + [bnd[3] for bnd in bounds]
                    lng_bnds = [bnd[0] for bnd in bounds] + [bnd[2] for bnd in bounds]
                    self.latitudinal_extent.update(lat_bnds)
                    self.longitudinal_extent.update(lng_bnds)

                FORMATTER.pop()

        # new location names
        # - test for duplicated names to already added values
        duped_names = [lc for lc in loc_names if lc in self.locations]
        if duped_names:
            LOGGER.error(
                "Location names already added to Location instance: ",
                extra={"join": duped_names},
            )

        self.locations.update(loc_names)

        for this_new_loc in locs:
            if wkt_provided and this_new_loc["wkt"] != "NA":
                geom = this_new_loc["wkt"]
            elif lonlat_provided and (
                (this_new_loc["latitude"] != "NA")
                and (this_new_loc["longitude"] != "NA")
            ):
                geom = "Point({longitude} {latitude})".format(**this_new_loc)
            else:
                geom = None

            self.location_index.append((this_new_loc["location name"], True, geom))

    @loggerinfo_push_pop("Checking known locations")
    def add_known_locations(self, loc_names: list):
        """Add known locations to a Locations instance.

        This method takes a list of values and tries to validate those values against
        known locations from the loaded resources. The values are expected to be
        strings.

        Args:
            loc_names: A list of known location names.
        """

        # Check for blanks
        loc_names_no_blanks = IsNotBlank(loc_names, keep_failed=False)
        if not loc_names_no_blanks:
            LOGGER.error("Location names contains empty cells or whitespace text")

        # Look for duplicated values in names - this includes all types and so
        # could give messy information but can't detect duplicates late on
        # cleaned names because we're dealing with sets by that point and there
        # can be no duplication
        dupes = HasDuplicates(loc_names_no_blanks)
        if dupes:
            LOGGER.error(
                "Added names contain duplicated values: ",
                extra={"join": dupes.duplicated},
            )

        # Validate and standardise types - strings or integer codes.
        loc_names_standardised = IsLocName(loc_names_no_blanks, keep_failed=False)

        if not loc_names_standardised:
            LOGGER.error(
                "Location names contains values that are not strings or integers: ",
                extra={"join": loc_names_standardised.failed},
            )

        # Enforce strings and check loc names exist
        loc_names_as_str = {str(v) for v in loc_names_standardised}
        unknown = loc_names_as_str - self.known_loc_names
        if unknown:
            LOGGER.error(
                "Unknown locations found: ",
                extra={
                    "join": unknown,
                },
            )

        # are aliases being used?
        aliased = loc_names_as_str & set(self.location_aliases.keys())
        if aliased:
            LOGGER.warning(
                "Locations aliases used. Maybe change to primary location names: ",
                extra={"join": aliased},
            )

        # Get the bounding box of known locations and aliased locations
        bbox_keys = (loc_names_as_str - (unknown | aliased)) | {
            self.location_aliases[ky] for ky in aliased
        }

        # get the extents of known unaliased locations
        if bbox_keys:
            bbox = [
                vl for ky, vl in list(self.valid_locations.items()) if ky in bbox_keys
            ]
            bbox = list(zip(*bbox))
            self.longitudinal_extent.update((min(bbox[0]), max(bbox[2])))
            self.latitudinal_extent.update((min(bbox[1]), max(bbox[3])))

        # Update location names and index
        # - test for duplicated names to already added values
        duped_names = [lc for lc in loc_names_as_str if lc in self.locations]
        if duped_names:
            LOGGER.error(
                "Location names already added to Location instance: ",
                extra={"join": duped_names},
            )

        self.locations.update(loc_names_as_str)
        index_entries = [(lc, False, None) for lc in loc_names_as_str]
        self.location_index.extend(index_entries)

    @property
    def is_empty(self) -> bool:
        """Reports if any locations have been loaded in a Locations instance."""
        return len(self.locations) == 0

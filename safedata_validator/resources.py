"""Load and check validation resources.

The `safedata_validator` package needs access to some local resources and
configuration to work. The three main resources for file validation are:

-   locations: A path to a locations data file, providing a gazetteer of known
    locations and their details.

-   gbif_database: The path to a local SQLite copy of the GBIF backbone database.

-   ncbi_database: The path to a local SQLite copy of the NCBI database files.

The [Resources][safedata_validator.resources.Resources] class is used to locate
and validate these resources, and then provide those validated resources to
other components of the package.

A configuration file can be passed as `cfg_path` when creating an instance,
but if no arguments are provided then an attempt is made to find and load
configuration files in the user and then site config locations defined by the
`appdirs` package. See [here][usage/usage] for details.

"""


import contextlib
import os
import sqlite3
from datetime import date
from typing import Union

import appdirs
import simplejson
from configobj import ConfigObj, flatten_errors
from dateutil.parser import isoparse
from dotmap import DotMap
from simplejson.errors import JSONDecodeError
from validate import Validator, VdtParamError, VdtValueError, is_list

from safedata_validator.logger import (
    FORMATTER,
    LOGGER,
    log_and_raise,
    loggerinfo_push_pop,
)

CONFIGSPEC = {
    "locations": "string()",
    "gbif_database": "string()",
    "ncbi_database": "string()",
    "extents": {
        "temporal_soft_extent": "date_list(min=2, max=2, default=None)",
        "temporal_hard_extent": "date_list(min=2, max=2, default=None)",
        "latitudinal_hard_extent": "float_list(min=2, max=2, default=list(-90, 90))",
        "latitudinal_soft_extent": "float_list(min=2, max=2, default=None)",
        "longitudinal_hard_extent": "float_list(min=2, max=2, default=list(-180, 180))",
        "longitudinal_soft_extent": "float_list(min=2, max=2, default=None)",
    },
    "zenodo": {
        "community_name": "string(default=safe)",
        "use_sandbox": "boolean(default=None)",
        "zenodo_sandbox_api": "string(default=None)",
        "zenodo_sandbox_token": "string(default=None)",
        "zenodo_api": "string(default=None)",
        "zenodo_token": "string(default=None)",
        "contact_name": "string(default=None)",
        "contact_affiliation": "string(default=None)",
        "contact_orcid": "string(default=None)",
    },
    "metadata": {"api": "string(default=None)", "token": "string(default=None)"},
}
"""dict: The safedata_validator package use the `configobj.ConfigObj`
package to handle resource configuration. This dict defines the basic expected
specification for the configuration and allows the ConfigObj.validate() method
to do basic validation and type conversions.
"""


def date_list(value: str, min: str, max: str) -> list[date]:
    """Validate config date lists.

    A configobj.Validator extension function to check configuration values
    containing a list of ISO formatted date strings and to return parsed
    values.

    Args:
        value: A string containing comma-separated ISO date strings
        min: The minimum allowed number of entries
        max: The maximum allowed number of entries
    """
    # min and max are supplied as a string, test conversion to int
    try:
        min_int = int(min)
    except ValueError:
        raise VdtParamError("min", min)
    try:
        max_int = int(max)
    except ValueError:
        raise VdtParamError("max", max)

    # Check the supplied value is a list, triggering any issues
    # with list formatting
    value = is_list(value, min=min_int, max=max_int)

    # Next, check every member in the list is an ISO date string
    # noting that this strips out time information
    out = []
    for entry in value:

        try:
            parsed_entry = isoparse(entry).date()
        except ValueError:
            raise VdtValueError(entry)

        out.append(parsed_entry)

    # Return parse values
    return out


@loggerinfo_push_pop("Configuring Resources")
class Resources:
    """Load and check validation resources.

    Creating an instance of this class locates and validate resources for using the
    `safedata_validator` package, either from the provided configuration details or from
    the user and then site config locations defined by the appdirs package.

    Args:
        config:
            A path to a configuration file, or a dict or list providing package
            configuration details. The list format should provide a list of strings,
            each representing a line in the configuration file. The dict format is a
            dictionary with the required nested dictionary structure and values.

    Attributes:
        config_type: The method used to specify the resources. One of
            'init_dict', 'init_list', 'init_file', 'user_config' or 'site_config'.
        locations: The path to the locations file
        gbif_database: The path to the GBIF database file
        ncbi_database: The path to the NCBI database file
        valid_locations: The locations defined in the locations file
        location_aliases: Location aliases defined in the locations file
        extents: A DotMap of extent data
        zenodo: A DotMap of Zenodo information
    """

    def __init__(self, config: Union[str, list, dict] = None) -> None:

        # User and site config paths
        user_cfg_file = os.path.join(
            appdirs.user_config_dir(), "safedata_validator", "safedata_validator.cfg"
        )
        site_cfg_file = os.path.join(
            appdirs.site_config_dir(), "safedata_validator", "safedata_validator.cfg"
        )

        # First try and populate from a config file.
        if config is not None:
            if isinstance(config, str):
                config_type = "init file"
            elif isinstance(config, list):
                config_type = "init list"
            elif isinstance(config, dict):
                config_type = "init dict"
        elif os.path.exists(user_cfg_file) and os.path.isfile(user_cfg_file):
            config = user_cfg_file
            config_type = "user file"
        elif os.path.exists(site_cfg_file) and os.path.isfile(site_cfg_file):
            config = site_cfg_file
            config_type = "site file"
        else:
            LOGGER.critical(f"No user config in {user_cfg_file}")
            LOGGER.critical(f"No site config in {site_cfg_file}")
            log_and_raise("No config files provided or found", RuntimeError)
            return

        # Report resource config location and type
        msg = f"Configuring resources from {config_type}"
        if "file" in config_type:
            msg += f": {config}"
        LOGGER.info(msg)

        # Try and load the found configuration
        config_loaded = self._load_config(config, config_type)

        # Set attributes -
        # HACK - this now seems clumsy - the ConfigObj instance is already a
        #        class containing the config attributes. Having a _function_
        #        that returns a modified ConfigObj instance seems more direct
        #        than having to patch this list of attributes.
        self.locations = config_loaded.locations
        self.gbif_database = config_loaded.gbif_database
        self.ncbi_database = config_loaded.ncbi_database
        self.config_type = config_loaded.config_type
        self.config_source = config_loaded.config_source

        self.extents = config_loaded.extents
        self.zenodo = config_loaded.zenodo
        self.metadata = config_loaded.metadata

        self.gbif_timestamp = None
        self.ncbi_timestamp = None

        # Valid locations is a dictionary keying string location names to tuples of
        # floats describing the location bounding box
        self.valid_location: dict[str, list[float]] = dict()
        # Location aliases is a dictionary keying a string to a key in valid locations
        self.location_aliases: dict[str, str] = dict()

        # Validate the resources
        self._validate_locations()
        self._validate_gbif()
        self._validate_ncbi()

    @staticmethod
    def _load_config(config: Union[str, list, dict], cfg_type: str) -> DotMap:
        """Load a configuration file.

        This private static method attempts to load a JSON configuration file
        from a path.

        Args:
            config: Passed from Resources.__init__()
            cfg_type: Identifies the route used to provide the configuration details

        Returns:
             If the file does not exist, the function returns None. Otherwise,
             it returns a DotMap of config parameters.
        """

        # Otherwise, there is a file, so try and use it and now raise if there
        # is a problem: don't skip over broken resource configurations.
        # - First, create a validator instance that handles lists of dates
        cf_validator = Validator({"date_list": date_list})

        # - Now load the config input and then apply the basic validation - are
        #   the values of the right type, right count etc.
        config_obj = ConfigObj(config, configspec=CONFIGSPEC)
        valid = config_obj.validate(cf_validator, preserve_errors=True)

        # If there are config file issues, then bail out.
        if isinstance(valid, dict):
            LOGGER.critical("Configuration issues: ")
            FORMATTER.push()
            for sec, key, err in flatten_errors(config_obj, valid):
                sec.append(key)
                LOGGER.critical(f"In config '{'.'.join(sec)}': {err}")
            FORMATTER.pop()
            raise RuntimeError("Configuration failure")

        # convert to a DotMap for ease
        config_obj = DotMap(config_obj)

        return config_obj

    def _validate_locations(self) -> None:
        """Validate and load a locations file.

        This private function checks whether a locations path: exists, is a JSON file,
        and contains location and alias data. It populates the instance attributes
        """

        if self.locations is None or self.locations == "":
            log_and_raise("Locations file missing in configuration", RuntimeError)

        LOGGER.info(f"Validating locations: {self.locations}")

        # Now check to see whether the locations file behaves as expected
        if not os.path.exists(self.locations) and not os.path.isfile(self.locations):
            log_and_raise("Local locations file not found", OSError)

        try:
            loc_payload = simplejson.load(open(self.locations, mode="r"))
        except (JSONDecodeError, UnicodeDecodeError):
            log_and_raise("Local locations file not JSON encoded.", OSError)

        # process the locations payload
        if {"locations", "aliases"} != loc_payload.keys():
            log_and_raise("Locations data malformed", RuntimeError)

        self.valid_locations = loc_payload["locations"]
        self.location_aliases = loc_payload["aliases"]

    def _validate_gbif(self) -> None:
        """Validate the GBIF settings.

        This private function checks whether to use the online API or a local backbone
        database and then validates the provided sqlite3 database file. It updates the
        instance with validated details.
        """

        self.gbif_timestamp = validate_taxon_db(
            self.gbif_database, "GBIF", ["backbone"]
        )

    def _validate_ncbi(self) -> None:
        """Validate the NCBI settings.

        This private function checks whether to use the online API or a local backbone
        database and then validates the provided sqlite3 database files. It updates the
        instance with validated details.
        """

        self.ncbi_timestamp = validate_taxon_db(
            self.ncbi_database, "NCBI", ["nodes", "merge", "names"]
        )


def validate_taxon_db(db_path: str, db_name: str, tables: list[str]) -> str:
    """Validate a local taxon database file.

    This helper function validates that a given path contains a valid taxonomy database:

    - the required tables are all present, automatically including the timestamp table.
    - the timestamp table contains a single ISO format date showing the database
      version.

    Args:
        db_path: Location of the SQLite3 database.
        db_name: A label for the taxonomy database - used in logger messages.
        tables: A list of table names expected to be present in the database.

    Returns:
        The database timestamp as an ISO formatted date string.
    """

    LOGGER.info(f"Validating local {db_name} database: {db_path}")

    # Does the provided path exist and is it a functional SQLite database
    # with a backbone table? Because sqlite3 can 'connect' to any path,
    # use a query attempt to reveal exceptions

    if not os.path.exists(db_path):
        log_and_raise(f"Local {db_name} database not found", OSError)

    # Connect to the file (which might or might not be a database containing the
    # required tables)
    with contextlib.closing(sqlite3.connect(db_path)) as conn:

        # Check that it is a database by running a query
        try:
            db_tables = conn.execute(
                "SELECT name FROM sqlite_master WHERE type ='table';"
            )
        except sqlite3.DatabaseError:
            log_and_raise(f"Local {db_name} database not an SQLite3 file.", OSError)

        # Check the required tables against found tables
        db_tables = set([rw[0] for rw in db_tables.fetchall()])
        required_tables = set(tables + ["timestamp"])
        missing = required_tables.difference(db_tables)

        if missing:
            log_and_raise(
                f"Local {db_name} database does not contain required tables: ",
                RuntimeError,
                extra={"join": missing},
            )

        # Check the timestamp table contains a single ISO date
        cursor = conn.execute("select * from timestamp;")
        timestamp = cursor.fetchall()

    # Is there one unique date in the table
    if len(timestamp) != 1:
        log_and_raise(
            f"Local {db_name} database timestamp table contains more than one entry.",
            RuntimeError,
        )

    try:
        # Extract first entry in first row
        timestamp = timestamp[0][0]
        isoparse(timestamp)
    except ValueError:
        log_and_raise(
            f"Local {db_name} database timestamp value is not an ISO date.",
            RuntimeError,
        )

    return timestamp

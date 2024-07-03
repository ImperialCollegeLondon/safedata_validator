"""The `safedata_validator` package needs access to some local resources and
configuration to work. The core resources for file validation are:

- gazetteer: A path to a GeoJSON formatted gazetteer of known locations and their
    details.

- location_aliases: A path to a CSV file containing known aliases of the location
    names provided in the gazetteer.

- gbif_database: The path to a local SQLite copy of the GBIF backbone database.

- ncbi_database: The path to a local SQLite copy of the NCBI database files.

- project_database: Optionally, a path to a CSV file providing valid project IDs.

The [Resources][safedata_validator.resources.Resources] class is used to locate and
validate these resources, and then provide those validated resources to other components
of the package.

A configuration file can be passed as `config` when creating an instance, but if no
arguments are provided then an attempt is made to find and load configuration files in
the user and then site config locations defined by the `appdirs` package. See
[here](../../data_managers/install/configuration.md#configuration-file-locations) for
details.
"""  # noqa D415

import contextlib
import os
import sqlite3
from csv import DictReader
from csv import Error as csvError
from datetime import date

import appdirs
import simplejson
from configobj import ConfigObj, flatten_errors
from dateutil.parser import isoparse
from dotmap import DotMap
from shapely.geometry import shape
from simplejson.errors import JSONDecodeError
from validate import Validator, VdtParamError, VdtValueError, is_list

from safedata_validator.logger import (
    FORMATTER,
    LOGGER,
    log_and_raise,
    loggerinfo_push_pop,
)

CONFIGSPEC = {
    "gazetteer": "string()",
    "location_aliases": "string()",
    "gbif_database": "string()",
    "ncbi_database": "string()",
    "project_database": "string(default=None)",
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
        "zenodo_sandbox_token": "string(default=None)",
        "zenodo_token": "string(default=None)",
        "contact_name": "string(default=None)",
        "contact_affiliation": "string(default=None)",
        "contact_orcid": "string(default=None)",
        "html_template": "string(default=None)",
    },
    "metadata": {
        "api": "string(default=None)",
        "token": "string(default=None)",
        "ssl_verify": "boolean(default=True)",
    },
    "xml": {
        # Additional elements required for XML
        "languageCode": "string(default=None)",
        "characterSet": "string(default=None)",
        "contactCountry": "string(default=None)",
        "contactEmail": "string(default=None)",
        "epsgCode": "integer(default=4326)",
        "projectURL": "string(default=None)",
        "topicCategories": "string_list(default=None)",
        "lineageStatement": "string(default=None)",
    },
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
    `safedata_validator` package. The resources can be located in several ways, which
    use the following order of priority:

    * configuration details provided directly via the ``config`` argument (see below),
    * a path to a configuration file set in the ``SAFEDATA_VALIDATOR_CONFIG``
      environment variable,
    * a configuration file in the standard user location, or
    * a configuration file in the standard system wide location.

    The standard locations follow the implementation of the ``appdirs`` package.
    Typically, end users will rely on the last two options, but the first two options
    are useful for testing and validation.

    Args:
        config:
            A path to a configuration file, or a dict or list providing package
            configuration details. The list format should provide a list of strings,
            each representing a line in the configuration file. The dict format is a
            dictionary with the required nested dictionary structure and values.

    Attributes:
        config_type: The method used to specify the resources. One of
            'init_dict', 'init_list', 'init_path', 'env_var_path', 'user_path' or
            'site_path'.
        gazetteer: The path to the gazetteer file
        location_aliases: The path to the location_aliases file
        gbif_database: The path to the GBIF database file
        ncbi_database: The path to the NCBI database file
        project_database: An optional path to a database of valid project IDs
        valid_locations: The locations defined in the locations file
        location_aliases: Location aliases defined in the locations file
        extents: A DotMap of extent data
        zenodo: A DotMap of Zenodo information
    """

    def __init__(self, config: str | list | dict | None = None) -> None:
        # Get the standard user and site config paths for the platform
        user_cfg_file = os.path.join(
            appdirs.user_config_dir(), "safedata_validator", "safedata_validator.cfg"
        )
        site_cfg_file = os.path.join(
            appdirs.site_config_dir(), "safedata_validator", "safedata_validator.cfg"
        )

        # Look for a config path as an environment variable
        config_env_path = os.getenv("SAFEDATA_VALIDATOR_CONFIG")

        # Now resolve what to use in order of priority
        if config is not None:
            if isinstance(config, str):
                if os.path.exists(config) and os.path.isfile(config):
                    config_type = "init_path"
                else:
                    log_and_raise(f"Config file path not found: {config}", RuntimeError)
                    return
            elif isinstance(config, list):
                config_type = "init_list"
            elif isinstance(config, dict):
                config_type = "init_dict"
        elif config_env_path is not None:
            config = config_env_path
            config_type = "env_var_path"
        elif os.path.exists(user_cfg_file) and os.path.isfile(user_cfg_file):
            config = user_cfg_file
            config_type = "user_path"
        elif os.path.exists(site_cfg_file) and os.path.isfile(site_cfg_file):
            config = site_cfg_file
            config_type = "site_path"
        else:
            LOGGER.critical(f"No user config in {user_cfg_file}")
            LOGGER.critical(f"No site config in {site_cfg_file}")
            log_and_raise("No config files provided or found", RuntimeError)
            return

        # Report resource config location and type
        msg = f"Configuring resources from {config_type}"
        if config_type.endswith("path"):
            msg += f": {config}"
        LOGGER.info(msg)

        # Try and load the found configuration
        config_loaded = self._load_config(config, config_type)

        # Set attributes -
        # HACK - this now seems clumsy - the ConfigObj instance is already a
        #        class containing the config attributes. Having a _function_
        #        that returns a modified ConfigObj instance seems more direct
        #        than having to patch this list of attributes.
        self.gaz_path = config_loaded.gazetteer
        self.localias_path = config_loaded.location_aliases
        self.gbif_database = config_loaded.gbif_database
        self.ncbi_database = config_loaded.ncbi_database
        self.project_database = (
            None
            if config_loaded.project_database == ""
            else config_loaded.project_database
        )
        self.config_type = config_loaded.config_type
        self.config_source = config_loaded.config_source

        self.extents = config_loaded.extents
        self.zenodo = config_loaded.zenodo
        self.metadata = config_loaded.metadata
        self.xml = config_loaded.xml

        self.gbif_timestamp: str | None = None
        self.ncbi_timestamp: str | None = None

        # Valid locations is a dictionary keying string location names to tuples of
        # floats describing the location bounding box
        self.valid_location: dict[str, list[float]] = dict()
        # Location aliases is a dictionary keying a string to a key in valid locations
        self.location_aliases: dict[str, str] = dict()
        # Projects are a dictionary keying project ID to a title.
        self.projects: dict[int, str] = dict()

        # Validate the resources
        self._validate_gazetteer()
        self._validate_location_aliases()
        self._validate_gbif()
        self._validate_ncbi()
        self._validate_projects()

    @staticmethod
    def _load_config(config: str | list | dict, cfg_type: str) -> DotMap:
        """Load a configuration file.

        This private static method attempts to load a JSON configuration file
        from a path.

        Args:
            config: Passed from Resources.__init__()
            cfg_type: Identifies the route used to provide the configuration details

        Raises:
            RuntimeError: If the file does not exist, or has issues.

        Returns:
            Returns a DotMap of config parameters.
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

    def _validate_gazetteer(self) -> None:
        """Validate and load a gazetteer file.

        This private function checks whether a gazetteer path: exists, is a JSON file,
        and contains location GeoJSON data. It populates the instance attributes
        """

        if self.gaz_path is None or self.gaz_path == "":
            log_and_raise("Gazetteer file missing in configuration", RuntimeError)

        LOGGER.info(f"Validating gazetteer: {self.gaz_path}")

        # Now check to see whether the locations file behaves as expected
        if not os.path.exists(self.gaz_path) and not os.path.isfile(self.gaz_path):
            log_and_raise("Gazetteer file not found", OSError)

        try:
            loc_payload = simplejson.load(open(self.gaz_path))
        except (JSONDecodeError, UnicodeDecodeError):
            log_and_raise("Gazetteer file not valid JSON", OSError)

        # Simple test for GeoJSON
        if (
            loc_payload.get("type") is None
            or loc_payload["type"] != "FeatureCollection"
        ):
            log_and_raise(
                "Gazetteer data not a GeoJSON Feature Collection", RuntimeError
            )

        try:
            self.valid_locations = {
                ft["properties"]["location"]: shape(ft["geometry"]).bounds
                for ft in loc_payload["features"]
            }
        except KeyError:
            log_and_raise(
                "Missing or incomplete location properties for gazetteer features",
                RuntimeError,
            )

    def _validate_location_aliases(self) -> None:
        """Validate and load location aliases.

        This private function checks whether a location_aliases path: exists, is a CSV
        file, and contains location_alias data. It populates the instance attributes
        """

        if self.localias_path is None or self.localias_path == "":
            log_and_raise(
                "Location aliases file missing in configuration", RuntimeError
            )

        LOGGER.info(f"Validating location aliases: {self.localias_path}")

        # Now check to see whether the locations file behaves as expected
        try:
            dictr = DictReader(open(self.localias_path))
        except FileNotFoundError:
            log_and_raise("Location aliases file not found", FileNotFoundError)
        except IsADirectoryError:
            log_and_raise("Location aliases path is a directory", IsADirectoryError)

        # Simple test for structure - field names only parsed when called, and this can
        # throw errors with bad file formats.
        try:
            if not dictr.fieldnames:
                log_and_raise("Location aliases file is empty", ValueError)
            else:
                fieldnames = set(dictr.fieldnames)
        except (UnicodeDecodeError, csvError):
            log_and_raise(
                "Location aliases file not readable as a CSV file with valid headers",
                ValueError,
            )

        if fieldnames != {"zenodo_record_id", "location", "alias"}:
            log_and_raise(
                "Location aliases file not readable as a CSV file with valid headers",
                ValueError,
            )

        # TODO - zenodo_record_id not being used here.
        self.location_aliases = {la["alias"]: la["location"] for la in dictr}

    def _validate_gbif(self) -> None:
        """Validate the GBIF settings.

        This private function validates the provided sqlite3 database file and updates
        the instance with validated details.
        """

        self.gbif_timestamp = validate_taxon_db(
            self.gbif_database, "GBIF", ["backbone"]
        )

    def _validate_ncbi(self) -> None:
        """Validate the NCBI settings.

        This private function validates the provided sqlite3 database files and updates
        the instance with validated details.
        """

        self.ncbi_timestamp = validate_taxon_db(
            self.ncbi_database, "NCBI", ["nodes", "merge", "names"]
        )

    def _validate_projects(self) -> None:
        """Validate and load a project database.

        This private function checks whether a project_database path: exists, is a CSV
        file, and contains project data. It populates the instance ``project_id``
        attribute.
        """

        if self.project_database is None:
            LOGGER.info("Configuration does not use project IDs.")
            return

        LOGGER.info(f"Validating project database: {self.project_database}")

        # Now check to see whether the project database behaves as expected
        try:
            dictr = DictReader(open(self.project_database, encoding="UTF-8"))
        except FileNotFoundError:
            log_and_raise("Project database file not found", FileNotFoundError)
        except IsADirectoryError:
            log_and_raise("Project database path is a directory", IsADirectoryError)

        # Simple test for structure - field names only parsed when called, and this can
        # throw errors with bad file formats.
        try:
            if not dictr.fieldnames:
                log_and_raise("Project database file is empty", ValueError)
            else:
                fieldnames = set(dictr.fieldnames)
        except (UnicodeDecodeError, csvError) as excep:
            LOGGER.critical(
                "Project database file not readable as a CSV file with valid headers"
            )
            raise excep

        required_names = {"project_id", "title"}
        if required_names.intersection(fieldnames) != required_names:
            log_and_raise(
                "Project database file does not contain project_id and title headers.",
                ValueError,
            )

        # Load the valid project ids
        try:
            self.projects = {int(prj["project_id"]): str(prj["title"]) for prj in dictr}
        except ValueError:
            log_and_raise(
                "Project database file values not integer IDs and text titles.",
                ValueError,
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

    LOGGER.info(f"Validating {db_name} database: {db_path}")

    if db_path is None or db_path == "":
        log_and_raise(f"{db_name} database not set in configuration", ValueError)

    # Does the provided path exist and is it a functional SQLite database
    # with a backbone table? Because sqlite3 can 'connect' to any path,
    # use a query attempt to reveal exceptions

    if not os.path.exists(db_path):
        log_and_raise(f"{db_name} database not found", FileNotFoundError)

    # Connect to the file (which might or might not be a database containing the
    # required tables)
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        # Check that it is a database by running a query
        try:
            db_tables = conn.execute(
                "SELECT name FROM sqlite_master WHERE type ='table';"
            )
        except sqlite3.DatabaseError:
            log_and_raise(f"Local {db_name} database not an SQLite3 file", ValueError)

        # Check the required tables against found tables
        db_tables_set = {rw[0] for rw in db_tables.fetchall()}
        required_tables = set([*tables, "timestamp"])
        missing = required_tables.difference(db_tables_set)

        if missing:
            log_and_raise(
                f"Local {db_name} database does not contain required tables: ",
                ValueError,
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
        timestamp_entry = timestamp[0][0]
        isoparse(timestamp_entry)
    except ValueError:
        log_and_raise(
            f"Local {db_name} database timestamp value is not an ISO date.",
            RuntimeError,
        )

    return timestamp_entry

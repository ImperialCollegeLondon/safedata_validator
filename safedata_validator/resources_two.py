"""The `safedata_validator` package needs access to some local resources and
configuration to work. The core resources for file validation are:

- gazetteer: A path to a GeoJSON formatted gazetteer of known locations and their
    details.

- location_aliases: A path to a CSV file containing known aliases of the location
    names provided in the gazetteer.

- gbif_database: The path to a local SQLite copy of the GBIF backbone database.

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
import tomllib
from csv import DictReader
from csv import Error as csvError
from datetime import date
from pathlib import Path
from typing import Self

import appdirs
import simplejson
from dateutil.parser import isoparse
from pydantic import (
    EmailStr,
    FilePath,
    HttpUrl,
    model_validator,
)
from shapely.geometry import shape
from simplejson.errors import JSONDecodeError

from safedata_validator.logger import (
    LOGGER,
    log_and_raise,
    loggerinfo_push_pop,
)
from safedata_validator.models import ForbiddenExtra


class Extents(ForbiddenExtra):
    temporal_hard_extent: tuple[date, date] | None = None
    temporal_soft_extent: tuple[date, date] | None = None
    latitudinal_hard_extent: tuple[float, float] | None = (-90, 90)
    latitudinal_soft_extent: tuple[float, float] | None = None
    longitudinal_hard_extent: tuple[float, float] | None = (-180, 180)
    longitudinal_soft_extent: tuple[float, float] | None = None


class Zenodo(ForbiddenExtra):
    community_name: str
    use_sandbox: bool
    zenodo_sandbox_token: str
    zenodo_token: str
    contact_name: str
    contact_affiliation: str
    contact_orcid: str
    project_url: HttpUrl
    html_template: FilePath | None = None


class Metadata(ForbiddenExtra):
    api: HttpUrl
    token: str
    ssl_verify: bool = True


class GeminiXML(ForbiddenExtra):
    languageCode: str
    characterSet: str
    contactCountry: str
    contactEmail: EmailStr
    epsgCode: int = 4326
    topicCategories: list[str]
    lineageStatement: str


class ResourcesPydantic(ForbiddenExtra):
    gazetteer: FilePath
    location_aliases: FilePath
    gbif_database: FilePath
    project_database: FilePath | None = None
    maximum_embargo_months: int = 24
    extents: Extents
    zenodo: Zenodo
    metadata: Metadata
    xml: GeminiXML

    # All of these fields are populated by model validation - this is a bit of hack,
    # technically users could try and populate these from the config but they will be
    # overwritten by the model validator functions
    valid_locations: dict = {}
    location_aliases: dict = {}
    project: dict = {}

    # @field_validator("gazetteer")
    # @classmethod
    @model_validator(mode="after")
    def validate_gazetteer(self) -> Self:
        """Validate and load a gazetteer file.

        This private function checks whether a gazetteer path: exists, is a JSON file,
        and contains location GeoJSON data. It populates the instance attributes
        """

        LOGGER.info(f"Validating gazetteer: {self.gazetteer}")

        try:
            loc_payload = simplejson.load(open(self.gazetteer))
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

        return self

    @model_validator(mode="after")
    def validate_location_aliases(self) -> Self:
        """Validate and load location aliases.

        This private function checks whether a location_aliases path: exists, is a CSV
        file, and contains location_alias data. It populates the instance attributes
        """

        LOGGER.info(f"Validating location aliases: {self.location_aliases}")

        # Now check to see whether the locations file behaves as expected
        try:
            dictr = DictReader(open(self.location_aliases))
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

        return self

    @model_validator(mode="after")
    def validate_gbif(self) -> Self:
        """Validate the GBIF settings.

        This private function validates the provided sqlite3 database file and updates
        the instance with validated details.
        """

        self.gbif_timestamp = validate_taxon_db(
            self.gbif_database, "GBIF", ["backbone"]
        )

        return self

    @model_validator(mode="after")
    def validate_projects(self) -> Self:
        """Validate and load a project database.

        This private function checks whether a project_database path: exists, is a CSV
        file, and contains project data. It populates the instance ``project_id``
        attribute.
        """

        if self.project_database is None:
            LOGGER.info("Configuration does not use project IDs.")
            return self

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
        return self


@loggerinfo_push_pop("Configuring Resources")
def load_resources(config: dict | Path | None = None) -> ResourcesPydantic:
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
    """

    cfg_path: Path | None = None
    config_found: bool = False
    standard_path: Path = Path("safedata_validator") / "safedata_validator.cfg"

    if isinstance(config, dict):
        config_found = True
    elif isinstance(config, Path):
        cfg_path = config
        config_found = True

    # Now resolve what to use in order of priority
    if not config_found:
        # Look for a path setting via an environment variable
        config_env_path = os.getenv("SAFEDATA_VALIDATOR_CONFIG")
        if config_env_path is not None:
            cfg_path = Path(config_env_path)
            config_found = True

    if not config_found:
        # Get the standard user config paths for the platform
        cfg_path = Path(appdirs.user_config_dir()) / standard_path
        if cfg_path.exists():
            config_found = True

    if not config_found:
        # Get the standard site config paths for the platform
        cfg_path = Path(appdirs.site_config_dir()) / standard_path
        if cfg_path.exists():
            config_found = True

    if not config_found:
        log_and_raise("No configuration data provided or found.", RuntimeError)

    if cfg_path is None:
        LOGGER.info("Configuring resources from dictionary")
    else:
        LOGGER.info(f"Configuring resources from file: {cfg_path}")
        with open(cfg_path, "rb") as fobj:
            config = tomllib.load(fobj)

    resources = ResourcesPydantic.model_validate(config)

    return resources


def validate_gbif_db(db_path: Path) -> str:
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

    LOGGER.info(f"Validating database: {db_path}")

    # Connect to the file (which might or might not be a database containing the
    # required tables)
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        # Check that it is a database by running a query
        try:
            db_tables = conn.execute(
                "SELECT name FROM sqlite_master WHERE type ='table';"
            )
        except sqlite3.DatabaseError:
            log_and_raise("GBIF database not an SQLite3 file", ValueError)

        # Check the required tables against found tables
        db_tables_set = {rw[0] for rw in db_tables.fetchall()}
        required_tables = set(["gbif_backbone", "timestamp"])
        missing = required_tables.difference(db_tables_set)

        if missing:
            log_and_raise(
                "GBIF database does not contain required tables: ",
                ValueError,
                extra={"join": missing},
            )

        # Check the timestamp table contains a single ISO date
        cursor = conn.execute("select * from timestamp;")
        timestamp = cursor.fetchall()

    # Is there one unique date in the table
    if len(timestamp) != 1:
        log_and_raise(
            "GBIF database timestamp table contains more than one entry.", RuntimeError
        )

    try:
        # Extract first entry in first row
        timestamp_entry = timestamp[0][0]
        isoparse(timestamp_entry)
    except ValueError:
        log_and_raise("GBIF database timestamp value is not an ISO date.", RuntimeError)

    return timestamp_entry

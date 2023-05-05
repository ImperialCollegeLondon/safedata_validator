"""General tests to check that config files are handled sensibly."""
import os
from contextlib import contextmanager
from logging import CRITICAL, ERROR, INFO, WARNING

import pytest

from safedata_validator.resources import Resources

from .conftest import FIXTURE_FILES, log_check


@contextmanager
def does_not_raise():
    yield


@pytest.fixture()
def config_file_list():

    return [
        f"gazetteer = {FIXTURE_FILES.rf.gaz_file}",
        f"location_aliases = {FIXTURE_FILES.rf.localias_file}",
        f"gbif_database = {FIXTURE_FILES.rf.gbif_file}",
        f"ncbi_database = {FIXTURE_FILES.rf.ncbi_file}",
        "use_project_ids = False",
        "[extents]",
        "temporal_soft_extent = 2002-02-01, 2030-02-01",
        "temporal_hard_extent = 2002-02-01, 2030-02-01",
        "latitudinal_hard_extent = -90, 90",
        "latitudinal_soft_extent = -4, 2",
        "longitudinal_hard_extent = -180, 180",
        "longitudinal_soft_extent = 110, 120",
        "[zenodo]",
        "community_name = safe",
        "use_sandbox = True",
        "zenodo_sandbox_api = https://sandbox.zenodo.org",
        "zenodo_sandbox_token = xyz",
        "zenodo_api = https://api.zenodo.org",
        "zenodo_token = xyz",
    ]


@pytest.fixture()
def config_file_dict():

    return {
        "gazetteer": FIXTURE_FILES.rf.gaz_file,
        "location_aliases": FIXTURE_FILES.rf.localias_file,
        "gbif_database": FIXTURE_FILES.rf.gbif_file,
        "ncbi_database": FIXTURE_FILES.rf.ncbi_file,
        "use_project_ids": True,
        "extents": {
            "temporal_soft_extent": ["2002-02-01", "2030-02-01"],
            "temporal_hard_extent": ["2002-02-01", "2030-02-01"],
            "latitudinal_hard_extent": [-90, 90],
            "latitudinal_soft_extent": [-4, 2],
            "longitudinal_hard_extent": [-180, 180],
            "longitudinal_soft_extent": [110, 120],
        },
        "zenodo": {
            "community_name": "safe",
            "use_sandbox": True,
            "zenodo_sandbox_api": "https://sandbox.zenodo.org",
            "zenodo_sandbox_token": "xyz",
            "zenodo_api": "https://api.zenodo.org",
            "zenodo_token": "xyz",
        },
    }


def nested_set(dic, keys, value):
    """Convenience function for setting dict config values."""
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value


@pytest.mark.parametrize(
    "config", ["config_file_dict", "config_file_list"], ids=("cfgdict", "cfglist")
)
@pytest.mark.parametrize(
    "dict_mod, list_mod, expected_exception, expected_log",
    [
        pytest.param(
            tuple(),
            tuple(),
            None,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (INFO, "Validating NCBI database: "),
            ),
            id="All correct",
        ),
        pytest.param(
            ((["gazetteer"], ""),),
            ((0, "gazetteer = "),),
            RuntimeError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (CRITICAL, "Gazetteer file missing in configuration"),
            ),
            id="Gazetteer not set",
        ),
        pytest.param(
            ((["gazetteer"], FIXTURE_FILES.mf),),
            ((0, f"gazetteer = {FIXTURE_FILES.mf}"),),
            OSError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (CRITICAL, "Gazetteer file not found"),
            ),
            id="Gazetteer path missing",
        ),
        pytest.param(
            ((["gazetteer"], FIXTURE_FILES.rf.gbif_file),),
            ((0, f"gazetteer = {FIXTURE_FILES.rf.gbif_file}"),),
            OSError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (CRITICAL, "Gazetteer file not valid JSON"),
            ),
            id="Gazetteer not JSON",
        ),
        pytest.param(
            ((["gazetteer"], FIXTURE_FILES.rf.json_not_locations),),
            ((0, f"gazetteer = {FIXTURE_FILES.rf.json_not_locations}"),),
            RuntimeError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (CRITICAL, "Gazetteer data not a GeoJSON Feature Collection"),
            ),
            id="Gazetteer not GeoJSON",
        ),
        pytest.param(
            ((["location_aliases"], ""),),
            ((1, "location_aliases = "),),
            RuntimeError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (CRITICAL, "Location aliases file missing in configuration"),
            ),
            id="Loc aliases not set",
        ),
        pytest.param(
            ((["location_aliases"], FIXTURE_FILES.mf),),
            ((1, f"location_aliases = {FIXTURE_FILES.mf}"),),
            OSError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (CRITICAL, "Location aliases file not found"),
            ),
            id="Loc aliases path missing",
        ),
        pytest.param(
            ((["location_aliases"], FIXTURE_FILES.rf.gbif_file),),
            ((1, f"location_aliases = {FIXTURE_FILES.rf.gbif_file}"),),
            ValueError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (
                    CRITICAL,
                    "Location aliases file not readable as a CSV file with valid "
                    "headers",
                ),
            ),
            id="Loc aliases not text",
        ),
        pytest.param(
            ((["location_aliases"], FIXTURE_FILES.rf.gaz_file),),
            ((1, f"location_aliases = {FIXTURE_FILES.rf.gaz_file}"),),
            ValueError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (
                    CRITICAL,
                    "Location aliases file not readable as a CSV file with valid "
                    "headers",
                ),
            ),
            id="Loc text not correct CSV",
        ),
        pytest.param(
            ((["location_aliases"], FIXTURE_FILES.rf.empty_localias_file),),
            ((1, f"location_aliases = {FIXTURE_FILES.rf.empty_localias_file}"),),
            ValueError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (CRITICAL, "Location aliases file is empty"),
            ),
            id="Loc alias file empty",
        ),
        pytest.param(
            ((["gbif_database"], ""),),
            ((2, "gbif_database = "),),
            ValueError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (CRITICAL, "GBIF database not set in configuration"),
            ),
            id="GBIF not set",
        ),
        pytest.param(
            ((["gbif_database"], FIXTURE_FILES.mf),),
            ((2, f"gbif_database = {FIXTURE_FILES.mf}"),),
            FileNotFoundError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (CRITICAL, "GBIF database not found"),
            ),
            id="GBIF path missing",
        ),
        pytest.param(
            ((["gbif_database"], FIXTURE_FILES.rf.gaz_file),),
            ((2, f"gbif_database = {FIXTURE_FILES.rf.gaz_file}"),),
            ValueError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (CRITICAL, "GBIF database not an SQLite3 file"),
            ),
            id="GBIF not SQLite",
        ),
        pytest.param(
            ((["gbif_database"], FIXTURE_FILES.rf.ncbi_file),),
            ((2, f"gbif_database = {FIXTURE_FILES.rf.ncbi_file}"),),
            ValueError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (CRITICAL, "GBIF database does not contain required tables"),
            ),
            id="GBIF wrong SQLite",
        ),
        pytest.param(
            ((["ncbi_database"], ""),),
            ((3, "ncbi_database = "),),
            ValueError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (INFO, "Validating NCBI database: "),
                (CRITICAL, "NCBI database not set in configuration"),
            ),
            id="NCBI not set",
        ),
        pytest.param(
            ((["ncbi_database"], FIXTURE_FILES.mf),),
            ((3, f"ncbi_database = {FIXTURE_FILES.mf}"),),
            FileNotFoundError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (INFO, "Validating NCBI database: "),
                (CRITICAL, "NCBI database not found"),
            ),
            id="NCBI path missing",
        ),
        pytest.param(
            ((["ncbi_database"], FIXTURE_FILES.rf.gaz_file),),
            ((3, f"ncbi_database = {FIXTURE_FILES.rf.gaz_file}"),),
            ValueError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (INFO, "Validating NCBI database: "),
                (CRITICAL, "NCBI database not an SQLite3 file"),
            ),
            id="NCBI not SQLite",
        ),
        pytest.param(
            ((["ncbi_database"], FIXTURE_FILES.rf.gbif_file),),
            ((3, f"ncbi_database = {FIXTURE_FILES.rf.gbif_file}"),),
            ValueError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (INFO, "Validating NCBI database: "),
                (CRITICAL, "NCBI database does not contain required tables"),
            ),
            id="NCBI wrong SQLite",
        ),
        pytest.param(
            ((["extents", "latitudinal_hard_extent"], ["-90deg", "90deg"]),),
            ((8, "latitudinal_hard_extent = -90deg, 90deg"),),
            RuntimeError,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (CRITICAL, "Configuration issues"),
                (CRITICAL, "In config 'extents.latitudinal_hard_extent':"),
            ),
            id="Testing configparser misconfig",
        ),
    ],
)
def test_load_resources_by_arg(
    config_filesystem,
    request,
    caplog,
    config,
    dict_mod,
    list_mod,
    expected_exception,
    expected_log,
):
    """Test to check that bad config inputs are logged correctly.

    This test uses the ability to load a config from a list or a dict to
    validate the behaviour of the config with bad inputs. Arguably overkill to
    run this using both list and dict inputs...
    """

    config = request.getfixturevalue(config)

    if expected_exception is None:
        ctxt_manager = does_not_raise()
    else:
        ctxt_manager = pytest.raises(expected_exception)

    with ctxt_manager:

        # Update the config, depending on whether the config is a list or a dict.
        if isinstance(config, list):
            for row, mod in list_mod:
                config[row] = mod
        elif isinstance(config, dict):
            for keys, mod in dict_mod:
                nested_set(config, keys, mod)

        Resources(config)

    log_check(caplog, expected_log)


@pytest.mark.parametrize(
    "filepath, expected_log",
    [
        (
            FIXTURE_FILES.vf.fix_config,
            (
                (INFO, "Configuring Resources"),
                (INFO, "Configuring resources from init "),
                (INFO, "Validating gazetteer: "),
                (INFO, "Validating location aliases: "),
                (INFO, "Validating GBIF database: "),
                (INFO, "Validating NCBI database: "),
            ),
        ),
    ],
)
def test_load_resources_by_file(config_filesystem, caplog, filepath, expected_log):
    """Loads a config from a file to validate the behaviour of the config."""

    Resources(filepath)

    log_check(caplog, expected_log)


# TODO  - there may be a way to combine these, but they rely on different fake
#         file systems as the first argument, so not straightforward


@pytest.mark.parametrize(
    "expected_log",
    [
        (
            (INFO, "Configuring Resources"),
            (CRITICAL, "No user config in"),
            (CRITICAL, "No site config in"),
            (CRITICAL, "No config files provided or found"),
        ),
    ],
)
def test_load_resources_from_missing_config(config_filesystem, caplog, expected_log):
    """This test checks failure modes when no config is provided."""

    with pytest.raises(RuntimeError):

        Resources()

        log_check(caplog, expected_log)


@pytest.mark.parametrize(
    "expected_log",
    [
        (
            (INFO, "Configuring Resources"),
            (INFO, "Configuring resources from user "),
            (INFO, "Validating gazetteer: "),
            (INFO, "Validating location aliases: "),
            (INFO, "Validating GBIF database: "),
            (INFO, "Validating NCBI database: "),
        ),
    ],
)
def test_load_resources_from_user_config(user_config_file, caplog, expected_log):
    """This test uses the ability to find a user config file."""

    Resources()

    log_check(caplog, expected_log)


@pytest.mark.parametrize(
    "expected_log",
    [
        (
            (INFO, "Configuring Resources"),
            (INFO, "Configuring resources from site "),
            (INFO, "Validating gazetteer: "),
            (INFO, "Validating location aliases: "),
            (INFO, "Validating GBIF database: "),
            (INFO, "Validating NCBI database: "),
        ),
    ],
)
def test_load_resources_from_site_config(site_config_file, caplog, expected_log):
    """This test uses the ability to find a site config file."""

    Resources()

    log_check(caplog, expected_log)

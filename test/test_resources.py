import pytest
import os
from contextlib import contextmanager
from logging import CRITICAL, ERROR, WARNING, INFO

from safedata_validator.resources import Resources
from .conftest import FIXTURE_FILES

@contextmanager
def does_not_raise():
    yield


@pytest.fixture()
def config_file_list():

    return [f'locations = {FIXTURE_FILES.rf.loc_file}',
            f'gbif_database = {FIXTURE_FILES.rf.gbif_file}',
            f'ncbi_database = {FIXTURE_FILES.rf.ncbi_file}',
            '[extents]',
            'temporal_soft_extent = 2002-02-01, 2030-02-01',
            'temporal_hard_extent = 2002-02-01, 2030-02-01',
            'latitudinal_hard_extent = -90, 90',
            'latitudinal_soft_extent = -4, 2',
            'longitudinal_hard_extent = -180, 180',
            'longitudinal_soft_extent = 110, 120',
            '[zenodo]',
            'community_name = safe',
            'use_sandbox = True',
            'zenodo_sandbox_api = https://sandbox.zenodo.org',
            'zenodo_sandbox_token = xyz',
            'zenodo_api = https://api.zenodo.org',
            'zenodo_token = xyz']


@pytest.fixture()
def config_file_dict():

    return {'locations': FIXTURE_FILES.rf.loc_file,
            'gbif_database': FIXTURE_FILES.rf.gbif_file,
            'ncbi_database': FIXTURE_FILES.rf.ncbi_file,
            'extents': {
               'temporal_soft_extent':  ['2002-02-01', '2030-02-01'],
               'temporal_hard_extent':  ['2002-02-01', '2030-02-01'],
               'latitudinal_hard_extent':  [-90, 90],
               'latitudinal_soft_extent':  [-4, 2],
               'longitudinal_hard_extent':  [-180, 180],
               'longitudinal_soft_extent':  [110, 120]},
            'zenodo': {
               'community_name':  'safe',
               'use_sandbox':  True,
               'zenodo_sandbox_api':  'https://sandbox.zenodo.org',
               'zenodo_sandbox_token':  'xyz',
               'zenodo_api':  'https://api.zenodo.org',
               'zenodo_token':  'xyz'}}


def nested_set(dic, keys, value):
    """Convenience function for setting dict config values"""
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value

@pytest.mark.parametrize(
    'config',
     ['config_file_dict','config_file_list'],
     ids=('cfgdict', 'cfglist'))
@pytest.mark.parametrize(
    'dict_mod, list_mod, expected_exception, expected_log',
    [ ( tuple(),
        tuple(),
        None,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (INFO, 'Validating local NCBI database: '))),
      # Locations is mandatory
      ( ( (['locations'], ''),),
        ( (0, 'locations = '),),
        RuntimeError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (CRITICAL, 'Locations file missing in configuration'))),
      # Locations is missing
      ( ( (['locations'], FIXTURE_FILES.mf ),),
        ( (0, f'locations = {FIXTURE_FILES.mf}'),),
        OSError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (CRITICAL, 'Local locations file not found'))),
      # Locations is a file but not text
      ( ( (['locations'], FIXTURE_FILES.rf.gbif_file ),),
        ( (0, f'locations = {FIXTURE_FILES.rf.gbif_file}'),),
        OSError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (CRITICAL, 'Local locations file not JSON encoded'))),
      # Locations is a json file but not locations
      ( ( (['locations'], FIXTURE_FILES.rf.json_not_locations),),
        ( (0, f'locations = {FIXTURE_FILES.rf.json_not_locations}'),),
        RuntimeError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (CRITICAL, 'Locations data malformed'))),
      # GBIF defaults to online
      ( ( (['gbif_database'], ''),),
        ( (1, 'gbif_database = '),),
        None,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Using GBIF online API to validate taxonomy'),
          (INFO, 'Validating local NCBI database: '))),
      # NCBI defaults to online
      ( ( (['ncbi_database'], ''),),
        ( (2, 'ncbi_database = '),),
        None,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (INFO, 'Using NCBI online API to validate taxonomy'))),
      # GBIF is missing
      ( ( (['gbif_database'], FIXTURE_FILES.mf),),
        ( (1, f'gbif_database = {FIXTURE_FILES.mf}'),),
        OSError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (CRITICAL, 'Local GBIF database not found'))),
      # NCBI is missing
      ( ( (['ncbi_database'], FIXTURE_FILES.mf),),
        ( (2, f'ncbi_database = {FIXTURE_FILES.mf}'),),
        OSError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (INFO, 'Validating local NCBI database: '),
          (CRITICAL, 'Local NCBI database not found'))),
      # GBIF database is not sqlite
      ( ( (['gbif_database'], FIXTURE_FILES.rf.loc_file),),
        ( (1, f'gbif_database = {FIXTURE_FILES.rf.loc_file}'),),
        OSError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (CRITICAL, 'Local SQLite database not valid'))),
      # NCBI database is not sqlite
      ( ( (['ncbi_database'], FIXTURE_FILES.rf.loc_file),),
        ( (2, f'ncbi_database = {FIXTURE_FILES.rf.loc_file}'),),
        OSError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (INFO, 'Validating local NCBI database: '),
          (CRITICAL, 'Local SQLite database not valid'))),
      # GBIF database is sqlite but not GBIF
      ( ( (['gbif_database'], FIXTURE_FILES.rf.sqlite_not_gbif),),
        ( (1, f'gbif_database = {FIXTURE_FILES.rf.sqlite_not_gbif}'),),
        RuntimeError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (CRITICAL, 'Local GBIF database does not contain the backbone table'))),
      # GBIF database is sqlite but not NCBI
      ( ( (['ncbi_database'], FIXTURE_FILES.rf.sqlite_not_gbif),),
        ( (2, f'ncbi_database = {FIXTURE_FILES.rf.sqlite_not_gbif}'),),
        RuntimeError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (INFO, 'Validating local NCBI database: '),
          (CRITICAL, 'Local NCBI database is missing either the nodes, names or merge table'))),
      # Test another misconfig
      ( ( (['extents', 'latitudinal_hard_extent'], ['-90deg', '90deg']),),
        ( (6, 'latitudinal_hard_extent = -90deg, 90deg'),),
        RuntimeError,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (CRITICAL, 'Configuration issues'),
          (CRITICAL, "In config 'extents.latitudinal_hard_extent':"))),
    ])
def test_load_resources_by_arg(config_filesystem, request, caplog, config,
                               dict_mod, list_mod, expected_exception, expected_log):
    """This test uses the ability to load a config from a list or a dict to
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

        res = Resources(config)

        assert len(expected_log) == len(caplog.records)

        assert all([exp[0] == rec.levelno
                    for exp, rec in zip(expected_log, caplog.records)])
        assert all([exp[1] in rec.message
                    for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize(
    'filepath, expected_log',
    [ ( FIXTURE_FILES.vf.fix_cfg_local,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (INFO, 'Using NCBI online API to validate taxonomy'))),
      # GBIF defaults to online
      ( FIXTURE_FILES.vf.fix_cfg_remote,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Using GBIF online API to validate taxonomy'),
          (INFO, 'Using NCBI online API to validate taxonomy'))),
      # NCBI online unless set otherwise
      ( FIXTURE_FILES.vf.fix_cfg_local_ncbi,
        ( (INFO, 'Configuring Resources'),
          (INFO, 'Configuring resources from init '),
          (INFO, 'Validating locations: '),
          (INFO, 'Validating local GBIF database: '),
          (INFO, 'Validating local NCBI database: '))),
    ])
def test_load_resources_by_file(config_filesystem, caplog, filepath, expected_log):
    """This test uses the ability to load a config from a file to
    validate the behaviour of the config.
    """

    res = Resources(filepath)

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


# TODO  - there may be a way to combine these, but they rely on different fake
#         file systems as the first argument, so not straightforward


@pytest.mark.parametrize(
    'expected_log',
    [ ( (INFO, 'Configuring Resources'),
        (CRITICAL, "No user config in"),
        (CRITICAL, "No site config in"),
        (CRITICAL, "No config files provided or found")),
    ])
def test_load_resources_from_missing_config(config_filesystem, caplog, expected_log):
    """This test checks failure modes when no config is provided
    """

    with pytest.raises(RuntimeError):

      res = Resources()

      assert len(expected_log) == len(caplog.records)

      assert all([exp[0] == rec.levelno
                  for exp, rec in zip(expected_log, caplog.records)])
      assert all([exp[1] in rec.message
                  for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize(
    'expected_log',
    [ ( (INFO, 'Configuring Resources'),
        (INFO, 'Configuring resources from user '),
        (INFO, 'Validating locations: '),
        (INFO, 'Validating local GBIF database: '),
        (INFO, 'Validating local NCBI database: ')),
    ])
def test_load_resources_from_user_config(user_config_file, caplog, expected_log):
    """This test uses the ability to find a user config file
    """

    res = Resources()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize(
    'expected_log',
    [ ( (INFO, 'Configuring Resources'),
        (INFO, 'Configuring resources from site '),
        (INFO, 'Validating locations: '),
        (INFO, 'Validating local GBIF database: '),
        (INFO, 'Validating local NCBI database: ')),
    ])
def test_load_resources_from_site_config(site_config_file, caplog, expected_log):
    """This test uses the ability to find a site config file
    """

    res = Resources()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

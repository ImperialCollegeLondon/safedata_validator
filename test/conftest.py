import os
from collections import OrderedDict

import simplejson
import openpyxl
import pytest

from safedata_validator.taxa import Taxa
from safedata_validator.locations import Locations
from safedata_validator.field import Dataset
from safedata_validator.resources import Resources
from safedata_validator.field import DataWorksheet

"""
This file contains fixtures that will be available to all test suites.
"""

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')

bad_file_path = os.path.join(FIXTURE_DIR, 'Test_format_bad.xlsx')
good_file_path = os.path.join(FIXTURE_DIR, 'Test_format_good.xlsx')

@pytest.fixture()
def config_filesystem(fs):
    """Create a config and resources setup for testing

    Testing requires access to the configuration files for the package resources
    and these are not going to be set up in default locations for testing. This
    fixture uses the pyfakefs plugin for pytest to create a fake file system
    containing config files and versions of the resources, cut down to save
    filespace.

    See also: https://gist.github.com/peterhurford/09f7dcda0ab04b95c026c60fa49c2a68

    Args:
        fs: The pyfakefs plugin object

    Returns:
        A fake filesystem containing a config file in the root, pointing to the
        actual fixture files in their correct location.
    """


    # Get paths to existing resource files
    loc_file = os.path.join(FIXTURE_DIR, 'locations.json')
    gbif_file = os.path.join(FIXTURE_DIR, 'gbif_backbone_truncated.sqlite')

    # Create some config files with correct paths - this doesn't feel right - could
    # be faking them but I can't work out how to do this correctly in fixtures

    local_config = os.path.join(FIXTURE_DIR, 'safedata_validator_local.json')
    remote_config = os.path.join(FIXTURE_DIR, 'safedata_validator_remote.json')

    # Point to real locations of test fixture example resources
    fs.add_real_file(loc_file)
    fs.add_real_file(gbif_file)

    # The config safedata_validator.json _cannot_ be a real fixture file
    # because it contains the test machine specific paths to the resources,
    # so fake it and store it in the fake fixture directory
    cfg_contents = {'locations': loc_file, 'gbif_database': gbif_file}
    cfg_contents = simplejson.dumps(cfg_contents)

    fs.create_file(local_config, contents=cfg_contents)

    cfg_contents = {'locations': loc_file, 'gbif_database': None}
    cfg_contents = simplejson.dumps(cfg_contents)

    fs.create_file(remote_config, contents=cfg_contents)

    yield fs


# ------------------------------------------
# Fixtures: a parameterized fixture containing examples of both local and
#     remote GBIF validators.
## ------------------------------------------

@pytest.fixture()
def resources_with_local_gbif_virtual(config_filesystem):
    """ Creates a Resource object configured to use a local GBIF database

    Returns:
        A safedata_validator.resources.Resources instance
    """

    # Get paths to existing resource files
    local_config = os.path.join(FIXTURE_DIR, 'fixtures', 
                               'safedata_validator_local.json')

    return Resources(cfg_file=local_config)



@pytest.fixture(scope='module')
def resources_with_local_gbif():
    """ Creates a Resource object configured to use a local GBIF database

    Returns:
        A safedata_validator.resources.Resources instance
    """
    
    local_config = os.path.join(FIXTURE_DIR, 'safedata_validator_local.json')
    return Resources(cfg_file=local_config)


@pytest.fixture(scope='module')
def resources_with_remote_gbif():
    """ Creates a Resource object configured to use the remote GBIF API

    Returns:
        A safedata_validator.resources.Resources instance
    """

    remote_config = os.path.join(FIXTURE_DIR, 'safedata_validator_remote.json')
    return Resources(cfg_file=remote_config)


@pytest.fixture(scope='module', params=['remote', 'local'])
def resources_local_and_remote(request, resources_with_local_gbif, resources_with_remote_gbif):
    """Parameterised fixture to run tests using both the local and remote GBIF.
    """

    if request.param == 'remote':
        return resources_with_remote_gbif
    elif request.param == 'local':
        return resources_with_local_gbif


@pytest.fixture(scope='module')
def example_excel_files(request):
    """This uses indirect parameterisation, to allow the shared fixture
    to be paired with request specific expectations rather than all pair
    combinations:

    https://stackoverflow.com/questions/70379640
    """
    if request.param == 'good':
        path = os.path.join(FIXTURE_DIR, 'Test_format_good.xlsx')
        wb = openpyxl.load_workbook(path, read_only=True)
        return wb
    elif request.param == 'bad':
        path = os.path.join(FIXTURE_DIR, 'Test_format_bad.xlsx')
        wb = openpyxl.load_workbook(path, read_only=True)
        return wb


# Fixtures to provide Taxon, Locations, Dataset, Dataworksheet 
# and field meta objects for testing

@pytest.fixture(scope='module')
def fixture_taxa(resources_with_local_gbif):
    """Fixture to provide a taxon object with a couple of names. These examples
    need to be in the cutdown local GBIF testing database in fixtures.
    """

    taxa = Taxa(resources_with_local_gbif)

    test_taxa = [
        ('C_born', 
            ['Crematogaster borneensis', 'Species', None, None], 
            None), 
        ('V_salv', 
            ['Varanus salvator', 'Species', None, None], 
            None),]
    
    for tx in test_taxa:
        taxa.validate_and_add_taxon(tx)
    
    return taxa


@pytest.fixture(scope='module')
def fixture_locations(resources_with_local_gbif):
    """Fixture to provide a taxon object with a couple of names. These examples
    need to be in the cutdown local GBIF testing database in fixtures.
    """

    locations = Locations(resources_with_local_gbif)

    test_locs = ['A_1', 'A_2', 1, 2]
    
    locations.add_known_locations(test_locs)
    
    return locations

@pytest.fixture()
def fixture_dataset(resources_with_local_gbif):
    """Fixture to provide a dataset that has been prepopulated with some taxon
    and location names for field tests.
    """

    dataset = Dataset(resources_with_local_gbif)
    
    test_taxa = [
        ('C_born', 
            ['Crematogaster borneensis', 'Species', None, None], 
            None), 
        ('V_salv', 
            ['Varanus salvator', 'Species', None, None], 
            None),]
    
    for tx in test_taxa:
        dataset.taxa.validate_and_add_taxon(tx)

    test_locs = ['A_1', 'A_2', 1, 2]
    
    dataset.locations.add_known_locations(test_locs)

    return dataset


@pytest.fixture(scope='module')
def fixture_field_meta():
    """field_meta object for use across tests
    """
    
    return OrderedDict(field_type = ['numeric', 'numeric', 'numeric'],
                       description = ['a', 'b', 'c'],
                       units = ['a', 'b', 'c'],
                       method = ['a', 'b', 'c'],
                       field_name = ['a', 'b', 'c'])


@pytest.fixture()
def fixture_dataworksheet(fixture_dataset):
    """field_meta object for use across tests
    """
    
    dws = DataWorksheet({'name': 'DF',
                         'title': 'My data table',
                         'description': 'This is a test data worksheet',
                         'external': None},
                         dataset=fixture_dataset)

    return dws

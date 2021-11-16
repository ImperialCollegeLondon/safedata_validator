import pytest
import os
import simplejson
from safedata_validator.resources import Resources

"""
This file contains fixtures that will be available to all test suites.
"""

# Get paths to resource files
fixture_dir = os.path.join(os.path.dirname(__file__), 'fixtures')
loc_file = os.path.join(fixture_dir, 'locations.json')
gbif_file = os.path.join(fixture_dir, 'gbif_backbone_truncated.sqlite')


@pytest.fixture(scope='module')
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

    # Point to real locations of test fixture example resources
    fs.add_real_file(loc_file)
    fs.add_real_file(gbif_file)

    # The config safedata_validator.json _cannot_ be a real fixture file
    # because it contains the test machine variable paths to the resources,
    # so fake it and store it in the fake fixture directory
    cfg_contents = {'locations': loc_file, 'gbif_database': gbif_file}
    cfg_contents = simplejson.dumps(cfg_contents)

    fs.create_file(os.path.join(fixture_dir, 'safedata_validator.json'),
                   contents=cfg_contents)

    yield fs


# ------------------------------------------
# Fixtures: a parameterized fixture containing examples of both local and
#     remote GBIF validators.
#
# TODO - compare local and remote outputs:
#  https://stackoverflow.com/questions/56558823/pytest-run-all-tests-twice-and-compare-results-bet-mock-and-real
# ------------------------------------------


@pytest.fixture(scope='module')
def resources_with_local_gbif():
    """ Creates a Resource object configured to use a local GBIF database

    Returns:
        A safedata_validator.resources.Resources instance
    """

    return Resources(locations=loc_file, gbif_database=gbif_file)


@pytest.fixture(scope='module')
def resources_with_remote_gbif():
    """ Creates a Resource object configured to use the remote GBIF API

    Returns:
        A safedata_validator.resources.Resources instance
    """

    return Resources(locations=loc_file, gbif_database=None)


@pytest.fixture(scope='module', params=['remote', 'local'])
def resources_local_and_remote(request):
    """Parameterised fixture to run tests using both the local and remote GBIF.
    """

    if request.param == 'remote':
        return resources_with_remote_gbif
    elif request.param == 'local':
        return resources_with_local_gbif

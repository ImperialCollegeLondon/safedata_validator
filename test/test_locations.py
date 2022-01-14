import pytest
from logging import ERROR, WARNING, INFO

from safedata_validator.locations import Locations

@pytest.fixture()
def locations_inst(resources_with_local_gbif):
    """Fixture to provide a taxon object with a couple of names. These examples
    need to be in the cutdown local GBIF testing database in fixtures.
    """

    return Locations(resources_with_local_gbif)

@pytest.mark.parametrize(
    'known_loc_names,expected_log',
      [({'A_1','A_2','A_3'},
        ((INFO, 'Checking known locations'),)), 
       ({'A_1', 2, 'A_3'},
        ((INFO, 'Checking known locations'),
         (WARNING, 'Locations aliases used'),)), 
          ])
def test_validate_known_locations(caplog, locations_inst, 
                                  known_loc_names, expected_log):
    # These tests check that bad inputs are caught correctly - the next test
    # handles anything involving checking against the GBIF database, so this
    # checks only using the local GBIF option for speed.

    locations_inst.validate_known_locations(known_loc_names)
 
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])
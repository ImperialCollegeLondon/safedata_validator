import pytest
from logging import ERROR, WARNING, INFO
from safedata_validator.locations import Locations


@pytest.fixture()
def locations_inst(resources_with_local_gbif):
    """Fixture to provide a set of locations and aliases.
    """
    return Locations(resources_with_local_gbif)


@pytest.mark.parametrize(
    'known_loc_names,expected_log',
      [(['A_1','A_2','A_3'],
        ((INFO, 'Checking known locations'),)), 
       (['A_1', 2, 'A_3'],
        ((INFO, 'Checking known locations'),
         (WARNING, 'Locations aliases used'),)), 
       (['A_1', None, 'A_3'],
        ((INFO, 'Checking known locations'),
         (ERROR, 'Location names contains empty cells or whitespace text'),)), 
       (['A_1', '    ', 'A_3'],
        ((INFO, 'Checking known locations'),
         (ERROR, 'Location names contains empty cells or whitespace text'),)), 
       (['A_1', 'A_2', 'Not a location'],
        ((INFO, 'Checking known locations'),
         (ERROR, 'Unknown locations found'),)), 
       (['A_1', 'A_2', 3.14159],
        ((INFO, 'Checking known locations'),
         (ERROR, 'Location names contains values that are not strings or integers'),)), 
       (['A_1', 'A_2', 'A_2'],
        ((INFO, 'Checking known locations'),
         (ERROR, 'Added names contain duplicated values'),)), 
       (['A_1', 'A_2', 'A_6'],
        ((INFO, 'Checking known locations'),
         (ERROR, 'Location names already added to Location instance'),)), 
          ])
def test_add_known_locations(caplog, locations_inst, 
                                  known_loc_names, expected_log):

    # Preload site names for testing duplicate capture and remove
    # repeated "Checking known locations" across all test outputs
    locations_inst.add_known_locations(['A_6'])
    caplog.clear()

    # Test the addition of the parameterised values
    locations_inst.add_known_locations(known_loc_names)
 
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])
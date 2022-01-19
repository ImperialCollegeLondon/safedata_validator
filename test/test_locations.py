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


@pytest.mark.parametrize(
    'new_loc_dicts,expected_log',
      [# Valid inputs -  WKT, LL and WKT + LL
       ([{'location name': 'my new loc',  'type': 'POINT', 'wkt': 'POINT(117 4)'}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating WKT data'))), 
       ([{'location name': 'my new loc', 'type': 'POINT', 'latitude': 4,  'longitude': 117}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating lat / long data'))),
       ([{'location name': 'my new loc', 'type': 'POINT', 'latitude': 4, 'longitude': 117, 'wkt': 'POINT(117 4)'}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating lat / long data'),
          (INFO, 'Validating WKT data'))),
         # Inconsistent keys
       ([{'location name': 'my new loc2',  'type': 'POINT', 'wkt': 'POINT(117 4)'},
         {'location name': 'my new loc', 'type': 'POINT', 'latitude': 4,  'longitude': 117}],
        ((INFO, 'Checking new locations'),
         (ERROR, 'Inconsistent keys in add_new_locations'))), 
         # Bad names
        ([{'type': 'POINT', 'wkt': 'POINT(117 4)'}],
          ((INFO, 'Checking new locations'),
         (ERROR, 'No location name entries in add_new_locations'))), 
        ([{'location name': None,  'type': 'POINT', 'wkt': 'POINT(117 4)'}],
        ((INFO, 'Checking new locations'),
         (ERROR, 'Location names contains empty cells or whitespace text'),
         (INFO, 'Validating WKT data'))), 
        ([{'location name': '    \t',  'type': 'POINT', 'wkt': 'POINT(117 4)'}],
        ((INFO, 'Checking new locations'),
         (ERROR, 'Location names contains empty cells or whitespace text'),
         (INFO, 'Validating WKT data'))),
        ([{'location name': 123,  'type': 'POINT', 'wkt': 'POINT(117 4)'}],
        ((INFO, 'Checking new locations'),
         (ERROR, 'New location names include non-string values'),
         (INFO, 'Validating WKT data'))),
        ([{'location name': 'my new loc', 'type': 'POINT', 'wkt': 'POINT(117 4)'},
          {'location name': 'my new loc', 'type': 'POINT', 'wkt': 'POINT(117 4)'}],
        ((INFO, 'Checking new locations'),
         (ERROR, 'New location names contain duplicated values'),
         (INFO, 'Validating WKT data'))), 
        ([{'location name': 'A_1', 'type': 'POINT', 'wkt': 'POINT(117 4)'}],
        ((INFO, 'Checking new locations'),
         (ERROR, 'New location names duplicate known names and aliases'),
         (INFO, 'Validating WKT data'))), 
        ([{'location name': 'already added', 'type': 'POINT', 'wkt': 'POINT(117 4)'}],
        ((INFO, 'Checking new locations'),
         (INFO, 'Validating WKT data'),
         (ERROR, 'Location names already added to Location instance')
         )), 
         # Bad location types
        ([{'location name': 'my new loc', 'wkt': 'POINT(117 4)'}],
         ((INFO, 'Checking new locations'),
          (ERROR, 'New locations do not provide the location type'),
          (INFO, 'Validating WKT data'))),
        ([{'location name': 'my new loc',  'type': None, 'wkt': 'POINT(117 4)'}],
         ((INFO, 'Checking new locations'),
          (ERROR,'Types for new locations contains blank or whitespace entries'),
          (INFO, 'Validating WKT data'))), 
        ([{'location name': 'my new loc',  'type': '    \n', 'wkt': 'POINT(117 4)'}],
         ((INFO, 'Checking new locations'),
          (ERROR,'Types for new locations contains blank or whitespace entries'),
          (INFO, 'Validating WKT data'))),
        ([{'location name': 'my new loc',  'type': 'POINTY', 'wkt': 'POINT(117 4)'}],
         ((INFO, 'Checking new locations'),
          (ERROR,'New locations include unknown location types'),
          (INFO, 'Validating WKT data'))),
       # Bad LL data
       ([{'location name': 'my new loc', 'type': 'POINT', 'latitude': 4, 'wkt': 'POINT(117 4)'}],
         ((INFO, 'Checking new locations'),
          (ERROR, 'New locations should either latitude _and_ longitude or neither'),
          (INFO, 'Validating WKT data'))),
       ([{'location name': 'my new loc',  'type': 'POINT'}],
         ((INFO, 'Checking new locations'),
          (ERROR, 'New locations reported: you must provide Lat/Long or WKT'))), 
        ([{'location name': 'my new loc', 'type': 'POINT', 'latitude': None,  'longitude': 117}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating lat / long data'),
          (ERROR, 'Blank latitude values for new locations'))),
        ([{'location name': 'my new loc', 'type': 'POINT', 'latitude': '   \t',  'longitude': 117}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating lat / long data'),
          (ERROR, 'Blank latitude values for new locations'))),
        ([{'location name': 'my new loc', 'type': 'POINT', 'latitude': '4',  'longitude': 117}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating lat / long data'),
          (ERROR, 'Non-numeric latitude values for new locations'))),
        # Bad WKT data
        ([{'location name': 'my new loc',  'type': 'POINT', 'wkt': None}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating WKT data'),
          (ERROR, 'Blank WKT values for new locations: use NA'))), 
        ([{'location name': 'my new loc',  'type': 'POINT', 'wkt': '\t   '}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating WKT data'),
          (ERROR, 'Blank WKT values for new locations: use NA'))), 
        ([{'location name': 'my new loc',  'type': 'POINT', 'wkt': 9}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating WKT data'),
          (ERROR, 'WKT values for new location not a string'))), 
        ([{'location name': 'my new loc',  'type': 'POINT', 'wkt': 'POINT(2)'}],
         ((INFO, 'Checking new locations'),
          (INFO, 'Validating WKT data'),
          (ERROR, "ParseException: Expected number but encountered ')'"),  # From shapely.geos
          (ERROR, 'WKT information badly formatted, not geometrically valid or 3D'))),
        ])
def test_add_new_locations(caplog, locations_inst, 
                           new_loc_dicts, expected_log):

    # Preload site names for testing duplicate capture and remove
    # repeated "Checking known locations" across all test outputs
    locations_inst.add_new_locations([{'location name': 'already added', 
                                       'type': 'POINT', 
                                       'wkt': 'POINT(117 4)'}])
    caplog.clear()

    # Test the addition of the parameterised values
    locations_inst.add_new_locations(new_loc_dicts)
 
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])
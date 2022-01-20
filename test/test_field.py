import pytest
from logging import ERROR, WARNING, INFO
import datetime

from safedata_validator.taxa import Taxa
from safedata_validator.locations import Locations
from safedata_validator.dataset import Dataset
from safedata_validator.field import (BaseField, CategoricalField, GeoField, NumericField,
                                      TaxaField, LocationsField)

# Fixtures to provide Taxon, Locations, Dataset and Dataworksheet 
# instances for testing

@pytest.fixture()
def field_test_taxa(resources_with_local_gbif):
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


@pytest.fixture()
def field_test_locations(resources_with_local_gbif):
    """Fixture to provide a taxon object with a couple of names. These examples
    need to be in the cutdown local GBIF testing database in fixtures.
    """

    locations = Locations(resources_with_local_gbif)

    test_locs = ['A_1', 'A_2', 1, 2]
    
    locations.add_known_locations(test_locs)
    
    return locations

@pytest.fixture()
def field_test_dataset(resources_with_local_gbif):
    """Fixture to provide a taxon object with a couple of names. These examples
    need to be in the cutdown local GBIF testing database in fixtures.
    """

    dataset = Dataset(resources_with_local_gbif)
    
    return dataset


# Checking the helper methods

@pytest.mark.parametrize(
    'data, exp_log, exp_levels, exp_desc',
    [
     ('A simple string', 
      ((INFO, 'Checking Column testy'),),
      ['A simple string'],
      [None]
      ),
     ('A simple string;Another string', 
      ((INFO, 'Checking Column testy'),),
      ['A simple string', 'Another string'],
      [None, None]
      ),
     ('A simple string:description;Another string:description', 
      ((INFO, 'Checking Column testy'),),
      ['A simple string', 'Another string'],
      ['description', 'description']
      ),
     ('level1:A description: with a colon in it', 
      ((INFO, 'Checking Column testy'),
       (ERROR, 'Extra colons in level description.')),
      ['level1'],
      ['A description']
       ),
     ('level1;level2:unbalanced descriptions', 
      ((INFO, 'Checking Column testy'),
       (ERROR, 'Provide descriptions for either all or none of the categories')),
      ['level1', 'level2'],
      [None, 'unbalanced descriptions']),
     ('level1;level1', 
      ((INFO, 'Checking Column testy'),
       (ERROR, 'Repeated level labels')),
      ['level1', 'level1'],
      [None, None]),
     ('1;2', 
      ((INFO, 'Checking Column testy'),
       (ERROR, 'Numeric level names not permitted')),
      ['1', '2'],
      [None, None])
    ])
def test_parse_levels(caplog, data, exp_log, exp_levels, exp_desc):
    """Testing behaviour of the NumericField class in using _validate_data
    """

    fld = BaseField({'field_name': 'testy',
                     'description': 'a test',
                      'field_type': 'irrelevant for this test'})
    
    obs_lev, obs_desc = fld._parse_levels(data)
    fld.report()

    assert len(exp_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(exp_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(exp_log, caplog.records)])

    assert len(obs_lev) == len(exp_levels)
    assert all([lo == le for lo, le in zip(obs_lev, exp_levels)])

    assert len(obs_desc) == len(exp_desc)
    assert all([do == de for do, de in zip(obs_desc, exp_desc)])

# BaseField behaviour

@pytest.mark.parametrize(
    'field_meta, expected_log',
    [
     ({'field_type': 'location',
       'description': 'SAFE 2nd order point number',
       'field_name': 'Location',
       'col_idx': 1}, 
      ((INFO, "Checking Column Location"),)),
     ({'description': 'SAFE 2nd order point number',
       'field_name': 'Location',
       'col_idx': 1}, 
      ((INFO, "Checking Column Location"),
       (ERROR, "field_type descriptor missing"),)),
     ({'field_type': None,
       'description': 'SAFE 2nd order point number',
       'field_name': 'Location',
       'col_idx': 1}, 
      ((INFO, "Checking Column Location"),
       (ERROR, "field_type descriptor is blank"),)),
     ({'field_type': '      ',
       'description': 'SAFE 2nd order point number',
       'field_name': 'Location',
       'col_idx': 1}, 
      ((INFO, "Checking Column Location"),
       (ERROR, "field_type descriptor is blank"),)),
     ({'field_type': 123,
       'description': 'SAFE 2nd order point number',
       'field_name': 'Location',
       'col_idx': 1}, 
      ((INFO, "Checking Column Location"),
       (ERROR, "field_type descriptor is not a string"),)),
     ({'field_type': '  padded   ',
       'description': 'SAFE 2nd order point number',
       'field_name': 'Location',
       'col_idx': 1}, 
      ((INFO, "Checking Column Location"),
       (ERROR, "field_type descriptor has whitespace padding"),)),
     ({'field_type': 'location',
       'description': 'SAFE 2nd order point number',
       'field_name': '123_not a valid name',
       'col_idx': 1}, 
      ((INFO, "Checking Column 123_not a valid name"),
       (ERROR, "Field name is not valid"),)),
     ({'field_type': 'location',
       'description': 'SAFE 2nd order point number',
       'field_name': None,
       'col_idx': 12}, 
      ((INFO, "Checking Column L"),
       (ERROR, "field_name descriptor is blank"),)),
       ])
def test_BaseField_init(caplog, field_meta, expected_log):
    """Testing behaviour of the BaseField class in handling bad descriptors
    via __init__ and testing within the _check_meta() method.
    """

    fld = BaseField(field_meta, None)
    fld.report()
    
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

@pytest.mark.parametrize(
    'data, expected_log',
    [
     ([1, 2, 3, 4, 5, 6, 7, 8, 9], 
      ((INFO, "Checking Column Location"),)),
     ([1, 'NA', 3, 4, 5, 6, 'NA', 8, 9], 
      ((INFO, "Checking Column Location"),
       (WARNING, "2 / 9 values missing"),)),
     ([1, None, 3, 4, 5, 6, '   ', 8, 9], 
      ((INFO, "Checking Column Location"),
       (ERROR, "2 cells are blank or contain only whitespace text"),)),
     ([1, 2, 3, 4, '#REF!', 6, '#N/A', 8, 9], 
      ((INFO, "Checking Column Location"),
       (ERROR, "2 cells contain Excel formula errors"),)),
       ])
def test_BaseField_validate_data(caplog, data, expected_log):
    """Testing behaviour of the BaseField class in using _validate_data
    """

    fld = BaseField({'field_type': 'location',
                    'description': 'SAFE 2nd order point number',
                    'field_name': 'Location'}, None)
    
    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

# NumericField - has BaseField init, just test overloaded validate_data

@pytest.mark.parametrize(
    'data, expected_log',
    [
     ([1, 2, 3, 4, 5, 6, 7, 8, 9], 
      ((INFO, "Checking Column tree_height"),)),
     ([1, 'NA', 3, 4, 5, 6, 'NA', 8, 9], 
      ((INFO, "Checking Column tree_height"),
       (WARNING, "2 / 9 values missing"))),
     ([1, None, 3, 4, 5, 6, '   ', 8, 9], 
      ((INFO, "Checking Column tree_height"),
       (ERROR, "2 cells are blank or contain only whitespace text"))),
     ([1, 2, 3, 4, '#REF!', 6, '#N/A', 8, 9], 
      ((INFO, "Checking Column tree_height"),
       (ERROR, "2 cells contain Excel formula errors"))),
     ([1, 2, 3, 4, 'wrong_type', 6, 7, 8, 9], 
      ((INFO, "Checking Column tree_height"),
       (ERROR, "Cells contain non-numeric values"))),
     ([1, 2, 'NA', 4, 'wrong_type', 6, None, 8, 9], 
      ((INFO, "Checking Column tree_height"),
       (ERROR, "Cells contain non-numeric values"), 
       (WARNING, "1 / 9 values missing"),
       (ERROR, "1 cells are blank or contain only whitespace text"))),
    ])
def test_NumericField_validate_data(caplog, data, expected_log):
    """Testing behaviour of the NumericField class in using _validate_data
    """

    fld = NumericField({'field_type': 'numeric',
                        'description': 'Tree height',
                        'field_name': 'tree_height',
                        'method': 'looking',
                        'units': 'metres'}, None)
    
    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

# CategoricalField - check init and validate_data (overloaded report just emits messages)

@pytest.mark.parametrize(
    'field_meta, expected_log',
    [
     ({'field_type': 'categorical',
       'description': 'a factor',
       'field_name': 'factor1',
       'col_idx': 1}, 
      ((INFO, "Checking Column factor1"),
       (ERROR, "levels descriptor missing"),)),
     ({'field_type': 'categorical',
       'description': 'a factor',
       'field_name': 'factor1',
       'levels': None,
       'col_idx': 1}, 
      ((INFO, "Checking Column factor1"),
       (ERROR, "levels descriptor is blank"),)),
     ({'field_type': 'categorical',
       'description': 'a factor',
       'field_name': 'factor1',
       'levels': '   ',
       'col_idx': 1}, 
      ((INFO, "Checking Column factor1"),
       (ERROR, "levels descriptor is blank"),)),
     ({'field_type': 'categorical',
       'description': 'a factor',
       'field_name': 'factor1',
       'levels': 123,
       'col_idx': 1}, 
      ((INFO, "Checking Column factor1"),
       (ERROR, "levels descriptor is not a string"),)),
     ({'field_type': 'categorical',
       'description': 'a factor',
       'field_name': 'factor1',
       'levels': ' level1;level2; ',
       'col_idx': 1}, 
      ((INFO, "Checking Column factor1"),
       (ERROR, "levels descriptor has whitespace padding"),
       (ERROR, "Categories found in levels descriptor not used in data:"),)),
       ])
def test_CategoricalField_init(caplog, field_meta, expected_log):
    """Testing behaviour of the BaseField class in handling bad descriptors
    via __init__ and testing within the _check_meta() method.
    """

    fld = CategoricalField(field_meta, None)
    fld.report()
    
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

@pytest.mark.parametrize(
    'data, expected_log',
    [
     (['level1', 'level2', 'level1', 'level2', 'level1', 'level2'], 
      ((INFO, "Checking Column factor1"),)),
     (['level1', 'NA', 'level1', 'level2', 'NA', 'level2'], 
      ((INFO, "Checking Column factor1"),
       (WARNING, "2 / 6 values missing"))),
     (['level1', '    ', 'level1', 'level2', None, 'level2'], 
      ((INFO, "Checking Column factor1"),
       (ERROR, "2 cells are blank or contain only whitespace text"))),
     (['level1', '#REF!', 'level1', '#N/A', 'level1', 'level2'], 
      ((INFO, "Checking Column factor1"),
       (ERROR, "2 cells contain Excel formula errors"))),
     (['level1', 'level2', 'level1', 1234, 'level1', 'level2'], 
      ((INFO, "Checking Column factor1"),
       (ERROR, "Cells contain non-text values"))),
     (['level1', 'level2', 'NA', None, 1234, 'level2'], 
      ((INFO, "Checking Column factor1"),
       (ERROR, "Cells contain non-text values"), 
       (WARNING, "1 / 6 values missing"),
       (ERROR, "1 cells are blank or contain only whitespace text"))),
     (['level1', 'level2', 'level3', 'level2', 'level1', 'level2'], 
      ((INFO, "Checking Column factor1"),
       (ERROR, "Categories found in data missing from levels descriptor"))),
     (['level1', 'level1', 'level1', 'level1', 'level1', 'level1'], 
      ((INFO, "Checking Column factor1"),
       (ERROR, "Categories found in levels descriptor not used in data"))),
    ])
def test_CategoricalField_validate_data(caplog, data, expected_log):
    """Testing behaviour of the CategoricalField class in using validate_data
    """

    fld = CategoricalField({'field_type': 'categorical',
                            'description': 'a factor',
                            'field_name': 'factor1',
                            'levels': 'level1;level2;',
                            'col_idx': 1}, None)

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

# Taxon field - check init and validate_data (overloaded report just emits messages)

@pytest.mark.parametrize(
    'provide_taxa_instance, expected_log',
    [
      ( True,
        ( (INFO, "Checking Column taxa_field"),
          (ERROR, 'No taxa loaded'))),
      ( False,
        ( (INFO, "Checking Column taxa_field"),
          (ERROR, 'No taxon details provided for dataset'),
          (ERROR, 'No taxa loaded'))),
       ])
def test_TaxaField_init(caplog, provide_taxa_instance, expected_log):
    """Testing behaviour of the TaxaField class in handling missing taxa.
    """

    if provide_taxa_instance:
        taxa = field_test_taxa
    else:
        taxa = None

    fld = TaxaField({'field_type': 'taxa',
                     'description': 'My taxa',
                     'field_name': 'taxa_field'}, taxa=taxa)
    fld.report()
    
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize(
    'data, expected_log',
    [
     (['C_born', 'V_salv', 'C_born', 'V_salv', 'C_born', 'V_salv'], 
      ((INFO, "Checking Column taxa_field"),)),
     (['C_born', 123, 'C_born', 'V_salv', 'C_born', 'V_salv'], 
      ((INFO, "Checking Column taxa_field"),
       (ERROR, 'Cells contain non-string values'))),
     ([], 
      ((INFO, "Checking Column taxa_field"),
       (ERROR, 'No taxa loaded'))),
     (['C_born', 'V_salv', 'C_born', 'V_salv', 'C_born', 'V_salv', 'P_leo'], 
      ((INFO, "Checking Column taxa_field"),
       (ERROR, 'Includes unreported taxa'))),
    ])
def test_TaxaField_validate_data(caplog, field_test_taxa, data, expected_log):
    """Testing behaviour of the TaxaField class in using validate_data
    """

    fld = TaxaField({'field_type': 'taxa',
                     'description': 'My taxa',
                     'field_name': 'taxa_field'}, taxa=field_test_taxa)
    
    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

# Locations field - check init and validate_data (overloaded report just emits messages)

@pytest.mark.parametrize(
    'provide_loc_instance, expected_log',
    [
      ( True,
        ( (INFO, "Checking Column locations"),
          (ERROR, 'No locations loaded'))),
      ( False,
        ( (INFO, "Checking Column locations"),
          (ERROR, 'No location details provided for dataset'),
          (ERROR, 'No locations loaded'))),
       ])
def test_LocationsField_init(caplog, field_test_locations, provide_loc_instance, expected_log):
    """Testing behaviour of the LocationsField class in handling missing locations.
    """

    if provide_loc_instance:
        locs = field_test_locations
    else:
        locs = None

    fld = LocationsField({'field_name': 'locations',
                          'field_type': 'locations',
                          'description': 'my locations'}, locations=locs)
    fld.report()
    
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize(
    'data, expected_log',
    [ # Good data
      ( ['A_1', 'A_2', 'A_1', 'A_2', 'A_1', 'A_2'], 
        ((INFO, "Checking Column locations"),)),
      ( ['A_1', 'A_2', 1, 'A_2', 'A_1', 2], 
        ((INFO, "Checking Column locations"),)),
      # Bad data
      ( ['A_1', 'A_2', 16.2, 'A_2', 'A_1', datetime.datetime.now()], 
        ((INFO, "Checking Column locations"),
         (ERROR, "Cells contain invalid location values"))),
      ( [], 
        ((INFO, "Checking Column locations"),
         (ERROR, "No locations loaded"))),
      ( ['A_1', 'A_2', 'A_3', 'A_2', 'A_1', 'A_3'], 
        ((INFO, "Checking Column locations"),
         (ERROR, "Includes unreported locations"))),
    ])
def test_LocationsField_validate_data(caplog, field_test_locations, data, expected_log):
    """Testing behaviour of the TaxaField class in using validate_data
    """

    fld = LocationsField({'field_name': 'locations',
                          'field_type': 'locations',
                          'description': 'my locations'}, 
                          locations=field_test_locations)

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

# Geographic coordinates field.

@pytest.mark.parametrize(
    'provide_ds_instance, expected_log',
    [
      ( True,
        ( (INFO, "Checking Column geocoords"),)),
      ( False,
        ( (INFO, "Checking Column geocoords"),
          (ERROR, 'No dataset object provided - cannot update extents'))),
       ])
def test_GeoField_init(caplog, field_test_dataset, provide_ds_instance, expected_log):
    """Testing behaviour of the GeoField class in handling missing dataset.
    """

    if provide_ds_instance:
        ds = field_test_dataset
    else:
        ds = None

    fld = GeoField({'field_name': 'geocoords',
                    'field_type': 'latitude',
                    'description': 'my gcs'}, 
                     dataset=ds)
    fld.report()
    
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize(
    'data, expected_log',
    [ # Good data
      ( [1,2,3,4,5,6], 
        ( (INFO, "Checking Column geocoords"),)),
      # Bad inputs
      ( [1,'2',3,4,5,6], 
        ( (INFO, "Checking Column geocoords"),
          (ERROR, 'Field contains non-numeric data'))),
      ( [1,'2°',3,4,5,6], 
        ( (INFO, "Checking Column geocoords"),
          (ERROR, 'Field contains non-numeric data'),
          (WARNING, 'Possible degrees minutes and seconds formatting?'))),
      ( [1,2,3,4,5,600], 
        ( (INFO, "Checking Column geocoords"),
          (ERROR, 'Values range (1, 600) exceeds hard bounds'))),
      ( [1,2,3,4,5,60], 
        ( (INFO, "Checking Column geocoords"),
          (WARNING, 'Values range (1, 60) exceeds soft bounds'))),
    ])
@pytest.mark.parametrize(
    'which', ['latitude', 'longitude'])
def test_GeoField_validate_data(caplog, field_test_dataset, data, expected_log, which):
    """Testing behaviour of the TaxaField class in using validate_data
    """

    fld = GeoField({'field_name': 'geocoords',
                    'field_type': which,
                    'description': 'my gcs'}, 
                    dataset=field_test_dataset)

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


















[
 {'field_type': 'date',
  'description': 'Collection date',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'Date',
  'col_idx': 2},
 {'field_type': 'latitude',
  'description': 'GPS Lat for sample point',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'Latitude',
  'col_idx': 3},
 {'field_type': 'longitude',
  'description': 'GPS Long for sample point',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'Longitude',
  'col_idx': 4},
 {'field_type': 'replicate',
  'description': 'Four replicates',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'QuadratNumber',
  'col_idx': 5},
 {'field_type': 'numeric',
  'description': 'Something',
  'method': 'Presence/absence of something',
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': 'binary',
  'field_name': 'Status',
  'col_idx': 6},
 {'field_type': 'categorical',
  'description': 'Something else',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': 'A:Treatment A;B:Treatment B;C:Treatment C;',
  'units': None,
  'field_name': 'Treat',
  'col_idx': 7},
 {'field_type': 'categorical',
  'description': 'Something else',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': 'G:Level G;V:Level V',
  'units': None,
  'field_name': 'Level',
  'col_idx': 8},
 {'field_type': 'categorical',
  'description': 'Something else',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': 'D:Day;N:Night',
  'units': None,
  'field_name': 'Time',
  'col_idx': 9},
 {'field_type': 'numeric',
  'description': 'Something else',
  'method': 'Bucket',
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': 'mm',
  'field_name': 'Rain',
  'col_idx': 10},
 {'field_type': 'numeric',
  'description': 'Something else',
  'method': 'Count',
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': 'count',
  'field_name': 'Group',
  'col_idx': 11},
 {'field_type': 'time',
  'description': 'A time',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'Obs_Time',
  'col_idx': 12},
 {'field_type': 'datetime',
  'description': 'A date time field',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'Datetime',
  'col_idx': 13},
 {'field_type': 'taxa',
  'description': 'The taxon studied',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'Taxon',
  'col_idx': 14},
 {'field_type': 'id',
  'description': 'Pit tags for the ants',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'TagNumber',
  'col_idx': 15},
 {'field_type': 'abundance',
  'description': 'The species abundance',
  'method': '50 by 50 cm quadrat',
  'taxon_field': 'Taxon',
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'Species_count',
  'col_idx': 16},
 {'field_type': 'numeric trait',
  'description': 'The average individual mass',
  'method': 'Average mass of a sample of 10 individual',
  'taxon_field': 'Taxon',
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': 'milligrams',
  'field_name': 'Ant_mass',
  'col_idx': 17},
 {'field_type': 'ordered categorical trait',
  'description': 'The caset of the individuals weighed',
  'method': None,
  'taxon_field': 'Taxon',
  'interaction_name': None,
  'interaction_field': None,
  'levels': 'Worker;Soldier;Drone;Queen',
  'units': None,
  'field_name': 'Caste',
  'col_idx': 18},
 {'field_type': 'numeric interaction',
  'description': 'Consumption of ants by lizards',
  'method': 'Number observed eaten in half hour observation',
  'taxon_field': None,
  'interaction_name': 'Water monitor:predator',
  'interaction_field': 'Taxon:prey',
  'levels': None,
  'units': 'count',
  'field_name': 'Lizard_predation',
  'col_idx': 19},
 {'field_type': 'comments',
  'description': 'My comments',
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': 'Comments',
  'col_idx': 20},
 {'field_type': None,
  'description': None,
  'method': None,
  'taxon_field': None,
  'interaction_name': None,
  'interaction_field': None,
  'levels': None,
  'units': None,
  'field_name': None,
  'col_idx': 21}]
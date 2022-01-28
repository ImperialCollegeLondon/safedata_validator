from collections import OrderedDict
import pytest
from logging import CRITICAL, ERROR, WARNING, INFO
import datetime

from safedata_validator.field import (BaseField, CategoricalField, GeoField, 
                                      NumericField, TaxaField, LocationsField, 
                                      NumericTaxonField, CategoricalTaxonField,
                                      NumericInteractionField, CategoricalInteractionField,
                                      DataWorksheet)
from test.conftest import fixture_taxa

# Checking the helper and private methods

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


@pytest.mark.parametrize(
    'tx_meta, has_taxa_object, has_dwsh_object, expected_log',
    [
      ( dict(taxon_name='foo', taxon_field='bar'),
        False, False,
        ( (ERROR, "Taxon name and taxon field both provided, use one only"),)),
      ( dict(),
        False, False,
        ( (ERROR, "One of taxon name or taxon field must be provided"),)),
      ( dict(taxon_name='foo'),
        False, False,
        ( (ERROR, "Taxon name provided but no taxa loaded"),)),
      ( dict(taxon_name='foo'),
        True, False,
        ( (ERROR, "Taxon name not found in the Taxa worksheet"),)),
      ( dict(taxon_name='C_born'),
        True, False,
        tuple()),
      ( dict(taxon_field='not_provided'),
        False, False,
        ( (CRITICAL, "Taxon field provided but no dataworksheet provided for this field"),)),
      ( dict(taxon_field='not_provided'),
        False, True,
        ( (ERROR, "Taxon field not found in this worksheet"),)),
      ( dict(taxon_field='my_taxon_field'),
        False, True,
        tuple()),
    ] 
)
def test_check_taxon_meta(caplog, fixture_taxa,
                          tx_meta, has_taxa_object, has_dwsh_object, expected_log):
    """Testing the use of the BaseField._check_taxon_meta() method
    """
    
    # Set up what information is available for taxon field validation
    tx_obj = fixture_taxa if has_taxa_object else None
    dwsh = DataWorksheet({'name': 'DF',
                          'title': 'My data table',
                          'description': 'This is a test data worksheet'})
    dwsh.taxa_fields = ['my_taxon_field']
    dwsh_obj = dwsh if has_dwsh_object else None

    # Technically, this violates the field_name last requirement, but that is
    # enforced at the dataworksheet level, not the field level. 
    field_meta = OrderedDict(field_type = 'numeric',
                             description = 'description',
                             field_name = 'field')
    field_meta.update(tx_meta)
    
    fld = BaseField(field_meta,
                    taxa = tx_obj,
                    dwsh = dwsh_obj)

    caplog.clear()

    # Test the logging from this private method.
    fld._check_taxon_meta()
    
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

@pytest.mark.parametrize(
    'iact_meta, has_taxa_object, has_dwsh_object, expected_log',
    [
      ( dict(),
        False, False,
        ( (ERROR, "At least one of interaction name or interaction field must be provided"),)),
      ( dict(interaction_name='foo:predator'),
        False, False,
        ( (ERROR, "Interaction name provided but no taxa loaded"),)),
      ( dict(interaction_name='foo:predator'),
        True, False,
        ( (ERROR, "Unknown taxa in interaction_name descriptor"),
          (ERROR, "At least two interacting taxon labels or fields must be identified"))),
      ( dict(interaction_name='C_born:predator'),
        True, False,
        ( (ERROR, "At least two interacting taxon labels or fields must be identified"),)),
      ( dict(interaction_name='C_born:predator;V_salv:prey'),  # Biologically not very plausible
        True, False,
        tuple()),
      ( dict(interaction_field='not_provided:predator'),
        False, False,
        ( (CRITICAL, "Interaction field provided but no dataworksheet provided for this field"),)),
      ( dict(interaction_field='foo:predator'),
        False, True,
        ( (ERROR, "Unknown taxon fields in interaction_field descriptor"),
          (ERROR, "At least two interacting taxon labels or fields must be identified"))),
      ( dict(interaction_field='predator:eats things'),
        False, True,
        ( (ERROR, "At least two interacting taxon labels or fields must be identified"), )),
      ( dict(interaction_field='predator:eats things;prey:gets eaten'),
        False, True,
        tuple()),
      ( dict(interaction_field='predator:eats things', 
             interaction_name='C_born:prey'),
        True, True,
        tuple()),
      ( dict(interaction_field='predator:eats things;prey:gets eaten', 
             interaction_name='C_born:decomposer'),  # Tritrophic example
        True, True,
        tuple()),
      ( dict(interaction_name='C_born;V_salv'),  
        True, False,
        ( (WARNING, "Label descriptions for interacting taxa incomplete or missing"), )),
      ( dict(interaction_field='predator;prey'),  
        False, True,
        ( (WARNING, "Label descriptions for interacting taxa incomplete or missing"), )),
      ( dict(interaction_field='predator', 
             interaction_name='C_born'),
        True, True,
        ( (WARNING, "Label descriptions for interacting taxa incomplete or missing"), )),
    ] 
)
def test_check_interaction_meta(caplog, fixture_taxa,
                                iact_meta, has_taxa_object, has_dwsh_object, expected_log):
    """Testing the use of the BaseField._check_interaction_meta() method
    """
    
    # Set up what information is available for taxon field validation
    tx_obj = fixture_taxa if has_taxa_object else None
    dwsh = DataWorksheet({'name': 'DF',
                          'title': 'My data table',
                          'description': 'This is a test data worksheet'})
    dwsh.taxa_fields = ['predator', 'prey']
    dwsh_obj = dwsh if has_dwsh_object else None

    # Technically, this violates the field_name last requirement, but that is
    # enforced at the dataworksheet level, not the field level. 
    field_meta = OrderedDict(field_type = 'numeric',
                             description = 'description',
                             field_name = 'field')
    field_meta.update(iact_meta)
    
    fld = BaseField(field_meta,
                    taxa = tx_obj,
                    dwsh = dwsh_obj)

    caplog.clear()

    # Test the logging from this private method.
    fld._check_interaction_meta()
    
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

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

# NumericField and derived classes
# - NumericField itself has BaseField init, just test overloaded validate_data 
# - Can reuse the same data to also check taxon and interaction classes
#   inheriting from NumericField for validate_data.
# - The __init__ testing duplicates testing of private methods above but this is
#   so that I can be sure the  class inheritance works as expected


@pytest.mark.parametrize('tx_meta, has_taxa_object, has_dwsh_object, expected_log',
    [
      ( dict(taxon_name='foo', taxon_field='bar'),
        False, False,
        ((ERROR, "Taxon name and taxon field both provided, use one only"),)),
      ( dict(),
        False, False,
        ( (ERROR, "One of taxon name or taxon field must be provided"),)),
      ( dict(taxon_name='foo'),
        False, False,
        ( (ERROR, "Taxon name provided but no taxa loaded"),)),
      ( dict(taxon_name='foo'),
        True, False,
        ( (ERROR, "Taxon name not found in the Taxa worksheet"),)),
      ( dict(taxon_name='C_born'),
        True, False,
        tuple()),
      ( dict(taxon_field='not_provided'),
        False, False,
        ( (CRITICAL, "Taxon field provided but no dataworksheet provided for this field"),)),
      ( dict(taxon_field='not_provided'),
        False, True,
        ( (ERROR, "Taxon field not found in this worksheet"),)),
      ( dict(taxon_field='my_taxon_field'),
        False, True,
        tuple()),
    ]
  )
def test_NumericTaxonField_init(caplog, fixture_taxa, 
                                tx_meta, has_taxa_object, has_dwsh_object, expected_log):
    """Testing behaviour of the NumericTaxonField class in using _validate_data
    """

    # Set up what information is available for taxon field validation
    tx_obj = fixture_taxa if has_taxa_object else None
    dwsh = DataWorksheet({'name': 'DF',
                          'title': 'My data table',
                          'description': 'This is a test data worksheet'})
    dwsh.taxa_fields = ['my_taxon_field']
    dwsh_obj = dwsh if has_dwsh_object else None

    caplog.clear()

    meta = {'field_type': 'abundance',
            'description': 'Number of ants',
            'field_name': 'ant_count',
            'method': 'quadrat',
            'units': 'individuals per m2'}
    
    meta.update(tx_meta)

    fld = NumericTaxonField(meta,
                            taxa = tx_obj,
                            dwsh = dwsh_obj)

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize('iact_meta, has_taxa_object, has_dwsh_object, expected_log',
    [
      ( dict(),
        False, False,
        ( (ERROR, "At least one of interaction name or interaction field must be provided"),)),
      ( dict(interaction_name='foo:predator'),
        False, False,
        ( (ERROR, "Interaction name provided but no taxa loaded"),)),
      ( dict(interaction_name='foo:predator'),
        True, False,
        ( (ERROR, "Unknown taxa in interaction_name descriptor"),
          (ERROR, "At least two interacting taxon labels or fields must be identified"))),
      ( dict(interaction_name='C_born:predator'),
        True, False,
        ( (ERROR, "At least two interacting taxon labels or fields must be identified"),)),
      ( dict(interaction_name='C_born:predator;V_salv:prey'),  # Biologically not very plausible
        True, False,
        tuple()),
      ( dict(interaction_field='not_provided:predator'),
        False, False,
        ( (CRITICAL, "Interaction field provided but no dataworksheet provided for this field"),)),
      ( dict(interaction_field='foo:predator'),
        False, True,
        ( (ERROR, "Unknown taxon fields in interaction_field descriptor"),
          (ERROR, "At least two interacting taxon labels or fields must be identified"))),
      ( dict(interaction_field='predator:eats things'),
        False, True,
        ( (ERROR, "At least two interacting taxon labels or fields must be identified"), )),
      ( dict(interaction_field='predator:eats things;prey:gets eaten'),
        False, True,
        tuple()),
      ( dict(interaction_field='predator:eats things', 
             interaction_name='C_born:prey'),
        True, True,
        tuple()),
      ( dict(interaction_field='predator:eats things;prey:gets eaten', 
             interaction_name='C_born:decomposer'),  # Tritrophic example
        True, True,
        tuple()),
      ( dict(interaction_name='C_born;V_salv'),  
        True, False,
        ( (WARNING, "Label descriptions for interacting taxa incomplete or missing"), )),
      ( dict(interaction_field='predator;prey'),  
        False, True,
        ( (WARNING, "Label descriptions for interacting taxa incomplete or missing"), )),
      ( dict(interaction_field='predator', 
             interaction_name='C_born'),
        True, True,
        ( (WARNING, "Label descriptions for interacting taxa incomplete or missing"), )),
    ] 
)
def test_NumericInteractionField_init(caplog, fixture_taxa,
                                      iact_meta, has_taxa_object, has_dwsh_object, expected_log):
    """Testing the use of the NumericInteractionField init
    """
    
    # Set up what information is available for taxon field validation
    tx_obj = fixture_taxa if has_taxa_object else None
    dwsh = DataWorksheet({'name': 'DF',
                          'title': 'My data table',
                          'description': 'This is a test data worksheet'})
    dwsh.taxa_fields = ['predator', 'prey']
    dwsh_obj = dwsh if has_dwsh_object else None

    # Technically, this violates the field_name last requirement, but that is
    # enforced at the dataworksheet level, not the field level. 
    field_meta = OrderedDict(field_type = 'numeric',
                             description = 'description',
                             field_name = 'field')
    field_meta.update(iact_meta)

    caplog.clear()

    # Test the __init__ method
    fld = NumericInteractionField(field_meta,
                                  taxa = tx_obj,
                                  dwsh = dwsh_obj)

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize('data, expected_log',
    [
      ( [1, 2, 3, 4, 5, 6, 7, 8, 9], 
        ( (INFO, "Checking Column num_data"),)),
      ( [1, 'NA', 3, 4, 5, 6, 'NA', 8, 9], 
        ( (INFO, "Checking Column num_data"),
          (WARNING, "2 / 9 values missing"))),
      ( [1, None, 3, 4, 5, 6, '   ', 8, 9], 
        ( (INFO, "Checking Column num_data"),
          (ERROR, "2 cells are blank or contain only whitespace text"))),
     (  [1, 2, 3, 4, '#REF!', 6, '#N/A', 8, 9], 
        ( (INFO, "Checking Column num_data"),
          (ERROR, "2 cells contain Excel formula errors"))),
     (  [1, 2, 3, 4, 'wrong_type', 6, 7, 8, 9], 
        ( (INFO, "Checking Column num_data"),
          (ERROR, "Cells contain non-numeric values"))),
     (  [1, 2, 'NA', 4, 'wrong_type', 6, None, 8, 9], 
        ( (INFO, "Checking Column num_data"),
          (ERROR, "Cells contain non-numeric values"), 
          (WARNING, "1 / 9 values missing"),
          (ERROR, "1 cells are blank or contain only whitespace text"))),
    ]
  )
@pytest.mark.parametrize('test_class, field_meta', 
  [ ( NumericField, 
      {'field_type': 'numeric',
       'description': 'Tree height',
       'field_name': 'num_data',
       'method': 'looking',
       'units': 'metres'}),
    ( NumericTaxonField,
      {'field_type': 'abundance',
       'description': 'Number of ants',
       'field_name': 'num_data',
       'method': 'quadrat',
       'units': 'individuals per m2',
       'taxon_name': 'C_born'}),
    ( NumericInteractionField,
      {'field_type': 'numeric interaction',
       'description': 'Number of prey eaten',
       'field_name': 'num_data',
       'method': 'cage experiment',
       'units': 'individuals per hour',
       'interaction_name': 'C_born:prey;V_salv:predator'}),
  ]
)
def test_NumericField_validate_data(caplog, fixture_taxa, test_class, field_meta,
                                         data, expected_log):
    """Testing behaviour of the NumericField and subclasses in using _validate_data
    """
    
    # Create and instance of the required class
    fld = test_class(field_meta, taxa = fixture_taxa)
    
    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


# CategoricalField and derived classes
# - Can reuse the same data to also check taxon and interaction classes
#   inheriting from CategoricalField for validate_data.
# - The __init__ testing duplicates testing of private methods above but this is
#   so that I can be sure the  class inheritance works as expected


@pytest.mark.parametrize('field_meta, expected_log',
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


@pytest.mark.parametrize('tx_meta, has_taxa_object, has_dwsh_object, expected_log',
    [
      ( dict(taxon_name='foo', taxon_field='bar'),
        False, False,
        ((ERROR, "Taxon name and taxon field both provided, use one only"),)),
      ( dict(),
        False, False,
        ( (ERROR, "One of taxon name or taxon field must be provided"),)),
      ( dict(taxon_name='foo'),
        False, False,
        ( (ERROR, "Taxon name provided but no taxa loaded"),)),
      ( dict(taxon_name='foo'),
        True, False,
        ( (ERROR, "Taxon name not found in the Taxa worksheet"),)),
      ( dict(taxon_name='C_born'),
        True, False,
        tuple()),
      ( dict(taxon_field='not_provided'),
        False, False,
        ( (CRITICAL, "Taxon field provided but no dataworksheet provided for this field"),)),
      ( dict(taxon_field='not_provided'),
        False, True,
        ( (ERROR, "Taxon field not found in this worksheet"),)),
      ( dict(taxon_field='my_taxon_field'),
        False, True,
        tuple()),
    ]
  )
def test_CategoricalTaxonField_init(caplog, fixture_taxa,
                                    tx_meta, has_taxa_object, has_dwsh_object, expected_log):
    """Testing behaviour of the NumericTaxonField class in using _validate_data
    """

    # Set up what information is available for taxon field validation
    tx_obj = fixture_taxa if has_taxa_object else None
    dwsh = DataWorksheet({'name': 'DF',
                          'title': 'My data table',
                          'description': 'This is a test data worksheet'})
    dwsh.taxa_fields = ['my_taxon_field']
    dwsh_obj = dwsh if has_dwsh_object else None

    caplog.clear()

    meta = {'field_type': 'categorical trait',
            'description': 'ant colours',
            'field_name': 'ant_colour',
            'levels': 'level1:level2'}
    
    meta.update(tx_meta)

    fld = CategoricalTaxonField(meta,
                                taxa = tx_obj,
                                dwsh = dwsh_obj)

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize('iact_meta, has_taxa_object, has_dwsh_object, expected_log',
    [
      ( dict(),
        False, False,
        ( (ERROR, "At least one of interaction name or interaction field must be provided"),)),
      ( dict(interaction_name='foo:predator'),
        False, False,
        ( (ERROR, "Interaction name provided but no taxa loaded"),)),
      ( dict(interaction_name='foo:predator'),
        True, False,
        ( (ERROR, "Unknown taxa in interaction_name descriptor"),
          (ERROR, "At least two interacting taxon labels or fields must be identified"))),
      ( dict(interaction_name='C_born:predator'),
        True, False,
        ( (ERROR, "At least two interacting taxon labels or fields must be identified"),)),
      ( dict(interaction_name='C_born:predator;V_salv:prey'),  # Biologically not very plausible
        True, False,
        tuple()),
      ( dict(interaction_field='not_provided:predator'),
        False, False,
        ( (CRITICAL, "Interaction field provided but no dataworksheet provided for this field"),)),
      ( dict(interaction_field='foo:predator'),
        False, True,
        ( (ERROR, "Unknown taxon fields in interaction_field descriptor"),
          (ERROR, "At least two interacting taxon labels or fields must be identified"))),
      ( dict(interaction_field='predator:eats things'),
        False, True,
        ( (ERROR, "At least two interacting taxon labels or fields must be identified"), )),
      ( dict(interaction_field='predator:eats things;prey:gets eaten'),
        False, True,
        tuple()),
      ( dict(interaction_field='predator:eats things', 
             interaction_name='C_born:prey'),
        True, True,
        tuple()),
      ( dict(interaction_field='predator:eats things;prey:gets eaten', 
             interaction_name='C_born:decomposer'),  # Tritrophic example
        True, True,
        tuple()),
      ( dict(interaction_name='C_born;V_salv'),  
        True, False,
        ( (WARNING, "Label descriptions for interacting taxa incomplete or missing"), )),
      ( dict(interaction_field='predator;prey'),  
        False, True,
        ( (WARNING, "Label descriptions for interacting taxa incomplete or missing"), )),
      ( dict(interaction_field='predator', 
             interaction_name='C_born'),
        True, True,
        ( (WARNING, "Label descriptions for interacting taxa incomplete or missing"), )),
    ] 
)
def test_CategoricalInteractionField_init(caplog, fixture_taxa,
                                          iact_meta, has_taxa_object, has_dwsh_object, expected_log):
    """Testing the use of the CategoricalInteractionField init
    """
    
    # Set up what information is available for taxon field validation
    tx_obj = fixture_taxa if has_taxa_object else None
    dwsh = DataWorksheet({'name': 'DF',
                          'title': 'My data table',
                          'description': 'This is a test data worksheet'})
    dwsh.taxa_fields = ['predator', 'prey']
    dwsh_obj = dwsh if has_dwsh_object else None

    # Technically, this violates the field_name last requirement, but that is
    # enforced at the dataworksheet level, not the field level. 
    field_meta = {'field_type': 'categorical interaction',
                  'description': 'outcome of competition',
                  'field_name': 'factor_data',
                  'levels': 'level1;level2'}
    field_meta.update(iact_meta)

    caplog.clear()

    # Test the __init__ method
    fld = NumericInteractionField(field_meta,
                                  taxa = tx_obj,
                                  dwsh = dwsh_obj)

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize('data, expected_log',
    [
     (['level1', 'level2', 'level1', 'level2', 'level1', 'level2'], 
      ((INFO, "Checking Column factor_data"),)),
     (['level1', 'NA', 'level1', 'level2', 'NA', 'level2'], 
      ((INFO, "Checking Column factor_data"),
       (WARNING, "2 / 6 values missing"))),
     (['level1', '    ', 'level1', 'level2', None, 'level2'], 
      ((INFO, "Checking Column factor_data"),
       (ERROR, "2 cells are blank or contain only whitespace text"))),
     (['level1', '#REF!', 'level1', '#N/A', 'level1', 'level2'], 
      ((INFO, "Checking Column factor_data"),
       (ERROR, "2 cells contain Excel formula errors"))),
     (['level1', 'level2', 'level1', 1234, 'level1', 'level2'], 
      ((INFO, "Checking Column factor_data"),
       (ERROR, "Cells contain non-text values"))),
     (['level1', 'level2', 'NA', None, 1234, 'level2'], 
      ((INFO, "Checking Column factor_data"),
       (ERROR, "Cells contain non-text values"), 
       (WARNING, "1 / 6 values missing"),
       (ERROR, "1 cells are blank or contain only whitespace text"))),
     (['level1', 'level2', 'level3', 'level2', 'level1', 'level2'], 
      ((INFO, "Checking Column factor_data"),
       (ERROR, "Categories found in data missing from levels descriptor"))),
     (['level1', 'level1', 'level1', 'level1', 'level1', 'level1'], 
      ((INFO, "Checking Column factor_data"),
       (ERROR, "Categories found in levels descriptor not used in data"))),
    ])
@pytest.mark.parametrize('test_class, field_meta', 
  [ ( CategoricalField, 
      {'field_type': 'categorical',
       'description': 'a factor',
       'field_name': 'factor_data',
       'levels': 'level1;level2;'}),
    ( CategoricalTaxonField,
      {'field_type': 'categorical trait',
       'description': 'ant colours',
       'field_name': 'factor_data',
       'levels': 'level1;level2',
       'taxon_name': 'C_born'}),
    ( CategoricalInteractionField,
      {'field_type': 'categorical interaction',
       'description': 'outcome of competition',
       'field_name': 'factor_data',
       'levels': 'level1;level2',
       'interaction_name': 'C_born:competitor 1;V_salv:competitor 2'}),
  ]
)
def test_CategoricalField_validate_data(caplog, fixture_taxa, 
                                        test_class, field_meta, data, expected_log):
    """Testing behaviour of the CategoricalField class in using validate_data
    """

    fld = test_class(field_meta, taxa=fixture_taxa)

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
def test_TaxaField_init(caplog, fixture_taxa, provide_taxa_instance, expected_log):
    """Testing behaviour of the TaxaField class in handling missing taxa.
    """

    if provide_taxa_instance:
        taxa = fixture_taxa
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
def test_TaxaField_validate_data(caplog, fixture_taxa, data, expected_log):
    """Testing behaviour of the TaxaField class in using validate_data
    """

    fld = TaxaField({'field_type': 'taxa',
                     'description': 'My taxa',
                     'field_name': 'taxa_field'}, taxa=fixture_taxa)
    
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
def test_LocationsField_init(caplog, fixture_locations, provide_loc_instance, expected_log):
    """Testing behaviour of the LocationsField class in handling missing locations.
    """

    if provide_loc_instance:
        locs = fixture_locations
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
def test_LocationsField_validate_data(caplog, fixture_locations, data, expected_log):
    """Testing behaviour of the TaxaField class in using validate_data
    """

    fld = LocationsField({'field_name': 'locations',
                          'field_type': 'locations',
                          'description': 'my locations'}, 
                          locations=fixture_locations)

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
def test_GeoField_init(caplog, fixture_dataset, provide_ds_instance, expected_log):
    """Testing behaviour of the GeoField class in handling missing dataset.
    """

    if provide_ds_instance:
        ds = fixture_dataset
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
def test_GeoField_validate_data(caplog, fixture_dataset, data, expected_log, which):
    """Testing behaviour of the TaxaField class in using validate_data
    """

    fld = GeoField({'field_name': 'geocoords',
                    'field_type': which,
                    'description': 'my gcs'}, 
                    dataset=fixture_dataset)

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])

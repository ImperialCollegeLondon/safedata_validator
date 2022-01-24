from collections import OrderedDict
from plistlib import load
import pytest
from logging import CRITICAL, ERROR, WARNING, INFO

from safedata_validator.field import DataWorksheet


@pytest.fixture()
def fixture_field_meta():
    """field_meta object for use across tests
    """
    
    return OrderedDict(field_type = ['numeric', 'numeric', 'numeric'],
                       description = ['a', 'b', 'c'],
                       units = ['a', 'b', 'c'],
                       method = ['a', 'b', 'c'],
                       field_name = ['a', 'b', 'c'])


@pytest.fixture()
def fixture_dataworksheet(fixture_field_meta):
    """field_meta object for use across tests
    """
    
    dws = DataWorksheet({'name': 'DF',
                         'title': 'My data table',
                         'description': 'This is a test data worksheet',
                         'external': None})

    return dws


@pytest.mark.parametrize(
    'sheet_meta,expected_log',
    [
      ({'name': 'DF',
       'title': 'My data table',
       'description': 'This is a test data worksheet',
       'external': None},
      ( (INFO, "Checking data worksheet"),)),
    ]
)
def test_DataWorksheet_init(caplog, fixture_field_meta, sheet_meta, expected_log):
    """Testing init errors - field_meta tested elsewhere
    TODO - currently no errors throw  in init - but might be added
    """

    dws = DataWorksheet(sheet_meta, fixture_field_meta)

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize(
    'field_meta,expected_log',
    [
      (OrderedDict(field_type = ['numeric', 'numeric', 'numeric'],
                   description = ['a', 'b', 'c'],
                   field_name = ['a', 'b', 'c']),
        ( (INFO, 'Validating field metadata'),)),
      (OrderedDict([('field_type ',  ['numeric', 'numeric', 'numeric']),
                    (' description', ['a', 'b', 'c']),
                    ('  field_name', ['a', 'b', 'c'])]),
        ( (INFO, 'Validating field metadata'),
          (ERROR, 'Whitespace padding in descriptor names'),)),
      (OrderedDict([('field_type',  ['numeric', 'numeric', 'numeric']),
                    ('field_name', ['a', 'b', 'c'])]),
        ( (INFO, 'Validating field metadata'),
          (ERROR, 'Mandatory field descriptors missing'),)),
      (OrderedDict([('field_type',  ['numeric', 'numeric', 'numeric']),
                    ('description', ['a', 'b', 'c']),
                    ('whatdoesthisdo', ['a', 'b', 'c']),
                    ('field_name', ['a', 'b', 'c'])]),
        ( (INFO, 'Validating field metadata'),
          (WARNING, 'Unknown field descriptors'),)),
      (OrderedDict([('field_type',  ['numeric', 'numeric', 'numeric']),
                    ('field_name', ['a', 'b', 'c']),
                    ('description', ['a', 'b', 'c'])]),
        ( (INFO, 'Validating field metadata'),
          (ERROR, 'Field_name row is not the last descriptor'),)),
      (OrderedDict([('field_type',  ['numeric', 'numeric', 'numeric']),
                    ('description', ['a', 'b', 'c']),
                    ('field_name', ['a', 'b', 'a'])]),
        ( (INFO, 'Validating field metadata'),
          (ERROR, 'Field names duplicated'),)),
    ]
)
def test_DataWorksheet_validate_field_meta(caplog, fixture_dataworksheet, field_meta, expected_log):
    """Testing errors created when field metadata is validated"""
    fixture_dataworksheet.validate_field_meta(field_meta)
    
    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize(
    'load_field_meta,data_rows,expected_log',
    [
      ( True,
        [[1, 1, 2, 3],
         [2, 1, 2, 3],
         [3, 1, 2, 3], ],
        ( (INFO, 'Validating field metadata'),)),
      ( False,
        [[1, 1, 2, 3],
         [2, 1, 2, 3],
         [3, 1, 2, 3], ],
        ((CRITICAL, 'No fields defined - use validate_field_meta'),)),
      ( True,
        [[1, 1, 2, 3, 4],
         [2, 1, 2, 3],
         [3, 1, 2, 3], ],
        ( (INFO, 'Validating field metadata'),
          (ERROR, 'Data rows of unequal length - cannot validate'),)),
      ( True,
        [[1, 1, 2, 3, 4],
         [2, 1, 2, 3, 4],
         [3, 1, 2, 3, 4], ],
        ( (INFO, 'Validating field metadata'),
          (ERROR, 'Data rows not of same length as field metadata - cannot validate'),)),
    ]
)
def test_DataWorksheet_validate_data_rows(caplog, fixture_dataworksheet, fixture_field_meta,
                                          load_field_meta, data_rows, expected_log):
    """Testing errors created as data rows are added.
    """
    caplog.clear()

    if load_field_meta:
        fixture_dataworksheet.validate_field_meta(fixture_field_meta)
    
    fixture_dataworksheet.validate_data_rows(data_rows)

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])


@pytest.mark.parametrize(
    'data_rows,expected_log',
    [
      ( [[1, 1, 2, 3],
         [2, 1, 2, 3],
         [3, 1, 2, 3], ],
        ( (INFO, "Validating field data"),
          (INFO, "Checking Column a"),
          (INFO, "Checking Column b"),
          (INFO, "Checking Column c"),
          (INFO, "Worksheet 'DF' contains 6 descriptors, 0 data rows and 3 fields"),
          (INFO, "Dataframe formatted correctly"))),
      ( [[11, 1, 2, 3],
         [12, 1, 2, 3],
         [13, 1, 2, 3], ],
        ( (INFO, "Validating field data"),
          (INFO, "Checking Column a"),
          (INFO, "Checking Column b"),
          (INFO, "Checking Column c"),
          (ERROR, "Row numbers not consecutive or do not start with 1"),
          (INFO, "Worksheet 'DF' contains 6 descriptors, 0 data rows and 3 fields"),
          (INFO, "Dataframe contains 1 errors"))),
    ]
)
def test_DataWorksheet_report(caplog, fixture_dataworksheet, fixture_field_meta,
                              data_rows, expected_log):
    """Testing report level errors 
      - emission of field errors and
      - dataset level issues
    """
    
    fixture_dataworksheet.validate_field_meta(fixture_field_meta)
    caplog.clear()

    fixture_dataworksheet.validate_data_rows(data_rows)
    fixture_dataworksheet.report()

    assert len(expected_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno 
                for exp, rec in zip(expected_log, caplog.records)])
    assert all([exp[1] in rec.message
                for exp, rec in zip(expected_log, caplog.records)])




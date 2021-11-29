import pytest
import os
import string
import openpyxl

from safedata_validator.logger import LOGGER
from safedata_validator.summary import *

# @pytest.fixture()
# def summary_ws():
#
#     wb = openpyxl.Workbook()
#     ws = wb.active
#
#     # Create a 11 by 38 simple block of data with A-K column headers
#     # containing a multiplication table.
#     for col in range(1, 12):
#         ws.cell(column=col, row=1).value = openpyxl.utils.get_column_letter(col)
#         for row in range(2, 40):
#             ws.cell(column=col, row=row).value = (row - 1) * col
#
#     # Put in an empty cell at 20, 50 which leads to the kind of extra rows
#     # and columns that often occur with Excel
#     ws.cell(column=20, row=50).value = None
#
#     return ws


# TODO - _read_block is being tested by repeated _read_`block` calls and
#        this could be simplified to reduce the code, but there are block
#        specific variations (types, mandatory etc) that make it easier
#        to test block by block

@pytest.mark.parametrize('alterations,should_log_error,expected_log', [
    (dict(),  # no amendments
     False,
     'Metadata for Authors found'),
    ({'author name': (None,)},
     True,
     'Missing metadata in mandatory field author name'),  # via _read_block
    ({'author name': ('   ',)},
     True,
     'Whitespace only cells in field author name'),  # via _read_block
    ({'author name': (123456789,)},
     True,
     'Field author name contains values of wrong type'),   # via _read_block
    ({'author name': ('David Orme',)},
     True,
     'Author names not formatted as last_name, first_names'),
    ({'author affiliation': (None,)},
     False,  # missing affiliation acceptable
     'Metadata for Authors found'),
    ({'author affiliation': ('   ',)},
     True,
     'Whitespace only cells in field author affiliation'),  # via _read_block
    ({'author affiliation': (123456789,)},
     True,
     'Field author affiliation contains values of wrong type'),  # via _read_block
    ({'author email': (None,)},
     False,  # missing email acceptable
     'Metadata for Authors found'),
    ({'author email': ('   ',)},
     True,
     'Whitespace only cells in field author email'),  # via _read_block
    ({'author email': (123456789,)},
     True,
     'Field author email contains values of wrong type'),  # via _read_block
    ({'author email': ('d.ormeatimperial.ac.uk',)},
     True,
     'Author emails not properly formatted'),
    ({'author orcid': (None,)},
     False,  # Blank ORCID is OK
     'Metadata for Authors found'),
    ({'author orcid': ('   ',)},
     True,
     'Whitespace only cells in field author orcid'),  # via _read_block
    ({'author orcid': (123456789,)},
     True,
     'Field author orcid contains values of wrong type'),  # via _read_block
    ({'author orcid': ('123456789',)},
     True,
     'Author ORCIDs not properly formatted')
    ])
def test_authors(caplog, alterations, should_log_error, expected_log):

    # Initialise a Summary instance.
    summary = Summary(None)

    # Valid set of information
    input = {'author name': ('Orme, David',),
             'author email': ('d.orme@imperial.ac.uk',),
             'author affiliation': ('Imperial College London',),
             'author orcid': ('0000-0002-7005-1394',)}

    # Update valid to test error conditions and populate _rows
    # directly (bypassing .load() and need to pack in worksheet object
    input.update(alterations)
    summary._rows = input

    # Test the block load
    summary._load_authors()

    if should_log_error:
        assert 'ERROR' in [r.levelname for r in caplog.records]

    assert expected_log in caplog.text


@pytest.mark.parametrize('alterations,should_log_error,expected_log', [
    (dict(),  # no amendments,
     False,
     'Metadata for Keywords found'),
    ({'keywords': (None, None, None)},
     True,
     'No Keywords metadata found'),  # via _read_block
    ({'keywords': ('   ', 'abc', 'def')},
     True,
     'Whitespace only cells in field keywords'),  # via _read_block
    ({'keywords': (123456789, 12.0)},
     True,
     'Field keywords contains values of wrong type'),  # via _read_block
    ({'keywords': ('abc,def', 'ghi;jkl')},
     True,
     'Put each keyword in a separate cell'),
    ])
def test_keywords(caplog, alterations, should_log_error, expected_log):

    # Initialise a Summary instance.
    summary = Summary(None)

    # Valid set of information
    input = {'keywords': ('abc', 'def')}

    # Update valid to test error conditions and populate _rows
    # directly (bypassing .load() and need to pack in worksheet object
    input.update(alterations)
    summary._rows = input

    # Test the block load
    summary._load_keywords()

    if should_log_error:
        assert 'ERROR' in [r.levelname for r in caplog.records]

    assert expected_log in caplog.text


@pytest.mark.parametrize('alterations,should_log_error,expected_log', [
    (dict(),  # no amendments,
     False,
     'Metadata for Permits found'),
    ({'permit type': (None,),
      'permit authority': (None,),
      'permit number': (None,)},
     False,  # Optional block
     ''),
    ({'permit type': (None,)},
     True,
     'Missing metadata in mandatory field permit type'),  # via _read_block
    ({'permit type': ('   ',)},
     True,
     'Whitespace only cells in field permit type'),  # via _read_block
    ({'permit type': (123456789,)},
     True,
     'Field permit type contains values of wrong type'),  # via _read_block
    ({'permit authority': (None,)},
     True,
     'Missing metadata in mandatory field permit authority'),  # via _read_block
    ({'permit authority': ('   ',)},
     True,
     'Whitespace only cells in field permit authority'),  # via _read_block
    ({'permit authority': (123456789,)},
     True,
     'Field permit authority contains values of wrong type'),  # via _read_block
    ({'permit number': (None,)},
     True,
     'Missing metadata in mandatory field permit number'),  # via _read_block
    ({'permit number': ('   ',)},
     True,
     'Whitespace only cells in field permit number'),  # via _read_block
    ({'permit number': (123456789,)},
     False,  # Allow numeric numbers!
     'Metadata for Permits found'),
    ({'permit number': (datetime.date(1901, 1, 1),)},
     True,
     'Field permit number contains values of wrong type'),  # via _read_block
    ])
def test_permits(caplog, alterations, should_log_error, expected_log):

    # Initialise a Summary instance.
    summary = Summary(None)

    # Valid set of information
    input = {'permit type': ('research',),
             'permit authority': ('Sabah Biodiversity Centre',),
             'permit number': ('ABC-123-456',)}

    # Update valid to test error conditions and populate _rows
    # directly (bypassing .load() and need to pack in worksheet object
    input.update(alterations)
    summary._rows = input

    # Test the block load
    summary._load_permits()

    if should_log_error:
        assert 'ERROR' in [r.levelname for r in caplog.records]

    assert expected_log in caplog.text


@pytest.mark.parametrize('alterations,should_log_error,expected_log,do_val_doi', [
    # (dict(),  # no amendments,
    #  False,
    #  'Metadata for DOI found',
    #  False),
    # ({'publication doi': (None,)},
    #  False,  # Optional block
    #  '',
    #  False),
    # ({'publication doi': ('   ',)},
    #  True,
    #  'Whitespace only cells in field publication doi',
    #  False),  # via _read_block
    # ({'publication doi': (123456789,)},
    #  True,
    #  'Field publication doi contains values of wrong type',
    #  False),  # via _read_block
    ({'publication doi': ('thisisnotaURL',)},
     True,
     'Publication DOIs not all in format',
     False),
    ({'publication doi': ('https://doi.org/this.does/not.exist',)},
     True,
     'DOI not found',
     True),
    ])
def test_doi(caplog, alterations, should_log_error, expected_log, do_val_doi):

    # Initialise a Summary instance.
    summary = Summary(None, validate_doi=do_val_doi)

    # Valid set of information
    input = {'publication doi': ('https://doi.org/10.1098/rstb.2011.0049',)}

    # Update valid to test error conditions and populate _rows
    # directly (bypassing .load() and need to pack in worksheet object
    input.update(alterations)
    summary._rows = input

    # Test the block load
    summary._load_doi()

    if should_log_error:
        assert 'ERROR' in [r.levelname for r in caplog.records]

    assert expected_log in caplog.text


@pytest.mark.parametrize('alterations,should_log_error,expected_log', [
    (dict(),  # no amendments
     False,
     'Metadata for Funding Bodies found'),
    ({'funding body': (None,),
      'funding type': (None,),
      'funding reference': (None,),
      'funding link': (None,)},
     False,  # Not a mandatory block
     ''),
    ({'funding body': (None,)},
     True,
     'Missing metadata in mandatory field funding body'),  # via _read_block
    ({'funding body': ('   ',)},
     True,
     'Whitespace only cells in field funding body'),  # via _read_block
    ({'funding body': (123456789,)},
     True,
     'Field funding body contains values of wrong type'),   # via _read_block
    ({'funding type': (None,)},
     True,
     'Missing metadata in mandatory field funding type'),  # via _read_block
    ({'funding type': ('   ',)},
     True,
     'Whitespace only cells in field funding type'),  # via _read_block
    ({'funding type': (123456789,)},
     True,
     'Field funding type contains values of wrong type'),   # via _read_block
    ({'funding reference': (None,)},
     False,
     'Metadata for Funding Bodies found'),  # Not mandatory
    ({'funding reference': ('   ',)},
     True,
     'Whitespace only cells in field funding reference'),  # via _read_block
    ({'funding reference': (123456789,)},
     False,
     'Metadata for Funding Bodies found'),   # Could be numeric
    ({'funding reference': (datetime.date(2011, 1, 1),)},
     True,
     'Field funding reference contains values of wrong type'),   # via _read_block
    ({'funding link': (None,)},
     False,
     'Metadata for Funding Bodies found'),  # Not mandatory
    ({'funding link': ('   ',)},
     True,
     'Whitespace only cells in field funding link'),  # via _read_block
    ({'funding link': (123456789,)},
     True,
     'Field funding link contains values of wrong type'),  # via _read_block
])
def test_funders(caplog, alterations, should_log_error, expected_log):

    # Initialise a Summary instance.
    summary = Summary(None)

    # Valid set of information
    input = {'funding body': ('NERC',),
             'funding type': ('Standard grant',),
             'funding reference': ('NE/K006339/1',),
             'funding link': ('https://gtr.ukri.org/projects?ref=NE%2FK006339%2F1',)}

    # Update valid to test error conditions and populate _rows
    # directly (bypassing .load() and need to pack in worksheet object
    input.update(alterations)
    summary._rows = input

    # Test the block load
    summary._load_funders()

    if should_log_error:
        assert 'ERROR' in [r.levelname for r in caplog.records]

    assert expected_log in caplog.text


@pytest.mark.parametrize('alterations,should_log_error,expected_log', [
    (dict(),  # no amendments
     False,
     'Metadata for Date Extents found'),
    ({'start date': (None,),
      'end date': (None,)},
     False,  # Not a mandatory block
     ''),
    ({'start date': (None,)},
     True,
     'Missing metadata in mandatory field start date'),  # via _read_block
    ({'start date': ('   ',)},
     True,
     'Whitespace only cells in field start date'),  # via _read_block
    ({'start date': ('2001-01-01',)},
     True,
     'Field start date contains values of wrong type'),  # via _read_block
    ({'end date': (None,)},
     True,
     'Missing metadata in mandatory field end date'),  # via _read_block
    ({'end date': ('   ',)},
     True,
     'Whitespace only cells in field end date'),  # via _read_block
    ({'end date': ('2001-01-01',)},
     True,
     'Field end date contains values of wrong type'),  # via _read_block
    ({'start date': (datetime.datetime(2001, 1, 1), datetime.datetime(2001, 1, 1)),
      'end date': (datetime.datetime(2011, 1, 1), datetime.datetime(2011, 1, 1))},
     True,
     'Only a single record should be present'),
    ({'start date': (datetime.datetime(2011, 1, 1),),
      'end date': (datetime.datetime(2001, 1, 1),)},
     True,
     'Start date is after end date'),
])
def test_temporal_extent(caplog, alterations, should_log_error, expected_log):

    # Initialise a Summary instance.
    summary = Summary(None)

    # Valid set of information - openpyxl loads dates as datetimes.
    input = {'start date': (datetime.datetime(2001, 1, 1),),
             'end date': (datetime.datetime(2011, 1, 1),)}

    # Update valid to test error conditions and populate _rows
    # directly (bypassing .load() and need to pack in worksheet object
    input.update(alterations)
    summary._rows = input

    # Test the block load
    summary._load_temporal_extent()

    if should_log_error:
        assert 'ERROR' in [r.levelname for r in caplog.records]

    assert expected_log in caplog.text




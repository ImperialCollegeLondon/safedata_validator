import pytest
from safedata_validator.field import DataWorksheet


@pytest.mark.parametrize(
    'example_excel_files, worksheet, sheet_meta, n_log_records',
    [('good', 'DF', 
      {'name': 'DF',
       'title': 'My data table',
       'description': 'This is a test data worksheet',
       'external': None},
      0), 
     ('bad', 'DF', 
      {'name': 'DF',
       'title': 'My data table',
       'description': 'This is a test data worksheet',
       'external': None},
      1)], 
    indirect = ['example_excel_files']  # take actual params from fixture
)
def test_dataworksheet_load_meta_and_size(
        caplog, example_excel_files, worksheet, sheet_meta, n_log_records):

    ws = example_excel_files[worksheet]
    dws = DataWorksheet(sheet_meta, ws)

    dws._load_meta_and_size()

    # TODO - assert something a bit more meaningful!
    assert len(caplog.records) == n_log_records
import pytest
from safedata_validator.field import Dataset
from .conftest import FIXTURE_FILES

@pytest.mark.parametrize(
    'file_key, n_errors',
    [
      ('good_excel_file' , 0), 
      ('bad_excel_file', 94), 
    ], 
)
def test_DataSet_load_from_file(config_filesystem,
                                resources_with_local_gbif, file_key, n_errors):
    """Test loading a dataworksheet from file - this duplicates a lot of 
    Dataset.load_from_workbook"""
    
    # Load the taxa and locations
    ds = Dataset(resources_with_local_gbif)
    ds.load_from_workbook(FIXTURE_FILES.rf[file_key])

    assert ds.n_errors == n_errors
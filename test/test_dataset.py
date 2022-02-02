import pytest
from safedata_validator.field import Dataset
from .conftest import good_file_path, bad_file_path

@pytest.mark.parametrize(
    'filepath, n_errors',
    [
      (good_file_path , 0), 
      (bad_file_path, 0), 
    ], 
)
def test_DataWorksheet_load_from_file(caplog, resources_with_local_gbif,
                                      filepath, n_errors):
    """Test loading a dataworksheet from file - this duplicates a lot of 
    Dataset.load_from_workbook"""
    
    # Load the taxa and locations
    ds = Dataset(resources_with_local_gbif)
    ds.load_from_workbook(filepath)

    assert ds.n_errors == n_errors
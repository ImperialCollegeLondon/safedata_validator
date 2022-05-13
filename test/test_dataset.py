import pytest
from safedata_validator.field import Dataset
from .conftest import FIXTURE_FILES

@pytest.mark.parametrize(
    'file_key, n_errors',
    [
      ('good_excel_file' , 0),
      ('bad_excel_file', 94),
      ('good_ncbi_file', 0),
      ('bad_ncbi_file', 107)
    ],
)
def test_DataSet_load_from_file(ncbi_resources_local_and_remote, file_key, n_errors):
    """Test loading a dataworksheet from file - this duplicates a lot of
    Dataset.load_from_workbook"""

    # Load the taxa and locations
    ds = Dataset(ncbi_resources_local_and_remote)
    ds.load_from_workbook(FIXTURE_FILES.rf[file_key])

    assert ds.n_errors == n_errors

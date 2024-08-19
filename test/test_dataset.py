"""Test loading a dataworksheet from file."""

import pytest

from safedata_validator.field import Dataset

from .conftest import FIXTURE_FILES


@pytest.mark.parametrize(
    "file_key, n_errors",
    [
        ("good_excel_file", 0),
        ("bad_excel_file", 94),
        ("good_ncbi_file", 0),
        ("bad_ncbi_file", 106),
    ],
)
def test_DataSet_load_from_file(fixture_resources, file_key, n_errors):
    # Load the taxa and locations
    ds = Dataset(fixture_resources)
    ds.load_from_workbook(FIXTURE_FILES.rf[file_key])

    assert ds.n_errors == n_errors


def test_example_dataset(fixture_resources):
    """Test that the example dataset we provide actually passes validation.

    It is important to note that this will only pass if the truncated databases are made
    using the database versions stated in their respective *_database_truncator.py
    scripts. Otherwise changes to the reference taxonomy could cause this test to fail.

    This test will also fail every two years when the embargo date in the example files
    comes around.
    """
    # Load the taxa and locations
    ds = Dataset(fixture_resources)
    ds.load_from_workbook(FIXTURE_FILES.rf["example_file"])

    assert ds.n_errors == 0

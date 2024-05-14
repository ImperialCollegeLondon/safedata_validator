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


# TODO - Work out if this requires a specific database version
def test_Example_dataset(fixture_resources):
    # Load the taxa and locations
    ds = Dataset(fixture_resources)
    ds.load_from_workbook(FIXTURE_FILES.rf["example_file"])

    assert ds.n_errors == 0

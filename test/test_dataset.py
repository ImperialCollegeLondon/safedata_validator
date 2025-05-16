"""Test loading a dataworksheet from file."""

import pytest

from safedata_validator.field import Dataset

from .conftest import FIXTURE_FILES


@pytest.mark.parametrize(
    "file_key, n_errors, n_taxa",
    [
        ("good_excel_file", 0, 20),
        ("bad_excel_file", 94, 13),
        ("good_seq_taxa_file", 0, 48),
        ("bad_seq_taxa_file", 6, 34),
    ],
)
def test_DataSet_load_from_file(fixture_resources, file_key, n_errors, n_taxa):
    # Load the taxa and locations
    ds = Dataset(fixture_resources)
    ds.load_from_workbook(FIXTURE_FILES.rf[file_key])

    assert ds.n_errors == n_errors
    assert len(ds.taxa.taxon_names) == n_taxa


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

"""Tests for the utility functions."""

from pathlib import Path

import pytest

from .conftest import FIXTURE_FILES


@pytest.mark.parametrize(
    "file_path,expected_result",
    [
        (Path(FIXTURE_FILES.rf.good_excel_file), True),
        (Path(FIXTURE_FILES.rf.good_ncbi_file_dataset_json), False),
    ],
)
def test_check_file_is_excel(file_path, expected_result):
    """Test that the function to check if a file is .xlsx works."""
    from safedata_validator.utilities import check_file_is_excel

    is_file_xlsx = check_file_is_excel(file_path)

    assert is_file_xlsx == expected_result


@pytest.mark.parametrize(
    "file_path,expected_result",
    [
        (Path(FIXTURE_FILES.rf.good_ncbi_file_dataset_json), True),
        (Path(FIXTURE_FILES.rf.good_excel_file), False),
    ],
)
def test_check_file_is_json(file_path, expected_result):
    """Test that the function to check if a file is JSON works."""
    from safedata_validator.utilities import check_file_is_json

    is_file_json = check_file_is_json(file_path)

    assert is_file_json == expected_result

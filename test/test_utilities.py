"""Tests for the utility functions."""

from pathlib import Path

import pytest

from .conftest import FIXTURE_FILES


@pytest.mark.parametrize(
    "file_path,expected_result",
    [
        (Path(FIXTURE_FILES.rf.good_excel_file), True),
        (Path(FIXTURE_FILES.rf.good_ncbi_file_dataset_json), False),
        (Path("safedata_validator"), False),
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
        (Path("safedata_validator"), False),
    ],
)
def test_check_file_is_json(file_path, expected_result):
    """Test that the function to check if a file is JSON works."""
    from safedata_validator.utilities import check_file_is_json

    is_file_json = check_file_is_json(file_path)

    assert is_file_json == expected_result


# PROVIDE EXCEL FILE AND ALL THREE JSON FILES, NEED TO TRACK DOWN NAMES
@pytest.mark.parametrize(
    "file_path,expected_result",
    [
        (Path(FIXTURE_FILES.rf.good_ncbi_file_dataset_json), False),
        (Path(FIXTURE_FILES.rf.good_ncbi_file_zenodo_json), True),
        (Path(FIXTURE_FILES.rf.json_not_locations), False),
        (Path(FIXTURE_FILES.rf.good_excel_file), False),
    ],
)
def test_check_file_is_zenodo_json(file_path, expected_result):
    """Test that the function to check if a file is JSON works."""
    from safedata_validator.utilities import check_file_is_zenodo_json

    is_file_json = check_file_is_zenodo_json(file_path)

    assert is_file_json == expected_result


@pytest.mark.parametrize(
    "file_path,expected_result",
    [
        (Path(FIXTURE_FILES.rf.good_ncbi_file_dataset_json), True),
        (Path(FIXTURE_FILES.rf.good_ncbi_file_zenodo_json), False),
        (Path(FIXTURE_FILES.rf.json_not_locations), False),
        (Path(FIXTURE_FILES.rf.good_excel_file), False),
    ],
)
def test_check_file_is_metadata_json(file_path, expected_result):
    """Test that the function to check if a file is JSON works."""
    from safedata_validator.utilities import check_file_is_metadata_json

    is_file_json = check_file_is_metadata_json(file_path)

    assert is_file_json == expected_result

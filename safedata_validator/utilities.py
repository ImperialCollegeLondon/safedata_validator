"""This module provides utility functions.

These utility functions are used by the high level functions defined in
`safedata_validator.zenodo` and `safedata_validator.server`. The functions are used to
check that the correct inputs have been supplied in the correct order.
"""

import mimetypes
from pathlib import Path

import simplejson


def check_file_is_excel(file_path: Path) -> bool:
    """Check that the file path provided points to a .xlsx file.

    Args:
        file_path: A path to the file to check

    Returns:
        A bool specifying whether or not the file is an .xlsx file.
    """

    file_type, _ = mimetypes.guess_type(file_path, strict=True)

    if file_type is None:
        return False
    elif (
        file_type == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    ):
        return True
    else:
        return False


def check_file_is_json(file_path: Path) -> bool:
    """Check that the file path provided points to a JSON file.

    Args:
        file_path: A path to the file to check

    Returns:
        A bool specifying whether or not the file is a JSON file.
    """

    file_type, _ = mimetypes.guess_type(file_path, strict=True)

    if file_type is None:
        return False

    extension = mimetypes.guess_extension(file_type, strict=True)

    if extension == ".json":
        return True
    else:
        return False


def check_file_is_zenodo_json(file_path: Path) -> bool:
    """Function to check if a file is a Zenodo JSON file.

    Args:
        file_path: A path to the file to check

    Returns:
        A bool specifying whether or not the file is a Zenodo JSON file.
    """

    # First check that the file is actually JSON
    if not check_file_is_json(file_path):
        return False

    # Minimal set of keys to identify zenodo JSON
    required_keys = ["created", "modified", "conceptrecid", "metadata"]

    with open(file_path) as file:
        contents = simplejson.load(file)

        # Check if all required keys are present
        missing_keys = [key for key in required_keys if key not in contents]

        if missing_keys:
            return False
        else:
            return True


def check_file_is_metadata_json(file_path: Path) -> bool:
    """Function to check if a file is a safedata metadata style JSON file.

    Args:
        file_path: A path to the file to check

    Returns:
        A bool specifying whether or not the file is a safedata metadata style JSON
        file.
    """

    # First check that the file is actually JSON
    if not check_file_is_json(file_path):
        return False

    # Minimal set of keys to identify metadata JSON
    required_keys = ["project_ids", "authors", "access", "funders"]

    with open(file_path) as file:
        contents = simplejson.load(file)

        # Check if all required keys are present
        missing_keys = [key for key in required_keys if key not in contents]

        if missing_keys:
            return False
        else:
            return True

"""This module provides utility functions.

These utility functions are used by the high level functions defined in
`safedata_validator.zenodo` and `safedata_validator.server`. The functions are used to
check that the correct inputs have been supplied in the correct order.
"""

from pathlib import Path

import simplejson


def check_file_is_excel(file_path: Path) -> bool:
    """Check that the file path provided points to a .xlsx file.

    Args:
        file_path: A path to the file to check

    Returns:
        A bool specifying whether or not the file is an .xlsx file.
    """

    if not file_path.is_file():
        return False

    if file_path.suffix != ".xlsx":
        return False

    return True


def check_file_is_json(file_path: Path) -> bool:
    """Check that the file path provided points to a JSON file.

    Args:
        file_path: A path to the file to check

    Returns:
        A bool specifying whether or not the file is a JSON file.
    """

    if not file_path.is_file():
        return False

    if file_path.suffix != ".json":
        return False

    return True


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

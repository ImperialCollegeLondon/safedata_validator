"""This module provides utility functions.

These utility functions are used by the high level functions defined in
`safedata_validator.zenodo` and `safedata_validator.server`. The functions are used to
check that the correct inputs have been supplied in the correct order.
"""

import mimetypes
from pathlib import Path


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

    extension = mimetypes.guess_extension(file_type, strict=True)

    if extension == ".xlsx":
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


# TODO - Function to distinguish dataset from zenodo json files (think this requires
# looking at tags)

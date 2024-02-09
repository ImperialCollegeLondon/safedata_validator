"""Test that Extent logs correctly."""

import importlib
import subprocess
from contextlib import contextmanager
from logging import CRITICAL, ERROR, WARNING

import pytest

from safedata_validator.extent import Extent

from .conftest import log_check


@contextmanager
def does_not_raise():
    yield


@pytest.mark.parametrize(
    argnames="entry_point",
    argvalues=[
        "safedata_validate",
        "safedata_zenodo",
        "safedata_metadata",
        "safedata_build_local_gbif",
        "safedata_build_local_ncbi",
    ],
)
def test_entry_points_run(entry_point):
    """This test checks the availability and basic function of the package entry points.

    The tests do not test complex function, just that the entry point is available and
    exits normally after getting the help text. For some reason, these tests fail in
    VSCode debug mode, something to do with multiprocessing and shell=True?
    """

    result = subprocess.run([entry_point, "-h"])

    assert result.returncode == 0


@pytest.mark.parametrize(
    argnames="entry_point",
    argvalues=[
        "_safedata_validator_cli",
        "_safedata_zenodo_cli",
        "_safedata_metadata_cli",
        "_build_local_gbif_cli",
        "_build_local_ncbi_cli",
    ],
)
def test_entry_points_via_call(entry_point):
    """This test checks that the entry points can all be invoked directly.

    The tests do not test complex function, just that the underlying entry point
    function can be called with a list of arguments. The test asks the function to print
    the function help, which results in a SystemExit that needs to be handled.
    """

    from safedata_validator import entry_points

    entry_point_function = getattr(entry_points, entry_point)

    with pytest.raises(SystemExit) as exit:
        entry_point_function(args_list=["-h"])

    assert exit.value.code == 0

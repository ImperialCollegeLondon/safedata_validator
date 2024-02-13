"""Test that Extent logs correctly."""

import subprocess
from contextlib import contextmanager

import pytest

from .conftest import FIXTURE_FILES


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

    The tests do not test complex functionality, just that the entry point is available
    and exits normally after getting the help text. For some reason, these tests fail in
    VSCode debug mode.
    """

    result = subprocess.run(f"{entry_point} -h", shell=True)

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


@pytest.mark.parametrize(
    argnames="command, file_exists, returns",
    argvalues=[
        pytest.param("generate_xml", False, 0, id="xml_ok"),
        pytest.param("generate_html", False, 0, id="html_ok"),
        pytest.param("generate_xml", True, 1, id="html_overwrite"),
        pytest.param("generate_html", True, 1, id="xml_overwrite"),
    ],
)
def test_sdv_zenodo_html_and_xml(user_config_file, command, file_exists, returns):
    """Checks that the safedata_zenodo XML and HTML generation commands work."""

    from safedata_validator.entry_points import _safedata_zenodo_cli

    if file_exists:
        with open("TMP", "w") as outfile:
            outfile.write("Don't overwrite me!")

    value = _safedata_zenodo_cli(
        args_list=[
            command,
            FIXTURE_FILES.rf.good_ncbi_file_zenodo_json,
            FIXTURE_FILES.rf.good_ncbi_file_dataset_json,
            "TMP",
        ]
    )

    assert value == returns

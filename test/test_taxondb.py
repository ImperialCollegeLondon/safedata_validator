"""Tests to check that taxondb functions are working as intended."""

from contextlib import nullcontext as does_not_raise

import pytest


@pytest.mark.parametrize(
    argnames="timestamp, raises",
    argvalues=[
        pytest.param(None, does_not_raise(), id="not provided"),
        pytest.param("2022-11-23", does_not_raise(), id="valid provided"),
        pytest.param("2022-11-29", pytest.raises(ValueError), id="unknown provided"),
        pytest.param("2022-29-23", pytest.raises(ValueError), id="invalid provided"),
    ],
)
def test_get_gbif_version(timestamp, raises):

    from safedata_validator.taxondb import get_gbif_version

    with raises:
        val = get_gbif_version(timestamp=timestamp)

        assert isinstance(val[0], str)
        assert isinstance(val[1], str)


@pytest.mark.parametrize(
    argnames="timestamp, raises",
    argvalues=[
        pytest.param(None, does_not_raise(), id="not provided"),
        pytest.param("2018-07-01", does_not_raise(), id="valid provided"),
        pytest.param("2001-01-01", pytest.raises(ValueError), id="unknown provided"),
        pytest.param("2022-May-05", pytest.raises(ValueError), id="invalid provided"),
    ],
)
def test_get_ncbi_version(timestamp, raises):

    from safedata_validator.taxondb import get_ncbi_version

    with raises:
        val = get_ncbi_version(timestamp=timestamp)

        assert isinstance(val[0], str)
        assert isinstance(val[1], str)
        assert isinstance(val[2], int)

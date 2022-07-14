"""Test that Extent logs correctly."""
from contextlib import contextmanager
from logging import CRITICAL, ERROR, INFO, WARNING

import pytest

from safedata_validator.extent import Extent

from .conftest import log_check


@contextmanager
def does_not_raise():
    yield


@pytest.mark.parametrize(
    "datatype, hbounds, sbounds, exp_exception, expected_log",
    [
        ((int,), (0, 120), (45, 60), does_not_raise(), tuple()),
        (
            (int,),
            ("a", 120),
            (45, 60),
            pytest.raises(AttributeError),
            ((CRITICAL, "Bounds are not all of type int"),),
        ),
        (
            (int,),
            (120, 0),
            (45, 60),
            pytest.raises(AttributeError),
            ((CRITICAL, "Bounds must be provided as (low, high) tuples"),),
        ),
        (
            (int,),
            (45, 60),
            (120, 0),
            pytest.raises(AttributeError),
            ((CRITICAL, "Hard bounds must lie outside soft bounds"),),
        ),
    ],
)
def test_extent_init(caplog, datatype, hbounds, sbounds, exp_exception, expected_log):

    with exp_exception:

        Extent(
            label="test", datatype=datatype, hard_bounds=hbounds, soft_bounds=sbounds
        )

        log_check(caplog, expected_log)


@pytest.mark.parametrize(
    "values, expected_log",
    [
        ([47, 49, 54, 60], tuple()),
        ([47, 49, 54.1, 60], ((ERROR, "Values are not all of type"),)),
        ([], ((ERROR, "No valid data in extent update"),)),
        (
            ["this", "is", "bad", "data", 2.0],
            (
                (ERROR, "Values are not all of type"),
                (ERROR, "No valid data in extent update"),
            ),
        ),
        ([-150, 47, 49, 54, 60], ((ERROR, "exceeds hard bounds"),)),
        ([15, 47, 49, 54, 60], ((WARNING, "exceeds soft bounds"),)),
    ],
)
def test_extent_update(caplog, values, expected_log):

    ext = Extent(
        label="test", datatype=(int,), hard_bounds=(0, 120), soft_bounds=(45, 60)
    )

    ext.update(values)
    log_check(caplog, expected_log)

"""Test the basic logging formatting."""

from logging import CRITICAL, DEBUG, ERROR, INFO, WARNING

import pytest

from safedata_validator.logger import LOGGER


@pytest.mark.parametrize(
    "level, message, extra, expected",
    [
        (INFO, "Test", {}, "- Test"),
    ],
)
def test_logger(caplog, level, message, extra, expected):
    """Testing the logger formatting behaviour."""

    LOGGER.log(level, message, extra=extra)

    assert caplog.records[0].message == expected

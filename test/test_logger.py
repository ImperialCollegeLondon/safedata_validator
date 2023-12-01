"""Test the basic logging formatting."""

from io import StringIO
from logging import CRITICAL, DEBUG, ERROR, INFO, WARNING, StreamHandler

import pytest

from safedata_validator.logger import FORMATTER, LOGGER


@pytest.mark.parametrize(
    "level, message, extra, depth, expected",
    [
        (DEBUG, "Test", {}, 0, "> Test"),
        (DEBUG, "Test: ", {"join": [1, 2, 3]}, 0, "> Test: 1, 2, 3"),
        (DEBUG, "Test", {}, 1, "    > Test"),
        (DEBUG, "Test: ", {"join": [1, 2, 3]}, 1, "    > Test: 1, 2, 3"),
        (DEBUG, "Test", {}, 2, "        > Test"),
        (DEBUG, "Test: ", {"join": [1, 2, 3]}, 2, "        > Test: 1, 2, 3"),
        (INFO, "Test", {}, 0, "- Test"),
        (INFO, "Test: ", {"join": [1, 2, 3]}, 0, "- Test: 1, 2, 3"),
        (WARNING, "Test", {}, 0, "? Test"),
        (WARNING, "Test: ", {"join": [1, 2, 3]}, 0, "? Test: 1, 2, 3"),
        (ERROR, "Test", {}, 0, "! Test"),
        (ERROR, "Test: ", {"join": [1, 2, 3]}, 0, "! Test: 1, 2, 3"),
        (CRITICAL, "Test", {}, 0, "X Test"),
        (CRITICAL, "Test: ", {"join": [1, 2, 3]}, 0, "X Test: 1, 2, 3"),
    ],
)
def test_logger(caplog, level, message, extra, depth, expected):
    """Testing the logger formatting behaviour."""

    # Capture _all_ logging events, not just WARNING and above
    caplog.set_level(DEBUG)

    # Truncate any existing log messages and reset the depth
    FORMATTER.depth = 0

    # Capture the formatted output as text
    log = StringIO()
    handler = StreamHandler(log)
    handler.setFormatter(FORMATTER)

    # Emit the message
    LOGGER.addHandler(handler)
    FORMATTER.push(depth)
    LOGGER.log(level, message, extra=extra)
    FORMATTER.pop(depth)

    # Look at the stdout - the value is accessed _directly_ from the LOG StringIO object
    # used by the CONSOLE_HANDLER. In theory, the capfd/capsys fixtures should be able
    # to recover the formatted text written to the command line, but they do not:
    # https://github.com/pytest-dev/pytest/issues/5997

    assert log.getvalue() == (expected + "\n")

    LOGGER.removeHandler(handler)

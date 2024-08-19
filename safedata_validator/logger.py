"""This submodule extends the standard logging setup to provide extra functionality
and to expose some global logging objects for use throughout the code.

1. The `logging.LogRecordFactory` is updated so that new records include a custom
   `levelcode` attribute to visually indicate log record severity in validation
   reports.

2. The [IndentFormatter][safedata_validator.logger.IndentFormatter] class then extends
   :class:`logging.Formatter` to provide compact messages with variable indentation to
   show nested sections of the validation process using the level codes as visual cues
   for problems.

3. The submodule then defines two `CounterHandler` classes which subclass
   `logging.StreamHandler` and `logging.FileHandler`. Both extend the basic handlers to
   add attributes that keep track of counts of different classes of records emitted
   through the handler.

4. The submodule provides the functions
   [use_file_logging][safedata_validator.logger.use_file_logging] and
   [use_stream_logging][safedata_validator.logger.use_stream_logging] to assign
   handlers to be used in the validation process. The
   [get_handler][safedata_validator.logger.get_handler] function is then used as a
   convenience function to retrieve the current handler to access counts of the various
   emitted records.

5. The functions [log_and_raise][safedata_validator.logger.log_and_raise] and
   [loggerinfo_push_pop][safedata_validator.logger.loggerinfo_push_pop] are convenience
   functions to minimise logging boilerplate code within the package.
"""  # noqa D415

import logging
from collections.abc import Callable
from functools import wraps
from pathlib import Path
from typing import Any

LOGGER_CODES = {
    "DEBUG": ">",
    "INFO": "-",
    "WARNING": "?",
    "ERROR": "!",
    "CRITICAL": "X",
}


# Extend the logging record factory to include custom attributes. This could also be
# achieved by changing all the level names: https://stackoverflow.com/questions/47035760
old_factory = logging.getLogRecordFactory()


def record_factory(*args, **kwargs):
    """Logging record factory with extended attributes.

    The representation of logging records uses single character codes to show logging
    levels. This record factory adds custom attributes to records to
    support that behaviour.

    """
    record = old_factory(*args, **kwargs)
    record.levelcode = LOGGER_CODES[record.levelname]

    return record


logging.setLogRecordFactory(record_factory)


class StreamCounterHandler(logging.StreamHandler):
    """Subclass of `logging.StreamHandler` counting calls emitted at each log level."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.counters = {"DEBUG": 0, "INFO": 0, "WARNING": 0, "ERROR": 0, "CRITICAL": 0}

    def emit(self, record: logging.LogRecord) -> None:
        """Emit a message and increment the counter for the message level.

        Args:
            record: A `logging.LogRecord` instance.
        """
        self.counters[record.levelname] += 1
        super().emit(record=record)

    def reset(self) -> None:
        """Reset the message counters to zero."""
        self.counters = {"DEBUG": 0, "INFO": 0, "WARNING": 0, "ERROR": 0, "CRITICAL": 0}


class FileCounterHandler(logging.FileHandler):
    """Subclass of `logging.FileHandler` counting calls emitted at each log level."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.counters = {"DEBUG": 0, "INFO": 0, "WARNING": 0, "ERROR": 0, "CRITICAL": 0}

    def emit(self, record: logging.LogRecord) -> None:
        """Emit a message and increment the counter for the message level.

        Args:
            record: A `logging.LogRecord` instance.
        """
        self.counters[record.levelname] += 1
        super().emit(record=record)

    def reset(self) -> None:
        """Reset the message counters to zero."""
        self.counters = {"DEBUG": 0, "INFO": 0, "WARNING": 0, "ERROR": 0, "CRITICAL": 0}


class IndentFormatter(logging.Formatter):
    """A logging record formatter with indenting.

    This record formatter tracks an indent depth that is used to nest messages, making
    it easier to track the different sections of validation in printed outputs. It also
    encodes logging levels as single character strings to make logging messages align
    vertically at different depths

    The depth of indenting can be set directly using `FORMATTER.depth = 1`  but it is
    more convenient to use the [push][safedata_validator.logger.IndentFormatter.push]
    and [pop][safedata_validator.logger.IndentFormatter.pop] methods to increase and
    decrease indenting depth.

    The `extra` argument to logger messages can be used to provide a dictionary and is
    used in this subclass to provide the ability to `join` a list of entries as comma
    separated list on to the end of the message.

    Args:
        fmt: The formatting used for emitted LogRecord instances
        datefmt: A date format string
        indent: This string is used for each depth of the indent.
    """

    def __init__(
        self,
        fmt: str = "%(levelcode)s %(message)s",
        datefmt: str | None = None,
        indent: str = "    ",
    ) -> None:
        logging.Formatter.__init__(self, fmt, datefmt)
        self.depth = 0
        self.indent = indent

    def pop(self, n: int = 1) -> None:
        """A convenience method to decrease the indentation of the formatter.

        Args:
            n: Decrease the indentation depth by n.
        """

        self.depth = max(0, self.depth - n)

    def push(self, n: int = 1) -> None:
        """A convenience method to increase the indentation of the formatter.

        Args:
            n: Increase the indentation depth by n.
        """
        self.depth = self.depth + n

    def format(self, rec: logging.LogRecord):
        """Format indented messages with encoded logger message levels.

        Args:
            rec: The logging record to be formatted.
        """

        # Format message
        msg = logging.Formatter.format(self, rec)

        # Add any joined values as repr
        if hasattr(rec, "join"):
            msg += ", ".join([repr(o) for o in getattr(rec, "join")])

        return self.indent * self.depth + msg


# Setup the logging instance
LOGGER = logging.getLogger(__name__)
"""logging.Logger: The safedata_validator Logger instance

This logger instance is used throughout the package for outputting validation
information and is customised to provide counts of error messages and an
indented logging style formatted.
"""


FORMATTER = IndentFormatter()
"""IndentFormatter: The safedata_validator message formatter

This formatter instance is used with the main logging stream handler for the
package and is exposed globally to make it easier to adjust indent depth using
the custom pop and push methods.
"""


def use_file_logging(filename: Path, level: int = logging.DEBUG) -> None:
    """Switch to file logging to a provided file path.

    This function adds a FileCounterHandler to :data:`~safedata_validator.logger.LOGGER`
    using the provided ``filename`` path. It will remove any other existing handlers
    first.

    Args:
        filename: The path to a file to use for logging.
        level: The lowest logging level to be recorded in the file.

    Raises:
        RuntimeError: If the file handler already exists. If the logging is to move to a
            new file, the existing handler needs to be explicitly removed first.
    """

    # Check for an existing file logger
    for handler in LOGGER.handlers:
        if isinstance(handler, FileCounterHandler) and handler.name == "sdv_file_log":
            raise RuntimeError(f"Already logging to file: {handler.baseFilename}")

    # Remove an existing stream logger.
    try:
        sdv_stream_log = next(
            handler for handler in LOGGER.handlers if handler.name == "sdv_stream_log"
        )
    except StopIteration:
        sdv_stream_log = None

    if sdv_stream_log:
        sdv_stream_log.close()
        LOGGER.removeHandler(sdv_stream_log)

    # Add a file handler
    handler = FileCounterHandler(filename=filename)
    handler.setFormatter(FORMATTER)
    handler.name = "sdv_file_log"
    LOGGER.addHandler(handler)
    LOGGER.setLevel(level)


def use_stream_logging(level: int = logging.DEBUG) -> None:
    """Switch to stream logging.

    This function attempts to remove the ``vr_logfile`` FileHandler that is added by
    :func:`~virtual_rainforest.core.logger.add_file_logger`. If that file handler is
    not found it simple exits, otherwise it removes the file handler and restores
    message propagation.
    """

    # Remove an existing file logger.
    try:
        sdv_file_log = next(
            handler for handler in LOGGER.handlers if handler.name == "sdv_file_log"
        )
    except StopIteration:
        sdv_file_log = None

    if sdv_file_log:
        sdv_file_log.close()
        LOGGER.removeHandler(sdv_file_log)

    # Check for an existing stream logger to avoid duplication
    for handler in LOGGER.handlers:
        if (
            isinstance(handler, StreamCounterHandler)
            and handler.name == "sdv_stream_log"
        ):
            return

    # Add a stream handler
    handler = StreamCounterHandler()
    handler.setFormatter(FORMATTER)
    handler.name = "sdv_stream_log"
    LOGGER.addHandler(handler)
    LOGGER.setLevel(level)


# Initialise with stream logging.
use_stream_logging()

#
# CONVENIENCE FUNCTIONS
#


def get_handler():
    """Helper function to get a reference to the current logging handler."""
    return next(hdlr for hdlr in LOGGER.handlers if hdlr.name.startswith("sdv"))


def log_and_raise(
    msg: str, exception: type[Exception], extra: dict | None = None
) -> None:
    """Emit a critical error message and raise an Exception.

    This convenience function adds a critical level message to the logger and
    then raises an exception with the same message. This is intended only for
    use in loading resources: the package cannot run properly with misconfigured
    resources but errors with the data checking should log and carry on.

    Args:
        msg: A message to add to the log
        exception: An exception type to be raised
        extra: A dictionary of extra information to be passed to the logger
    """

    LOGGER.critical(msg, extra=extra)
    raise exception(msg)


# See
# https://stackoverflow.com/questions/10176226/how-do-i-pass-extra-arguments-to-a-python-decorator


def loggerinfo_push_pop(wrapper_message: str) -> Callable:
    """Wrap a callable with an Info logging message and indentation.

    This decorator is used to reduce boilerplate logger code within functions. It emits
    a message and then increases the indentation depth while the wrapped function is
    running.

    Args:
        wrapper_message: The test to use in the info logging message.
    """

    def decorator_func(function: Callable) -> Callable:
        @wraps(function)
        def wrapped_func(*args, **kwargs: Any):
            # Emit the logger info and step in a level
            LOGGER.info(wrapper_message)
            FORMATTER.push()

            # Invoke the wrapped function
            retval = function(*args, **kwargs)

            # Step back out
            FORMATTER.pop()

            return retval

        return wrapped_func

    return decorator_func

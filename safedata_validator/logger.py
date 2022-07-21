""" ## Logging setup for safedata_validator

This submodule extends the standard logging setup to provide extra functionality
and to expose some global logging objects for use throughout the code.

The submodule defines CounterHandler as a subclass of logging.StreamHandler,
extended to maintain a count of messages emitted at the different logging
levels. An instance of this (`COUNTER_HANDLER`) is used with a StringIO
instance (`LOG`) to maintain a logging history and provided error counts.

A second default logging.StreamHandler instance (`CONSOLE_HANDLER`) is used to
emit records to the command line for use in scripts.

The submodule also defines IndentFormatter as a subclass of logging.Formatter,
which is used in both `COUNTER_HANDLER` and `CONSOLE_HANDLER` to provide compact
messages with variable indentation to show different sections of the validation
process.

Note that the handlers are created when the module is loaded, so when running
behind a web server, the content of the handlers persist between runs of the
code. To avoid constant concatenation of outputs, the logger should be cleared
when a new Dataset is being validated.
"""

import logging
from functools import wraps
from io import StringIO


class CounterHandler(logging.StreamHandler):
    def __init__(self, *args, **kwargs):
        """This is an subclass of `logging.StreamHandler` that maintains a count of
        calls at each log level.
        """
        logging.StreamHandler.__init__(self, *args, **kwargs)
        self.counters = {"DEBUG": 0, "INFO": 0, "WARNING": 0, "ERROR": 0, "CRITICAL": 0}

    def emit(self, record: logging.LogRecord):
        """The emit method for CounterHandler emits the message but also increments the
        counter for the message level.

        Args:
            record: A `logging.LogRecord` instance.
        """
        self.counters[record.levelname] += 1

        msg = self.format(record=record)
        self.stream.write(msg)
        self.flush()

    def reset(self):
        """This method is used to reset the counters to zero"""
        self.counters = {"DEBUG": 0, "INFO": 0, "WARNING": 0, "ERROR": 0, "CRITICAL": 0}


class IndentFormatter(logging.Formatter):
    def __init__(
        self,
        fmt: str = "%(indent)s%(levelcode)s %(message)s",
        datefmt: str = None,
        indent: str = "    ",
    ) -> None:
        """A logging record formatter with indenting

        This record formatter tracks an indent depth that is used to nest
        messages, making it easier to track the different sections of validation
        in printed outputs. It also encodes logging levels as single character
        strings to make logging messages align vertically at different depths

        The depth of indenting can be set directly using:

            FORMATTER.depth = 1

        but it is more convenient to use the
        [push][safedata_validator.logger.IndentFormatter.push] and
        [pop][safedata_validator.logger.IndentFormatter.pop] methods to increase
        and decrease indenting depth.

        The `extra` argument to logger messages can be used to provide a
        dictionary and is used in this subclass to provide the ability to `join`
        a list of entries as comma separated list on to the end of the message.

        Args:
            fmt: The formatting used for emitted LogRecord instances
            datefmt: A date format string
            indent: This string is used for each depth of the indent.
        """
        logging.Formatter.__init__(self, fmt, datefmt)
        self.depth = 0
        self.indent = indent

    def pop(self, n: int = 1) -> None:
        """A convenience method to increase the indentation of the formatter.

        Args:
            n: Increase the indentation depth by n.
        """

        self.depth = max(0, self.depth - n)

    def push(self, n: int = 1) -> None:
        """A convenience method to decrease the indentation of the formatter.

        Args:
            n: Decrease the indentation depth by n.
        """
        self.depth = self.depth + n

    def format(self, rec: logging.LogRecord):
        """A custom format method to output indented messages with encoded
        logger message levels.
        """
        rec.indent = self.indent * self.depth

        # encode level
        codes = {
            "DEBUG": ">",
            "INFO": "-",
            "WARNING": "?",
            "ERROR": "!",
            "CRITICAL": "X",
        }
        rec.levelcode = codes[rec.levelname]

        # format message
        msg = logging.Formatter.format(self, rec)

        # add any joined values as repr
        if hasattr(rec, "join"):
            msg += ", ".join([repr(o) for o in rec.join])

        return msg


# Setup the logging instance
LOGGER = logging.getLogger(__name__)
"""logging.Logger: The safedata_validator Logger instance

This logger instance is used throughout the package for outputting validation
information and is customised to provide counts of error messages and an
indented logging style formatted.
"""

LOG = StringIO()
"""io.StringIO: The safedata_validator message log

This StringIO object is attached to a stream handler for LOGGER and is used
to keep a programmatically accessible log of the messages emitted during validation.
It should be truncated between validation runs to remove messages from previous
runs.
"""

COUNTER_HANDLER = CounterHandler(LOG)
"""CounterHandler: The safedata_validator counting handler

This handler is used to track the number of messages emitted by the logger in
different logging levels and emits log messages to the LOG StringIO instance to
keep a record of the messages. It is exposed globally to make it easy to access
counts and the validation log programmatically.
"""

CONSOLE_HANDLER = logging.StreamHandler()
"""logging.StreamHandler: A logger outputting to the console

This handler is used to write logging messages to the console, and is exposed
globally to make it easier to mute it.
"""

FORMATTER = IndentFormatter()
"""IndentFormatter: The safedata_validator message formatter

This formatter instance is used with the main logging stream handler for the
package and is exposed globally to make it easier to adjust indent depth using
the custom pop and push methods.
"""

# Combine those instances into the full LOGGER setup with 2 handlers:
# - COUNTER_HANDLER - logs records into LOG to keep a history of the validation and
#        maintains counts of message levels handler
# - CONSOLE_HANDLER - logs records to the console for command line use and can be
#     muted

COUNTER_HANDLER.setFormatter(FORMATTER)
CONSOLE_HANDLER.setFormatter(FORMATTER)

LOGGER.setLevel(logging.DEBUG)
LOGGER.addHandler(COUNTER_HANDLER)
LOGGER.addHandler(CONSOLE_HANDLER)

#
# CONVENIENCE FUNCTIONS
#


def log_and_raise(msg: str, exception: Exception) -> None:
    """
    This convenience function adds a critical level message to the logger and
    then raises an exception with the same message. This is intended only for
    use in loading resources: the package cannot run properly with misconfigured
    resources but errors with the data checking should log and carry on.

    Args:
        msg: A message to add to the log
        exception: An exception type to be raised

    Returns:
        None
    """

    LOGGER.critical(msg)
    raise exception(msg)


# See
# https://stackoverflow.com/questions/10176226/how-do-i-pass-extra-arguments-to-a-python-decorator


def loggerinfo_push_pop(wrapper_message):
    """
    This decorator is used to reduce boilerplate logger code within functions.
    It emits a message and then increases the indentation depth while the
    wrapped function is running.
    """

    def decorator_func(function):
        @wraps(function)
        def wrapped_func(*args, **kwargs):

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

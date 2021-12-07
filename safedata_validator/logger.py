import logging
from io import StringIO

"""
Logger setup - setup the standard logger to provide:
  1) A handler providing a counter of the number of calls to each log level
  2) A formatter to provide user controlled indentation and to
     code levels as symbols for compactness
  3) A handler writing to a global log file written to a StringIO container
  4) Optionally a handler also writing the log to stdout can be added when
     a Dataset instance is created.

Note that the handlers are created when the module is loaded, so when running
behind a web server, the content of the handlers persist between runs of the 
code. To avoid constant concatenation of outputs, Dataset.__init__() empties 
them.

The depth of indenting can either be set when a message is logged:

    LOGGER.info('Configuring Resources', extra=dict(depth=0))
    
or alternatively set directly to set the indent depth for subsequent messages:
 
    FORMATTER.depth = 1
"""


class CounterHandler(logging.Handler):
    """
    Handler instance that maintains a count of calls at each log level
    """

    def __init__(self, *args, **kwargs):
        logging.Handler.__init__(self, *args, **kwargs)
        self.counters = {'DEBUG': 0, 'INFO': 0, 'WARNING': 0, 'ERROR': 0, 'CRITICAL': 0}

    def emit(self, rec):
        self.counters[rec.levelname] += 1

    def reset(self):
        self.counters = {'DEBUG': 0, 'INFO': 0, 'WARNING': 0, 'ERROR': 0, 'CRITICAL': 0}


class IndentFormatter(logging.Formatter):
    """
    A formatter that provides an indentation depth and encodes the logging
    levels as single character strings. The extra argument to logger messages
    can be used to provide a dictionary to set:
        - 'join': a list of entries to join as comma separated list on
          to the end of the message.
        - 'quote': a flag to set joined entries to be quoted to show
          whitespace around values.
    """

    def __init__(self, fmt=None, datefmt=None):
        logging.Formatter.__init__(self, fmt, datefmt)
        self.depth = 0

    def pop(self, n=1):

        self.depth = max(0, self.depth - n)

    def push(self, n=1):

        self.depth = self.depth + n

    def format(self, rec):

        rec.indent = '    ' * self.depth

        # encode level
        codes = {'DEBUG': '>', 'INFO': '-', 'WARNING': '?', 'ERROR': '!', 'CRITICAL': 'X'}
        rec.levelcode = codes[rec.levelname]

        # format message
        msg = logging.Formatter.format(self, rec)
        # add any joined values as repr
        if hasattr(rec, 'join'):
            msg += ', '.join([repr(o) for o in rec.join])

        return msg


# Setup the logging instance
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)

# Add the counter handler
CH = CounterHandler()
CH.setLevel(logging.DEBUG)
LOGGER.addHandler(CH)

# Create the formatter
FORMATTER = IndentFormatter("%(indent)s%(levelcode)s %(message)s")

# Create a StringIO object to hold the log, set a stream handler to use that,
# attach the custom formatter and add it to the logger. LOG will then contain
# a complete record of logging messages.
LOG = StringIO()
LOG_SH = logging.StreamHandler(LOG)
LOG_SH.setFormatter(FORMATTER)
LOGGER.addHandler(LOG_SH)

# A) REPORT TRACKING - clear out previous content in the logger. When the function runs
# on a webserver, the same logger instance is used for all runs, so the log contents
# persist unless they are tidied up.
CH.reset()
LOG.truncate(0)

verbose = True

if verbose:
    # Look to see if the logger already has a verbose handler writing to the console and
    # add one if it doesn't. This is primarly to provide output for command line usage.
    handler_names = [handler.get_name() for handler in LOGGER.handlers]
    if 'console_log' not in handler_names:
        console_log = logging.StreamHandler()
        console_log.setFormatter(FORMATTER)
        console_log.set_name('console_log')
        LOGGER.addHandler(console_log)
elif not verbose:
    # Remove the console_log handler if one exists
    handler_names = [handler.get_name() for handler in LOGGER.handlers]
    if 'console_log' in handler_names:
        console_log_handler = LOGGER.handlers[handler_names.index('console_log')]
        LOGGER.removeHandler(console_log_handler)

#
# CONVENIENCE FUNCTIONS
#

def log_and_raise(logger, msg, raise_type):
    """ A convenience function that adds a message and a critical entry
    to the logger and then raises an exception with the same message.

    Args:
        logger: A logging.Logger instance
        msg: A message to add to the log and error
        raise_type: An exception type

    Returns:
        None
    """

    logger.critical(msg)
    raise raise_type(msg)


def loggerinfo_push_pop(wrapper_message):
    """This is a convenience decorator to allow reduce boilerplate logger code
    around functions."""

    def decorator_func(func):
        def wrapper_func(*args, **kwargs):

            # Emit the logger info and step in a level
            LOGGER.info(wrapper_message)
            FORMATTER.push()

            # Invoke the wrapped function 
            retval = func(*args, **kwargs)

            # Step back out
            FORMATTER.pop()

            return retval
        return wrapper_func
    return decorator_func

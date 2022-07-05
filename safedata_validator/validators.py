import re
from collections import Counter
from collections.abc import Iterable
from numbers import Number
from typing import Any

from safedata_validator import NA_type
from safedata_validator.logger import LOGGER

"""
Some regular expressions used within validators are defined at the module
level, rather than compiling for each instance
"""

RE_WSPACE_ONLY = re.compile(r"^\s*$")
RE_CONTAINS_WSPACE = re.compile(r"\s")
RE_CONTAINS_PUNC = re.compile(r"[,;:]")

RE_WSPACE_AT_ENDS = re.compile(r"^\s+.+|.+\s+$")
RE_DMS = re.compile(r'[°\'"dms’”]+')

RE_R_ELLIPSIS = re.compile(r"^\\.{2}[0-9]+$|^\\.{3}$")
RE_R_NAME_CHARS = re.compile(r"^[\w\.]+$")
RE_R_NAME_BAD_START = re.compile(r"^_|^\\.[0-9]")

# String values returned by openpyxl for errors when worksheets opened using
# data_only=True
EXCEL_ERRORS = set(
    [
        "#DIV/0!",
        "#NAME?",
        "#N/A",
        "#NUM!",
        "#VALUE!",
        "#REF!",
        "#NULL!",
        "#SPILL!",
        "#CALC!",
    ]
)

# First some simple single value validators


def blank_value(value):
    """Simple check for 'empty' cell values (None, zero length string or whitespace
    only)"""
    return (value is None) or (str(value).isspace()) or not len(str(value))


def valid_r_name(string):
    """Check that a field name is a valid r name"""

    reserved_words = [
        "if",
        "else",
        "repeat",
        "while",
        "function",
        "for",
        "in",
        "next",
        "break",
        "TRUE",
        "FALSE",
        "NULL",
        "Inf",
        "NaN",
        "NA",
        "NA_integer_",
        "NA_real_",
        "NA_complex_",
        "NA_character_",
    ]

    valid = True

    # length (not worrying about pre 2.13.0)
    if len(string) > 10000:
        valid = False

    # ellipsis reserved words ('...' and '..#')
    if RE_R_ELLIPSIS.match(string):
        valid = False

    # is it a reserved word?
    if string in reserved_words:
        valid = False

    if not RE_R_NAME_CHARS.match(string):
        valid = False

    if RE_R_NAME_BAD_START.match(string):
        valid = False

    return valid


# Now a class of validator to handle sets of values and track
# invalid values


class Filter:
    def __init__(self, values, keep_failed=True):
        """Filter an iterable of input values

        This class is a general purpose class to check the values in an iterable.
        The class has two static methods, which can be overloaded in subclasses
        to yield different checking behaviour. These methods are:

            * `tfunc`, which should perform a boolean test on a value and,
            * `rfunc`, which should return a value for failing inputs.

        The __iter__ method is set to iterate over the resulting checked values
        and the __bool__ method is set to report whether any values failed the
        checks. The original failing values themselves are saved in the `failed`
        attribute. If `keep_failed` is True, the output of `rfunc` for the failed
        value is added to the list of output values - acting as a cleaner - but
        if `keep_failed` is False, the value is dropped from the output values -
        acting as a filter.

        Args:
            values: An iterable of values to be filtered
            keep_failed: A boolean showing whether to keep failed values in
                the output iterable after passing through tfunc, or reduce
                the values to only passing values.

        Attributes:
            failed: A list of failing values
            values: An iterable of filtered values
        """

        self.failed = []
        self.keep_failed = keep_failed
        self.values = [v for v in self._filter(values)]

    def _filter(self, values):

        for val in values:
            if self.tfunc(val):
                yield val
            else:
                self.failed.append(val)
                if self.keep_failed:
                    yield self.rfunc(val)

    @staticmethod
    def tfunc(val):
        return True

    @staticmethod
    def rfunc(val):
        return val

    def __bool__(self):
        return not self.failed

    def __repr__(self):
        return str(not self.failed)

    def __iter__(self):
        yield from self.values


class IsString(Filter):
    @staticmethod
    def tfunc(val):
        return isinstance(val, str)

    @staticmethod
    def rfunc(val):
        return val


class IsNumber(Filter):
    @staticmethod
    def tfunc(val):
        return isinstance(val, (float, int))

    @staticmethod
    def rfunc(val):
        return val


class IsNotNumericString(Filter):
    @staticmethod
    def tfunc(val):
        try:
            float(val)
            return False
        except ValueError:
            return True

    @staticmethod
    def rfunc(val):
        return val


class IsLocName(Filter):
    @staticmethod
    def tfunc(val):

        if isinstance(val, (str, int)) or (isinstance(val, float) and val.is_integer()):
            return True

        return False

    @staticmethod
    def rfunc(val):
        return str(val)


class IsNotBlank(Filter):
    @staticmethod
    def tfunc(val):
        return not ((isinstance(val, str) and val.isspace()) or (val is None))

    @staticmethod
    def rfunc(val):
        return val


class IsNotNone(Filter):
    @staticmethod
    def tfunc(val):
        return val is not None

    @staticmethod
    def rfunc(val):
        return val


class IsNotSpace(Filter):

    # Only fails strings that contain whitespace. None is OK

    @staticmethod
    def tfunc(val):
        return not (isinstance(val, str) and val.isspace())

    @staticmethod
    def rfunc(val):
        return None


class IsNotPadded(Filter):
    @staticmethod
    def tfunc(val):
        # Deliberately excludes whitespace cells
        return not (isinstance(val, str) and (val != val.strip()) and not val.isspace())

    @staticmethod
    def rfunc(val):
        return val.strip()


class IsNotNA(Filter):
    @staticmethod
    def tfunc(val):
        return not ((isinstance(val, str) and val == "NA"))

    @staticmethod
    def rfunc(val):
        return val


class IsNotExcelError(Filter):
    @staticmethod
    def tfunc(val):
        return not (isinstance(val, str) and (val in EXCEL_ERRORS))

    @staticmethod
    def rfunc(val):
        return val


class IsLower(Filter):
    @staticmethod
    def tfunc(val):
        return not isinstance(val, str)

    @staticmethod
    def rfunc(val):
        return val.lower()


class NoPunctuation(Filter):
    @staticmethod
    def tfunc(val):
        return isinstance(val, str) and (not RE_CONTAINS_PUNC.search(val))

    @staticmethod
    def rfunc(val):
        return val


# HasDuplicates acts at the level of the set of values, not individual values,
# so is not a Filter subclass


class HasDuplicates:
    def __init__(self, values):
        """Find any duplicates in a set of values

        Args:
            values: An iterable of values to be checked

        Attributes:
            duplicated: A list of duplicated values
        """

        self.duplicated = list(self._get_duplicates(values))

    @staticmethod
    def _get_duplicates(values):
        # https://stackoverflow.com/questions/46554866
        c = Counter()
        seen = set()
        for i in values:
            c[i] += 1
            if c[i] > 1 and i not in seen:
                seen.add(i)
                yield i

    def __bool__(self):
        return bool(len(self.duplicated))

    def __repr__(self):
        return str(bool(len(self.duplicated)))


class IsInSet:
    def __init__(self, values, test_set):
        """Test whether all values are found in test_set.

        Args:
            values: An iterable of values to be filtered
            test_set: Another iterable containing accepted values

        Attributes:
            failed: A list of failing values
            values: An iterable of values found in set
        """

        self.failed = []
        self.test_set = test_set
        self.values = [v for v in self._filter(values)]

    def _filter(self, values):

        for val in values:
            if val in self.test_set:
                yield val
            else:
                self.failed.append(val)

    def __bool__(self):
        return not self.failed

    def __repr__(self):
        return str(not self.failed)

    def __iter__(self):
        yield from self.values


class TypeCheck:
    def __init__(self, values: Iterable, types: tuple):
        """Filter an iterable of input values

        This class is used to check whether all the values are members of a
        given set of types.

        The __iter__ method is set to iterate over the resulting checked values
        and the __bool__ method is set to report whether any values failed the
        checks. The failing values themselves are saved in the `failed`
        attribute.

        Args:
            values: An iterable of values to be filtered
            types: A tuple of accepted types

        Attributes:
            failed: A list of failing values
            values: An iterable of filtered values
        """

        self.failed = []
        self.values = [v for v in self._filter(values, types)]

    def _filter(self, values, types):

        for val in values:
            if isinstance(val, types):
                yield val
            else:
                self.failed.append(val)

    @staticmethod
    def tfunc(val):
        return True

    @staticmethod
    def rfunc(val):
        return val

    def __bool__(self):
        return not self.failed

    def __repr__(self):
        return str(not self.failed)

    def __iter__(self):
        yield from self.values


# INFO - pandas has Excel reading including vectorised string operations?
#        It uses openpyxl for xlsx but doesn't seem to use the memory optimised
#        read_only workbook. Speed vs memory use? Also less flexible
#      - datatable also interesting, but less flexible and (10/2021) relies
#        on xlrd, which is no longer recommended for xslx files.


class GetDataFrame:
    def __init__(self, ws, header_row=1):
        """This class takes a worksheet and extracts a data frame starting at a
        given header row. It then cleans off any blank terminal rows and columns,
        which are sometimes loaded from Excel, and returns a dictionary of columns.

        Args:
            ws: An openpyxl Worksheet instance
            header_row: The row containing the field headers.

        Attributes:
            headers: A list of header values
            data_columns: A list of list containing column data
            bad_headers: A dictionary of malformed header values
        """

        self.headers = []
        self.bad_headers = dict()
        self.data_columns = []

        # Load the rows and pop off the top tuple for headers, converting it to
        # a list so that terminal Nones can be trimmed off
        data = ws.iter_rows(min_row=header_row, values_only=True)
        data_rows = [rw for rw in data]

        # If there are fewer than 2 rows, then there is either no data or only a
        # header.
        if len(data_rows) < 2:
            return

        headers = list(data_rows.pop(0))

        # Strip off terminal empty rows from the data
        # - empty rows in the middle of a dataframe should remain and trigger
        #   an error elsewhere.
        last_is_empty = True

        while last_is_empty:
            last_is_empty = set(data_rows[-1]) == {None}
            if last_is_empty:
                del data_rows[-1]

        # Pivot to columns
        data_columns = [cl for cl in zip(*data_rows)]
        del data_rows

        # Strip off terminal columns where header and values are none
        last_is_empty = True

        while last_is_empty:
            data_is_empty = set(data_columns[-1]) == {None}
            header_is_empty = headers[-1] is None
            last_is_empty = header_is_empty & data_is_empty

            if last_is_empty:
                del data_columns[-1]
                del headers[-1]

        # Check the headers
        # - should not be blank
        headers = IsNotBlank(headers)
        if not headers:
            LOGGER.error("Headers contain empty or whitespace cells")
            self.bad_headers["blank"] = headers.failed

        # - should be strings
        headers = IsString(headers)
        if not headers:
            LOGGER.error(
                "Headers contain non-string values: ", extra={"join": headers.failed}
            )
            self.bad_headers["non_string"] = headers.failed

        # - should have no padding
        headers = IsNotPadded(headers)
        if not headers:
            LOGGER.error(
                "Headers contain whitespace padded strings: ",
                extra={"join": headers.failed},
            )
            self.bad_headers["padded"] = headers.failed

        # - should be unique (after removing whitespace)
        dupes = HasDuplicates(headers)
        if dupes:
            LOGGER.error(
                "Headers contain duplicated values: ", extra={"join": dupes.duplicated}
            )
            self.bad_headers["duplicated"] = dupes.duplicated

        # Do not package headers and columns into a dict - will lose duplicated
        # keys and can still validate fields with dupe keys
        self.headers = headers.values
        self.data_columns = data_columns

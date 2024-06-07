"""The validators module.

The module contains a set of classes providing data validation for use throughout the
safedata_validator package
"""

import re
from collections import Counter
from collections.abc import Iterable
from typing import Any

from openpyxl.worksheet.worksheet import Worksheet

from safedata_validator.logger import LOGGER

RE_WSPACE_ONLY = re.compile(r"^\s*$")
RE_CONTAINS_WSPACE = re.compile(r"\s")
RE_CONTAINS_PUNC = re.compile(r"[,;:]")

RE_WSPACE_AT_ENDS = re.compile(r"^\s+.+|.+\s+$")
RE_DMS = re.compile(r'[°\'"dms’”]+')  # noqa: RUF001

RE_R_ELLIPSIS = re.compile(r"^\\.{2}[0-9]+$|^\\.{3}$")
RE_R_NAME_CHARS = re.compile(r"^[\w\.]+$")
RE_R_NAME_BAD_START = re.compile(r"^_|^\\.[0-9]")
"""Constants defining global regular expressions for the module"""

EXCEL_ERRORS = {
    "#DIV/0!",
    "#NAME?",
    "#N/A",
    "#NUM!",
    "#VALUE!",
    "#REF!",
    "#NULL!",
    "#SPILL!",
    "#CALC!",
}
"""A global set of openpyxl error strings for Worksheets using data_only=True"""


def blank_value(value) -> bool:
    """Check if a single value from a data table is empty.

    Empty values include None, zero length strings and whitespace only strings.
    """
    return (value is None) or (str(value).isspace()) or not len(str(value))


def valid_r_name(string) -> bool:
    """Check that string is a valid R object name."""

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


class Filter:
    """A base class to filter an iterable of input values.

    This base class provides the core functionality to check the values in an iterable.
    The class has two static methods, which can be overloaded in subclasses to yield
    different checking behaviour. These methods are:

        * `tfunc`, which should perform a boolean test on a value and,
        * `rfunc`, which should return a value for failing inputs.

    Values which fail `tfunc` are saved in the `failed` attribute. If `keep_failed` is
    True, the output of `rfunc` for the failed value is added to the list of output
    values - acting as a cleaner - but if `keep_failed` is False, the value is dropped
    from the output values - acting as a filter.

    Two magic methods are set for the class:

    * __iter__  is set to iterate over the resulting filtered values
    * __bool__  reports whether any values failed the checks in `tfunc`

    Args:
        values: An iterable of values to be filtered
        keep_failed: A boolean showing whether to keep failed values in the output
            iterable after passing through tfunc, or reduce the values to only passing
            values.

    Attributes:
        failed: A list of failing values
        values: An iterable of filtered/cleaned values
    """

    def __init__(self, values, keep_failed=True):
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
    def tfunc(val) -> bool:
        """Apply a test to a given value."""
        return True

    @staticmethod
    def rfunc(val) -> Any:
        """Clean a value that has failed testing."""
        return val

    def __bool__(self) -> bool:
        """Filter instances compare True when no values have failed testing."""
        return not self.failed

    def __repr__(self) -> str:
        """Filter instance representation reports the pass or fail status of testing."""
        return str(not self.failed)

    def __iter__(self) -> Any:
        """Iterating over a Filter instance yields the cleaned/filtered values."""
        yield from self.values


class IsString(Filter):
    """A Filter subclass for string values.

    The `tfunc` method overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method to check for string
    values. Failing values are kept unchanged in the instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test for string values."""
        return isinstance(val, str)

    @staticmethod
    def rfunc(val) -> Any:
        """Return failing values unchanged."""
        return val


class IsNumber(Filter):
    """A Filter subclass for numeric values.

    The `tfunc` method overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method to check if values are a
    float or int. Failing values are kept unchanged in the instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test for float or int values."""
        return isinstance(val, float | int)

    @staticmethod
    def rfunc(val) -> Any:
        """Return failing values unchanged."""
        return val


class IsNotNumericString(Filter):
    """A Filter subclass to trap numeric strings.

    The `tfunc` method overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method to return False if the
    value can be cast to a numeric value. Failing values are kept unchanged in the
    instance values.
    """

    @staticmethod
    def tfunc(val):
        """Test that a value does not represent a numeric value."""
        try:
            float(val)
            return False
        except ValueError:
            return True

    @staticmethod
    def rfunc(val):
        """Return failing values unchanged."""
        return val


class IsLocName(Filter):
    """A Filter subclass to check location names.

    Location names are typically strings but can also be numeric codes, which might load
    as integer or float values. The `tfunc` method therefore overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method to return False if the
    value is not a string, an integer or a float that is an integer. Failing values are
    kept unchanged in the instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test that a value is a string or integer."""
        if isinstance(val, str | int) or (isinstance(val, float) and val.is_integer()):
            return True

        return False

    @staticmethod
    def rfunc(val) -> str:
        """Convert failing values to a string."""
        return str(val)


class IsNotBlank(Filter):
    """A Filter subclass to catch blank values.

    Blank values include None, empty strings and whitespace only strings. The `tfunc`
    method therefore overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method to detect blank values.
    Failing values are kept unchanged in the instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test if a value is blank."""
        return not ((isinstance(val, str) and val.isspace()) or (val is None))

    @staticmethod
    def rfunc(val) -> Any:
        """Return failing values unchanged."""
        return val


class IsNotNone(Filter):
    """A Filter subclass to catch None values.

    The `tfunc` method overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method simply to test for None
    values. Failing values are kept unchanged in the instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test for None."""
        return val is not None

    @staticmethod
    def rfunc(val) -> Any:
        """Return failing values unchanged."""
        return val


class IsNotSpace(Filter):
    """A Filter subclass to catch whitespace strings.

    The `tfunc` method overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method simply to test for values
    that only contain whitespace. Note that None values or *empty* strings do not fail.
    Failing values are kept unchanged in the instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test for whitespace strings."""
        return not (isinstance(val, str) and val.isspace())

    @staticmethod
    def rfunc(val) -> None:
        """Replace failing values with None."""
        return None


class IsNotPadded(Filter):
    """A Filter subclass to catch padded strings.

    The `tfunc` method overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method simply to test for values
    that include whitespace padding. Note that *empty* strings and whitespace only
    strings are not considered to be padded. Failing values are kept unchanged in the
    instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test for padded strings."""
        # Deliberately excludes whitespace cells
        return not (isinstance(val, str) and (val != val.strip()) and not val.isspace())

    @staticmethod
    def rfunc(val) -> str:
        """Replace padded strings with padding removed."""
        return val.strip()


class IsNotNA(Filter):
    """A Filter subclass to catch NA values.

    The safedata_validator package uses the string 'NA' to explicitly identify missing
    data. The `tfunc` method overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method simply to test for 'NA'
    strings in data values. Failing values are kept unchanged in the instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test for NA strings."""
        return not (isinstance(val, str) and val == "NA")

    @staticmethod
    def rfunc(val) -> Any:
        """Return failing values unchanged."""
        return val


class IsNotExcelError(Filter):
    """A Filter subclass to catch Excel error strings.

    The safedata_validator package loads data from Excel spreadsheets using the
    data_only=True option, which renders Excel errors as strings. The `tfunc` method
    overrides the base [tfunc][safedata_validator.validators.Filter.tfunc] method to
    catch those strings in the list of values. Failing values are kept unchanged in the
    instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test for Excel error strings."""
        return not (isinstance(val, str) and (val in EXCEL_ERRORS))

    @staticmethod
    def rfunc(val) -> Any:
        """Return failing values unchanged."""
        return val


class IsLower(Filter):
    """A Filter subclass to filter values to lowercase strings.

    The `tfunc` method overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method to catch string values
    and convert them to lowercase.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test for string values."""
        return not isinstance(val, str)

    @staticmethod
    def rfunc(val) -> str:
        """Return lowercase values for strings."""
        return val.lower()


class NoPunctuation(Filter):
    """A Filter subclass to catch values containing punctuation.

    The `tfunc` method overrides the base
    [tfunc][safedata_validator.validators.Filter.tfunc] method to catch strings
    containing punctuation. Failing values are kept unchanged in the instance values.
    """

    @staticmethod
    def tfunc(val) -> bool:
        """Test for punctuation in strings."""
        return isinstance(val, str) and (not RE_CONTAINS_PUNC.search(val))

    @staticmethod
    def rfunc(val):
        """Return failing values unchanged."""
        return val


class HasDuplicates:
    """Find any duplicates in a set of values.

    The class checks a provided iterable of valies and identifies whether any values are
    duplicated. The __bool__ magic method returns False if no duplicates are found.

    Args:
        values: An iterable of values to be checked

    Attributes:
        duplicated: A list of duplicated values
    """

    def __init__(self, values: Iterable):
        self.duplicated = list(self._get_duplicates(values))

    @staticmethod
    def _get_duplicates(values) -> Any:
        # https://stackoverflow.com/questions/46554866
        c: Counter = Counter()
        seen = set()
        for i in values:
            c[i] += 1
            if c[i] > 1 and i not in seen:
                seen.add(i)
                yield i

    def __bool__(self) -> bool:
        """Report if any values are duplicated."""
        return bool(len(self.duplicated))

    def __repr__(self) -> str:
        """Class representation shows if values are duplicated."""
        return str(bool(len(self.duplicated)))


class IsInSet:
    """Test values are within a test set.

    A base class that tests if all values in an iterable against a tuple of test
    values. The base class simply tests if the iterable values are found in the test
    values.

    The `filter` method can be overridden to alter the test behaviour. The filter method
    should return passing values and append failing values to the `failed` attribute.

    Args:
        values: An iterable of values to be filtered
        test_values: A tuple containing accepted values

    Attributes:
        failed: A list of failing values not found in the test set.
        values: An iterable of values found in the test set
        test_values: The iterable of test values.
    """

    def __init__(self, values: Iterable, test_values: tuple):
        self.failed: list = []
        self.test_values = test_values
        self.values = [v for v in self.filter(values)]

    def filter(self, values) -> Any:
        """Filter the input values.

        This function returns values found in the test set but appends values not found
        in the set to the self.failed attribute.
        """
        for val in values:
            if val in self.test_values:
                yield val
            else:
                self.failed.append(val)

    def __bool__(self) -> bool:
        """Class compares as False when all values in test set."""
        return not self.failed

    def __repr__(self) -> str:
        """Class represent as 'False' when all values in test set."""
        return str(not self.failed)

    def __iter__(self) -> Any:
        """Class iterates over values found in set."""
        yield from self.values


class TypeCheck(IsInSet):
    """Check the types of a set of values.

    This class is used to check whether all the values are members of a given set of
    types. It overrides the base [filter][safedata_validator.validators.IsInSet.filter]
    method to check that the types of `values` are all instances of test_values.
    """

    def filter(self, values) -> Any:
        """Test the types of values."""
        for val in values:
            if isinstance(val, self.test_values):
                yield val
            else:
                self.failed.append(val)


# INFO - pandas has Excel reading including vectorised string operations?
#        It uses openpyxl for xlsx but doesn't seem to use the memory optimised
#        read_only workbook. Speed vs memory use? Also less flexible
#      - datatable also interesting, but less flexible and (10/2021) relies
#        on xlrd, which is no longer recommended for xslx files.


class GetDataFrame:
    """Get a data frame from an Excel worksheet.

    This class takes a worksheet and extracts a data frame starting at a
    given header row. It then cleans off any blank terminal rows and columns,
    which are sometimes loaded from Excel, and returns a dictionary of columns.

    Args:
        ws: An openpyxl Worksheet instance
        header_row: The row containing the field headers.

    Attributes:
        headers: A list of header values
        data_columns: A list of tuples containing column data
        bad_headers: A dictionary of malformed header values
    """

    def __init__(self, ws: Worksheet, header_row: int = 1):
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
        headers_chk: Filter = IsNotBlank(headers)
        if not headers_chk:
            LOGGER.error("Headers contain empty or whitespace cells")
            self.bad_headers["blank"] = headers_chk.failed

        # - should be strings
        headers_chk = IsString(headers_chk)
        if not headers_chk:
            LOGGER.error(
                "Headers contain non-string values: ",
                extra={"join": headers_chk.failed},
            )
            self.bad_headers["non_string"] = headers_chk.failed

        # - should have no padding
        headers_chk = IsNotPadded(headers_chk)
        if not headers_chk:
            LOGGER.error(
                "Headers contain whitespace padded strings: ",
                extra={"join": headers_chk.failed},
            )
            self.bad_headers["padded"] = headers_chk.failed

        # - should be unique (after removing whitespace)
        dupes = HasDuplicates(headers_chk)
        if dupes:
            LOGGER.error(
                "Headers contain duplicated values: ", extra={"join": dupes.duplicated}
            )
            self.bad_headers["duplicated"] = dupes.duplicated

        # Do not package headers and columns into a dict - will lose duplicated
        # keys and can still validate fields with dupe keys
        self.headers = headers_chk.values
        self.data_columns = data_columns

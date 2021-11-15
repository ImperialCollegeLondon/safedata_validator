import re
from numbers import Number
from collections import Counter
import openpyxl
from typing import Any

"""
Some regular expressions used within validators are defined at the module 
level, rather than compiling for each instance
"""

RE_ORCID = re.compile(r'[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]')
RE_EMAIL = re.compile(r'\S+@\S+\.\S+')
RE_NAME = re.compile(r'[^,]+,[ ]?[^,]+')
RE_WSPACE_ONLY = re.compile(r'^\s*$')
RE_CONTAINS_WSPACE = re.compile(r'\s')
RE_CONTAINS_PUNC = re.compile(r'[,;:]')
RE_DOI = re.compile('https?://(dx.)?doi.org/')
RE_WSPACE_AT_ENDS = re.compile(r'^\s+.+|.+\s+$')
RE_DMS = re.compile(r'[°\'"dms’”]+')


def _screen_values(obj):
    """
    Both Standardiser and Validator objects test their inputs to check
    the types of all values and optionally assess if any are None. This
    shared routine runs this testing and raises an Error if any check fails.

    :param obj: An instance of either a Standardiser or Validator
    :return: None
    """

    if not isinstance(obj.values, (list, tuple)):
        obj.values = [obj.values]

    if not obj.none_ok and any([v is None for v in obj.values]):
        raise TypeError(f'Values contain None')

    if obj.none_ok:
        right_type = [v is None or isinstance(v, obj.valid_types) for v in obj.values]
        if not all(right_type):
            raise TypeError(f'Values not all in types {obj.valid_types} or None')
    else:
        right_type = [isinstance(v, obj.valid_types) for v in obj.values]
        if not all(right_type):
            raise TypeError(f'Values not all in types {obj.valid_types}')


class Standardiser:

    """The base Standardiser class.

    This is initialised with a value or list of values and applies some kind of
    standardisation to each member of that list. To make the class easier to use,
    the __iter__ method is defined to return an iterator over the standardised
    values.

    The conversion applied is defined in the standardise static method, which can
    be overridden to create different Standardisers. The read only property
    valid_types can also be overridden to modify what types of value are accepted.
    Note that the standardise method must be able to handle None values, since these
    are optionally permitted.
    """
    def __init__(self, values, none_ok=False):

        # Look for incoming Standardisers to allow them to be chained
        if issubclass(type(values), Standardiser):
            values = values.values

        self.none_ok = none_ok
        self.values = values
        _screen_values(self)
        standard_values = [self.standardise(v) for v in values]

        if len(values) != len(standard_values):
            raise RuntimeError('Standardiser is not generating a 1 to 1 conversion.')
        else:
            self.values = standard_values

    def __repr__(self):
        return repr(self.values)

    def __iter__(self):
        return iter(self.values)

    @property
    def valid_types(self):
        return int

    @staticmethod
    def standardise(value):
        return None if value is None else value


class ToLower(Standardiser):

    @property
    def valid_types(self):
        return str

    @staticmethod
    def standardise(value):
        return None if (value is None or RE_WSPACE_ONLY.match(value) is not None) else value.lower()


class CellToValue(Standardiser):

    """Converts openpyxl cells to their values
    """

    @property
    def valid_types(self):
        return (openpyxl.cell.read_only.ReadOnlyCell,
                openpyxl.cell.read_only.EmptyCell)

    @staticmethod
    def standardise(value):
        return None if isinstance(value, openpyxl.cell.read_only.EmptyCell) else value.value


class BlankToNone(Standardiser):

    """Converts whitespace cells to None, accepting pretty much any basic value type.
    """

    @property
    def valid_types(self):
        return str, int, float

    @staticmethod
    def standardise(value):
        return None if (value is None or (isinstance(value, str) and
                                          RE_WSPACE_ONLY.match(value) is not None)) else value


class Validator:

    """The base Validator class.

    This is initialised with a value or list of values and applies a test to that
    value or each member of the list. The class instance then populates lists of
    valid and invalid values and the all_valid property indicates whether there
    are any invalid values. To make instances easy to use in logical tests, the
    __bool__ method is defined to return the content of all_valid.

    The test applied is defined in the test static method, which can be overridden
    to create Validators with different criteria. The read only property valid_types
    can also be overridden to modify what types of value are accepted. Note that
    the test method must be able to handle None values, since these are optionally
    permitted.
    """

    def __init__(self, values, none_ok=False):
        """

        :param values:
        :param none_ok:
        """
        self.valid = []
        self.invalid = []
        self.none_ok = none_ok
        self.all_valid = True

        # Look for incoming Standardisers to allow them to be chained
        if issubclass(type(values), Standardiser):
            values = values.values

        self.values = values

        _screen_values(self)

        test_results = self.test(self.values)

        for value, is_valid in test_results:
            if is_valid:
                self.valid.append(value)
            else:
                self.invalid.append(value)

        if self.invalid:
            self.all_valid = False

    def __bool__(self):
        return self.all_valid

    def __repr__(self):
        return str(self.all_valid)

    @property
    def valid_types(self):
        return int

    @staticmethod
    def test(value):
        return value, value > 0

# Implemented


class AllNone(Validator):

    def __init__(self, values, none_ok=True):
        """It makes no use this Validator with none_ok = False, so the defaults
        are changed here and none_ok = False raises an error.
        """
        if not none_ok:
            raise RuntimeError('AllNone should not be called with none_ok=False')
        super(AllNone, self).__init__(values=values, none_ok=none_ok)

    @property
    def valid_types(self):
        return object  # Lets _anything_ through

    @staticmethod
    def test(values):
        return [(v, v is None) for v in values]


class IsBlank(Validator):

    @property
    def valid_types(self):
        return str

    @staticmethod
    def test(values):
        return [(v, (v is None) or (RE_WSPACE_ONLY.match(v) is not None)) for v in values]

class IsNotBlank(Validator):

    @property
    def valid_types(self):
        return str

    @staticmethod
    def test(values):
        return [(v, (v is not None) or (RE_WSPACE_ONLY.match(v) is None)) for v in values]


class IsUnique(Validator):

    @property
    def valid_types(self):
        return str

    @staticmethod
    def test(values):
        return [(c[0], c[1] == 1) for c in Counter(values).items()]


class IsOrcID(Validator):

    @property
    def valid_types(self):
        return str

    @staticmethod
    def test(values):
        return [(v, (v is None or RE_ORCID.match(v))) for v in values]


class IsNotPadded(Validator):

    @property
    def valid_types(self):
        return str, float, int

    @staticmethod
    def test(values):
        return [(v, not isinstance(v, str) or RE_WSPACE_AT_ENDS.match(v) is None) for v in values]


class IsPadded(Validator):

    @property
    def valid_types(self):
        return str, float, int

    @staticmethod
    def test(values):
        return [(v, v is not None and isinstance(v, str) and RE_WSPACE_AT_ENDS.match(v) is not None) for v in values]



# In Progress
#
#
# class IsNotBlankString(Validator):
#
#     @staticmethod
#     def test(values):
#         not_blank = [not(is_blank_string(v)) for v in values]
#         return zip(values, not_blank)
#
#
# class IsNumeric(Validator):
#
#     @staticmethod
#     def test(values):
#
#         return zip(values, [isinstance(v, Number) for v in values])
#
#
# class IsNumericString(Validator):
#
#     @staticmethod
#     def floatable(value):
#         try:
#             float(value)
#             return True
#         except ValueError:
#             return False
#
#     def test(self, values):
#
#         return zip(values, [self.floatable(v) for v in values])
#

#
# def is_blank_string(value):
#
#     return (value is None) or (RE_WSPACE_ONLY.match(value) is not None)
#
#
# def is_blank(value):
#
#     return (value is None) or (isinstance(value, str) and RE_WSPACE_ONLY.match(value) is not None)
#
#
# def blank_to_none(value):
#
#     return None if is_blank(value) else value

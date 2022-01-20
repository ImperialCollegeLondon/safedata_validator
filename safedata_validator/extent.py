from safedata_validator.logger import LOGGER

class Extent:

    def __init__(self, label, datatype, hard_bounds=None, soft_bounds=None):
        """
        This class holds the extent of a particular variable across a dataset.
        It is designed to track the extents required by GEMINI 2: latitude,
        longitude and date, but the implementation is general. Values are fed
        to the class using the update method, which adjusts the extent as
        necessary.

        The Extent instance is created by providing a datatype and optionally
        any hard and soft bounds to be applied. Hard bounds log an error when
        exceeded and soft bounds log a warning.

        Args:
            label: A label for the extent, used in reporting
            datatype: A type or tuple of types for input checking
            hard_bounds: A 2 tuple of hard bounds
            soft_bounds: A 2 tuple of soft bounds
        """

        # The extent is stored internally as a list for ease of update
        # but only accessible via the property as a tuple to avoid it
        # being modifiable by reference. All other variables are similarly
        # protected against adjustment.

        self.label = label
        self._datatype = datatype
        self._extent = [None, None]
        self._hard_bounds = None
        self._soft_bounds = None

        if hard_bounds is not None:
            self._check_bounds(hard_bounds)

        if soft_bounds is not None:
            self._check_bounds(soft_bounds)

        if ((soft_bounds is not None and hard_bounds is not None) and
                (soft_bounds[0] <= hard_bounds[0] or soft_bounds[1] >= hard_bounds[1])):
            raise AttributeError('Hard bounds must lie outside soft bounds')

        self._hard_bounds = hard_bounds
        self._soft_bounds = soft_bounds

    def __repr__(self):
        return f'Extent: {self.label} {self.extent}'

    @property
    def datatype(self):
        return self._datatype

    @property
    def extent(self):
        return tuple(self._extent)

    @property
    def hard_bounds(self):
        return self._hard_bounds

    @property
    def soft_bounds(self):
        return self._soft_bounds

    def _check_bounds(self, bounds):
        """
        Private function to validate hard and soft bounds, these are set at
        dataset initialisation and so raise an error, rather than logging.

        Args:
            bounds: Expecting an iterable of length 2 with low, high values
        """

        values_match_datatype = [isinstance(v, self.datatype) for v in bounds]

        if not all(values_match_datatype):
            raise AttributeError(f'Bounds are not all of type {self.datatype}')

        if len(bounds) != 2 or bounds[1] <= bounds[0]:
            raise AttributeError('Bounds must be provided as (low, high) tuples')

    def update(self, values):
        """
        Takes a set of values and updates the extent of values provided to the
        instance.

        Args:
             values: An iterable of values, which should all be of the datatype(s)
                specified when creating the Extent instance.
        """

        values_match_datatype = [isinstance(v, self.datatype) for v in values]

        if not all(values_match_datatype):
            invalid = [v for v, b in zip(values, values_match_datatype) if not b]
            values = [v for v, b in zip(values, values_match_datatype) if b]
            LOGGER.error(f'Values are not all of type {self.datatype}: ',
                         extra={'join': invalid, 'quote': True})

        values = (min(values), max(values))

        if self.hard_bounds and (self.hard_bounds[0] > values[0] or
                                 self.hard_bounds[1] < values[1]):
            err = f'Values range {values} exceeds hard bounds {self.hard_bounds}'
            LOGGER.error(err)

        if self.soft_bounds and (self.soft_bounds[0] > values[0] or
                                 self.soft_bounds[1] < values[1]):
            wrn = f'Values range {values} exceeds soft bounds {self.soft_bounds}'
            LOGGER.warning(wrn)

        if self.extent[0] is None or values[0] < self.extent[0]:
            self._extent[0] = values[0]

        if self.extent[1] is None or values[1] > self.extent[1]:
            self._extent[1] = values[1]

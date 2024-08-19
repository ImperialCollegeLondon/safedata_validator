"""The extent submodule defines the Extent class to track the extent of a
particular variable across a dataset. It is designed to track the extents
required by GEMINI 2: latitude, longitude and date, but the implementation is
general. Values are fed to a class instance using the
[update][safedata_validator.extent.Extent.update] method, which adjusts the
extent as necessary.

Typical usage:

    ```python
    ext = Extent('latitude', (int, float), hard_bounds=(-90, 90))
    ext.update([1,2,3,4,5,6])
    ```
"""  # noqa D415

from collections.abc import Iterable

from safedata_validator.logger import LOGGER, log_and_raise
from safedata_validator.validators import TypeCheck


class Extent:
    """Track the extent of data.

    An Extent instance is created by providing a datatype and optionally
    any hard and soft bounds to be applied. When an Extent instance is updated,
    values outside hard bounds will generate an error in logging and values
    outside soft bounds will log a warning.

    Args:
        label: A label for the extent, used in reporting
        datatype: A type or tuple of types for input checking
        hard_bounds: A 2 tuple of hard bounds
        soft_bounds: A 2 tuple of soft bounds
    """

    def __init__(
        self,
        label: str,
        datatype: tuple[type, ...],
        hard_bounds: tuple | None = None,
        soft_bounds: tuple | None = None,
    ):
        # The extent is stored internally as a list for ease of update
        # but only accessible via the property as a tuple to avoid it
        # being modifiable by reference. All other variables are similarly
        # protected against adjustment.

        self.label = label
        self._datatype = datatype
        self._extent = [None, None]
        self._hard_bounds = None
        self._soft_bounds = None
        self._populated = False

        if hard_bounds is not None:
            self._check_bounds(hard_bounds)

        if soft_bounds is not None:
            self._check_bounds(soft_bounds)

        if (soft_bounds is not None and hard_bounds is not None) and (
            soft_bounds[0] <= hard_bounds[0] or soft_bounds[1] >= hard_bounds[1]
        ):
            log_and_raise(
                f"Hard bounds must lie outside soft bounds in {label}", AttributeError
            )

        self._hard_bounds = hard_bounds
        self._soft_bounds = soft_bounds

    def __repr__(self):
        """Provide a simple representation of the class."""
        return f"Extent: {self.label} {self.extent}"

    @property
    def datatype(self) -> tuple[type, ...]:
        """Returns the data types accepted by the Extent object."""
        return self._datatype

    @property
    def extent(self) -> tuple:
        """Returns a tuple showing the current extent."""
        return tuple(self._extent)

    @property
    def hard_bounds(self) -> tuple | None:
        """Returns a tuple showing the hard bounds of the Extent object."""
        return self._hard_bounds

    @property
    def soft_bounds(self) -> tuple | None:
        """Returns a tuple showing the hard bounds of the Extent object."""
        return self._soft_bounds

    @property
    def populated(self) -> bool:
        """Returns a boolean showing if the extent has been populated."""
        return self._populated

    def _check_bounds(self, bounds: tuple):
        """Private function to validate hard and soft bounds.

        These are set at dataset initialisation and so raise an error, rather than
        logging.

        Args:
            bounds: Expecting an iterable of length 2 with low, high values
        """

        valid_types = TypeCheck(bounds, self.datatype)

        if not valid_types:
            log_and_raise(f"Bounds are not all of type {self.datatype}", AttributeError)

        if len(bounds) != 2 or bounds[1] <= bounds[0]:
            log_and_raise(
                "Bounds must be provided as (low, high) tuples", AttributeError
            )

    def update(self, values: Iterable) -> None:
        """Update extent of instance based on values contained in an iterable.

        Args:
            values: An iterable of values, which should all be of the
                datatype(s) specified when creating the Extent instance.
        """

        valid_types = TypeCheck(values, self.datatype)

        if not valid_types:
            LOGGER.error(
                f"Values are not all of type {self.datatype}: ",
                extra={"join": valid_types.failed},
            )

        if len(valid_types.values) == 0:
            LOGGER.error("No valid data in extent update")
            return

        minv = min(valid_types.values)
        maxv = max(valid_types.values)

        if self.hard_bounds and (
            self.hard_bounds[0] > minv or self.hard_bounds[1] < maxv
        ):
            LOGGER.error(
                f"Values (range {minv}, {maxv}) exceeds hard bounds {self.hard_bounds}"
            )

        elif self.soft_bounds and (
            self.soft_bounds[0] > minv or self.soft_bounds[1] < maxv
        ):
            LOGGER.warning(
                f"Values (range {minv}, {maxv}) exceeds soft bounds {self.soft_bounds}"
            )

        # Update the bounds, handling None from __init__
        self._extent[0] = min(minv, self._extent[0]) if self._extent[0] else minv
        self._extent[1] = max(maxv, self._extent[1]) if self._extent[1] else maxv
        self._populated = True

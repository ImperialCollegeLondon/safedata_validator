"""Load and check validation resources

The `safedata_validator` package needs access to some local resources and
configuration to work. The three main resources for file validation are:

-   locations: A path to a locations data file, providing a gazetteer of known
    locations and their details.

-   gbif_database: Optionally, a path to a local SQLite copy of the GBIF
    backbone database. If this is not provided, then the slower GBIF online API
    will be used instead.

-   ncbi_database: Optionally, a path to a local SQLite copy of the NCBI
    database files. If this is not provided, then the slower BioEntrez package
    will be used to access the online APU instead.

The [Resources][safedata_validator.resources.Resources] class is used to locate
and validate these resources, and then provide those validated resources to
other components of the package.

A configuration file can be passed as `cfg_path` when creating an instance,
but if no arguments are provided then an attempt is made to find and load
configuration files in the user and then site config locations defined by the
`appdirs` package. See [usage/usage][here] for details.

"""


import os
import sqlite3
import appdirs
import simplejson
from dateutil.parser import isoparse
from validate import Validator, VdtParamError, VdtValueError, is_list
from configobj import ConfigObj, flatten_errors
from dotmap import DotMap
from typing import Union

from safedata_validator.logger import LOGGER, FORMATTER, log_and_raise, loggerinfo_push_pop


CONFIGSPEC = {
    'locations':  'string()',
    'gbif_database':  'string(default=None)',
    'ncbi_database':  'string(default=None)',
    'extents': {
        'temporal_soft_extent':  'date_list(min=2, max=2, default=None)',
        'temporal_hard_extent':  'date_list(min=2, max=2, default=None)',
        'latitudinal_hard_extent':  'float_list(min=2, max=2, default=list(-90, 90))',
        'latitudinal_soft_extent':  'float_list(min=2, max=2, default=None)',
        'longitudinal_hard_extent':  'float_list(min=2, max=2, default=list(-180, 180))',
        'longitudinal_soft_extent':  'float_list(min=2, max=2, default=None)'},
    'zenodo': {
        'community_name':  'string(default=safe)',
        'use_sandbox':  'boolean(default=None)',
        'zenodo_sandbox_api':  'string(default=None)',
        'zenodo_sandbox_token':  'string(default=None)',
        'zenodo_api':  'string(default=None)',
        'zenodo_token':  'string(default=None)'}}
"""dict: The safedata_validator package use the `configobj.ConfigObj`
package to handle the loading and initial validation of resource configuration.
This dict defines the basic expected specification for the configuration and
allows the ConfigObj.validate() method to do basic validation and type
conversions.
"""

def date_list(value, min, max):
    """A configobj.Validator extension function to check configuration values
    containing a list of ISO formatted date strings and to return parsed
    values.

    Args:
        value: A string containing comma-separated ISO date strings
        min: The minimum allowed number of entries
        max: The maximum allowed number of entries
    """
    # min and max are supplied as a string, test conversion to int
    try:
        min = int(min)
    except ValueError:
        raise VdtParamError('min', min)
    try:
        max = int(max)
    except ValueError:
        raise VdtParamError('max', max)

    # Check the supplied value is a list, triggering any issues
    # with list formatting
    value = is_list(value, min=min, max=max)

    # Next, check every member in the list is an ISO date string
    # noting that this strips out time information
    out = []
    for entry in value:

        try:
            parsed_entry = isoparse(entry).date()
        except ValueError:
            raise VdtValueError(entry)

        out.append(parsed_entry)

    # Return parse values
    return out


@loggerinfo_push_pop('Configuring Resources')
class Resources:

    def __init__(self, config: Union[str, list, dict] = None) -> None:
        """Load and check validation resources

        Creating an instance of this class locates and validate resources for
        using the `safedata_validator` package, either from the provided
        configuration details or from the user and then site config locations
        defined by the appdirs package.

        Args:
            config: A path to a configuration file, or a dict or list
                providing package configuration details. The list format
                should provide a list of strings, each representing a
                line in the configuration file. The dict format is a
                dictionary with the required nested dictionary structure
                and values

        Attributes:
            config_type: The method used to specify the resources. One of
                'init_dict', 'init_list', 'init_file', 'user_config' or 'site_config'.
            locations: The path to the locations file
            gbif_database: The path to the GBIF database file or None
            use_local_gbif: Is a local file used or should the GBIF API be used
            ncbi_database: The path to the NCBI database file or None
            use_local_ncbi: Is a local file used or should the NCBI API be used
            valid_locations: The locations defined in the locations file
            location_aliases: Location aliases defined in the locations file
            extents: A DotMap of extent data
            zenodo: A DotMap of Zenodo information
        """

        # User and site config paths
        user_cfg_file = os.path.join(appdirs.user_config_dir(),
                                'safedata_validator', 'safedata_validator.cfg')
        site_cfg_file = os.path.join(appdirs.site_config_dir(),
                                'safedata_validator', 'safedata_validator.cfg')

        # First try and populate from a config file.
        if config is not None:
            if isinstance(config, str):
                config_type = 'init file'
            elif isinstance(config, list):
                config_type = 'init list'
            elif isinstance(config, dict):
                config_type = 'init dict'
        elif os.path.exists(user_cfg_file) and os.path.isfile(user_cfg_file):
            config = user_cfg_file
            config_type = 'user file'
        elif os.path.exists(site_cfg_file) and os.path.isfile(site_cfg_file):
            config = site_cfg_file
            config_type = 'site file'
        else:
            LOGGER.critical(f'No user config in {user_cfg_file}')
            LOGGER.critical(f'No site config in {site_cfg_file}')
            log_and_raise(f'No config files provided or found', RuntimeError)

        # Report resource config location and type
        msg = f'Configuring resources from {config_type}'
        if 'file' in config_type:
            msg += f": {config}"
        LOGGER.info(msg)

        # Try and load the found configuration
        config = self._load_config(config, config_type)

        # Set attributes
        self.locations = config.locations
        self.gbif_database = config.gbif_database
        self.ncbi_database = config.ncbi_database
        self.config_type = config.config_type
        self.config_source = config.config_source

        self.extents = config.extents
        self.zenodo = config.zenodo

        self.use_local_gbif = None
        self.use_local_ncbi = None
        self.valid_locations = None
        self.location_aliases = None

        # Validate the resources
        self._validate_locations()
        self._validate_gbif()
        self._validate_ncbi()

    @staticmethod
    def _load_config(config: Union[str, list, dict], cfg_type: str):
        """Load a configuration file

        This private static method attempts to load a JSON configuration file
        from a path.

        Args:
            config: Passed from Resources.__init__()
            cfg_type: Identifies the route used to provide the configuration details

        Returns:
             If the file does not exist, the function returns None. Otherwise,
             it returns a DotMap of config parameters.
        """

        # Otherwise, there is a file, so try and use it and now raise if there
        # is a problem: don't skip over broken resource configurations.
        # - First, create a validator instance that handles lists of dates
        cf_validator = Validator({'date_list': date_list})

        # - Now load the config input and then apply the basic validation - are
        #   the values of the right type, right count etc.
        config_obj = ConfigObj(config, configspec=CONFIGSPEC)
        valid = config_obj.validate(cf_validator, preserve_errors=True)

        # If there are config file issues, then bail out.
        if  isinstance(valid, dict):
            print(config)
            print(config_obj)
            # print(open(config))
            LOGGER.critical("Configuration issues: ")
            FORMATTER.push()
            for sec, key, err in flatten_errors(config_obj, valid):
                sec.append(key)
                LOGGER.critical(f"In config '{'.'.join(sec)}': {err}")
            FORMATTER.pop()
            raise RuntimeError('Configuration failure')

        # convert to a DotMap for ease
        config_obj = DotMap(config_obj)

        return config_obj

    def _validate_locations(self):
        """Validate and load a locations file

        This private function checks whether a locations path: exists, is a
        JSON file, and contains location and alias data. It populates the
        instance attributes

        Returns:
            None - updates instance.
        """

        if self.locations is None or self.locations == '':
            log_and_raise(f'Locations file missing in configuration',
                          RuntimeError)

        LOGGER.info(f'Validating locations: {self.locations}')

        # Now check to see whether the locations file behaves as expected
        if not os.path.exists(self.locations) and not os.path.isfile(self.locations):
            log_and_raise('Local locations file not found',
                          OSError)

        try:
            loc_payload = simplejson.load(open(self.locations, mode='r'))
        except (simplejson.errors.JSONDecodeError, UnicodeDecodeError):
            log_and_raise('Local locations file not JSON encoded.',
                          OSError)

        # process the locations payload
        if {'locations', 'aliases'} != loc_payload.keys():
            log_and_raise('Locations data malformed',
                          RuntimeError)

        self.valid_locations = loc_payload['locations']
        self.location_aliases = loc_payload['aliases']

    def _validate_gbif(self):
        """Validate the GBIF settings

        This private function checks whether to use the online API or a local
        backbone database and then validates the provided sqlite3 database file.

        Returns:
            None - updates instance.
        """
        if self.gbif_database is None or self.gbif_database == '':
            LOGGER.info('Using GBIF online API to validate taxonomy')
            self.use_local_gbif = False
        else:
            LOGGER.info(f'Validating local GBIF database: {self.gbif_database}')

            # Does the provided path exist and is it a functional SQLite database
            # with a backbone table? Because sqlite3 can 'connect' to any path,
            # use a query attempt to reveal exceptions

            if not os.path.exists(self.gbif_database):
                log_and_raise('Local GBIF database not found',
                              OSError)

            try:
                conn = sqlite3.connect(self.gbif_database)
                _ = conn.execute('select count(*) from backbone;')
            except sqlite3.OperationalError:
                log_and_raise('Local GBIF database does not contain the backbone table',
                              RuntimeError)
            except sqlite3.DatabaseError:
                log_and_raise('Local SQLite database not valid',
                              OSError)
            else:
                self.use_local_gbif = True
            finally:
                conn.close()

    def _validate_ncbi(self):
        """Validate the NCBI settings

        This private function checks whether to use the online API or a local
        backbone database and then validates the provided sqlite3 database files.

        Returns:
            None - updates instance.
        """
        if self.ncbi_database is None or self.ncbi_database == '':
            LOGGER.info('Using NCBI online API to validate taxonomy')
            self.use_local_ncbi = False
        else:
            LOGGER.info(f'Validating local NCBI database: {self.ncbi_database}')

            # Does the provided path exist and is it a functional SQLite database
            # with a backbone table? Because sqlite3 can 'connect' to any path,
            # use a query attempt to reveal exceptions

            if not os.path.exists(self.ncbi_database):
                log_and_raise('Local NCBI database not found',
                              OSError)

            try:
                conn = sqlite3.connect(self.ncbi_database)
                _ = conn.execute('select count(*) from nodes;')
                _ = conn.execute('select count(*) from names;')
                _ = conn.execute('select count(*) from merge;')
            except sqlite3.OperationalError:
                log_and_raise('Local NCBI database is missing either the nodes, '
                              'names or merge table',
                              RuntimeError)
            except sqlite3.DatabaseError:
                log_and_raise('Local SQLite database not valid',
                              OSError)
            else:
                self.use_local_ncbi = True
            finally:
                conn.close()

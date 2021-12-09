import os
import sqlite3
import simplejson
import appdirs
from .logger import LOGGER, FORMATTER, log_and_raise


class Resources:
    """Load and check validation resources

    The following resources need to be available to use the safedata_validator
    package:

        - locations: A path to a locations data file.
        - gbif_database: A path to an local SQLite copy of the GBIF backbone
            database. If this is None, the GBIF online API will be used instead.

    This class is used to locate and validate these resources. For local GBIF
    databases, the {meth}`~safedata_validator.config.Resources.get_get_local_gbif_conn`
    can be used to obtain a sqlite3 connection to the database. The returned
    connection should be closed when it is no longer needed.

    The paths can be set as arguments to initialise a Resource instance or stored
    in a JSON configuration file providing `location` and `gbif_database` keys.
    This file can be passed as `cfg_path` when creating an instance, but if no
    arguments are provided then an attempt is made to find and load configuration
    files in the user and then site config locations defined by the appdirs package.

    Args:
        locations: A path to a JSON formatted file containing location names
            and aliases
        gbif_database: A path to a local sqlite3 file containing the GBIF
            backbone database or None to use the GBIF online API.
        cfg_file: A path to a configuration file.

    Attributes:
        config_type: The method used to specify the resources. One of 'init_params',
            'init_file', 'user_config' or 'site_config'.
        self.config_file = config['config_file']
        locations: The path to the locations file
        gbif_database: The path to the GBIF database file or None
        use_local_gbif: Is a local file used or should the GBIF API be used
        valid_locations: The locations defined in the locations file
        location_aliases: Location aliases defined in the locations file
    """

    def __init__(self, locations=None, gbif_database=None, cfg_file=None):

        # First step is to try and populate locations and gbif_database
        # attributes either from the init params or from a config file.
        LOGGER.info('Configuring Resources', extra=dict(depth=0))
        FORMATTER.depth = 1

        if locations is not None:
            config = dict(locations=locations,
                          gbif_database=gbif_database,
                          config_type='init params',
                          config_file=None)
        elif cfg_file is not None:
            config = self._load_config(cfg_file, 'init file')
        else:
            # Look for user config
            cfg_file = os.path.join(appdirs.user_config_dir(),
                                    'safedata_validator', 'safedata_validator.json')
            config = self._load_config(cfg_file, 'user config')

            # If the user file does not exist, try a site config
            if config is None:
                cfg_file = os.path.join(appdirs.site_config_dir(),
                                        'safedata_validator', 'safedata_validator.json')
                config = self._load_config(cfg_file, 'site config')

            # If there is still nothing, then there is nothing left to try.
            if config is None:
                log_and_raise(LOGGER,
                              f'No init params or config files provided.',
                              RuntimeError)

        # Set attributes
        self.locations = config['locations']
        self.gbif_database = config['gbif_database']
        self.config_type = config['config_type']
        self.config_file = config['config_file']

        self.use_local_gbif = None
        self.valid_locations = None
        self.location_aliases = None

        # Report resource config location and type
        msg = f'Configuring resources from {self.config_type}'
        if self.config_file is not None:
            msg += self.config_file

        LOGGER.info(msg)

        # Validate the resources
        self._validate_locations()
        self._validate_gbif()

    @staticmethod
    def _load_config(cfg_path, cfg_type):
        """Load a configuration file

        This private static method attempts to load a JSON configuration file
        from a path.

        Returns:
             If the file does not exist, the function returns None. Otherwise,
             it returns a dictionary of config parameters.
        """

        # If the file does not exist, return None gracefully
        if not (os.path.exists(cfg_path) and os.path.isfile(cfg_path)):
            return None

        # Otherwise, there is a file, so try and use it and now raise if there
        # is a problem: don't skip over broken resource configurations
        try:
            with open(cfg_path, 'r') as json:
                config = simplejson.load(json)
        except simplejson.errors.JSONDecodeError:
            log_and_raise(LOGGER,
                          f'File {cfg_path} not JSON encoded.',
                          IOError)
        except IOError:
            log_and_raise(LOGGER,
                          f'Could not read {cfg_path}.',
                          IOError)

        # Check JSON keys and that locations is not null
        if {'locations', 'gbif_database'} != config.keys():
            log_and_raise(LOGGER,
                          f'Config keys not found in {cfg_path}.',
                          RuntimeError)

        if config['locations'] is None:
            log_and_raise(LOGGER,
                          f'Locations cannot be null in {cfg_path}.',
                          RuntimeError)

        config['config_file'] = cfg_path
        config['config_type'] = cfg_type

        return config

    def _validate_locations(self):
        """Validate and load a locations file

        This private function checks whether a locations path: exists, is a
        JSON file, and contains location and alias data. It populates the
        instance attributes

        Returns:
            None - updates instance.
        """

        LOGGER.info(f'Validating locations from {self.locations}')

        # Now check to see whether the locations file behaves as expected
        if not os.path.exists(self.locations) and not os.path.isfile(self.locations):
            log_and_raise(LOGGER,
                          'Local locations file not found',
                          IOError)

        try:
            loc_payload = simplejson.load(open(self.locations, mode='r'))
        except simplejson.errors.JSONDecodeError:
            log_and_raise(LOGGER,
                          'Local locations file not JSON encoded.',
                          IOError)

        # process the locations payload
        if {'locations', 'aliases'} != loc_payload.keys():
            log_and_raise(LOGGER,
                          'Locations data malformed',
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
        if self.gbif_database is None:
            LOGGER.info('Using GBIF online API to validate taxonomy')
            self.use_local_gbif = False
        else:
            LOGGER.info(f'Validating local GBIF database: {self.gbif_database}')

            # Does the provided path exist and is it a functional SQLite database
            # with a backbone table? Because sqlite3 can 'connect' to any path,
            # use a query attempt to reveal exceptions

            if not os.path.exists(self.gbif_database):
                log_and_raise(LOGGER,
                              'Local GBIF database not found',
                              IOError)

            try:
                conn = sqlite3.connect(self.gbif_database)
                _ = conn.execute('select count(*) from backbone;')
            except sqlite3.OperationalError:
                log_and_raise(LOGGER,
                              'Local GBIF database does not contain the backbone table',
                              RuntimeError)
            except sqlite3.DatabaseError:
                log_and_raise(LOGGER,
                              'Local SQLite database not valid',
                              IOError)
            else:
                self.use_local_gbif = True
            finally:
                conn.close()

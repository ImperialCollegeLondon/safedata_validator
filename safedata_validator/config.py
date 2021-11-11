from .logger import LOGGER
import simplejson
import appdirs
import os
import dataclasses
from enforce_typing import enforce_types


class config():
    """Load safedata_validator config

    The following safedata_validator variables are configurable:
    - locations: A path to a locally stored locations data file or the URL of a web service
         that provides the same data.
    - gbif_database: A path to an local SQLite copy of the GBIF backbone database.
    """

    def __init__(self, config_file=None):

        self.locations = None
        self.gbif_database = None

        if config_file is not None:
            LOGGER.info(f'Using provided config: {config_file}')
            self._load_cfg(config_file, 'provided')
            return

        user_config_dir = appdirs.user_config_dir('safedata_validator')
        user_config_file = os.path.join(user_config_dir, 'safedata_validator.json')
        if os.path.exists(user_config_file) and os.path.isfile(user_config_file):
            LOGGER.info(f'Using user config: {user_config_file}')
            self._load_cfg(user_config_file)
            return

        site_config_dir = appdirs.site_config_dir('safedata_validator')
        site_config_file = os.path.join(site_config_dir, 'safedata_validator.json')
        if os.path.exists(site_config_file) and os.path.isfile(site_config_file):
            LOGGER.info(f'Using site config: {site_config_file}')
            self._load_cfg(user_config_file, 'site')
            return
        
        LOGGER.error('No config found.')
        return

    def _load_cfg(self, config_file, file_desc='user'):
        try:
            with open(config_file, 'r') as json:
                config = simplejson.load(json)
            LOGGER.info(f'Loaded {file_desc} config')
        except IOError:
            LOGGER.error(f'Found {file_desc} config but could not load.')
            return

        if 'locations' in config:
            if os.path.exists(config['locations']) and os.path.isfile(config['locations']):
                self.locations = config['locations']
            else:
                LOGGER.error(f'Config locations path not valid.')

        if 'gbif_database' in config:
            if os.path.exists(config['gbif_database']) and os.path.isfile(config['gbif_database']):
                self.gbif_database = config['gbif_database']
            else:
                LOGGER.error(f'Config gbif_database path not valid.')






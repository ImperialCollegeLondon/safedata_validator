from openpyxl import worksheet
from shapely import wkt
from shapely.errors import WKTReadingError

from safedata_validator.logger import LOGGER, FORMATTER, CH, loggerinfo_push_pop
from safedata_validator.validators import (GetDataFrame, IsLower, IsNotBlank, IsNumber,
                                           IsString, HasDuplicates, IsLocName)
from safedata_validator.extent import Extent

from safedata_validator.resources import Resources


class Locations:

    def __init__(self, resources: Resources) -> None:
        
        """
        Attempts to load and check the contents of the Locations worksheet and
        compile the geographic extent of the locations used. The values in the
        data file are validated against the locations data loaded when the Dataset
        was initialised.

        Args:
            worksheet:
            resources:

        Returns:
            A Locations instance
        """

        self.locations = set()
        self.location_index = []
        self.extents = {'latitude': Extent('latitude', (float, int)),
                        'longitude': Extent('longitude', (float, int))}
        self.locations_used = set()

        self.valid_locations = resources.valid_locations
        self.location_aliases = resources.location_aliases
        self.known_loc_names = set(list(resources.valid_locations.keys()) +
                                   list(resources.location_aliases.keys()))

    @loggerinfo_push_pop('Loading Locations worksheet')
    def load(self, worksheet: worksheet):

        """
        Attempts to load and check the contents of a Locations worksheet and
        compile the geographic extent of the locations used. 

        Args:
            worksheet:
        """

        start_errors = CH.counters['ERROR']

        # Load the locations data frame - which runs header checks
        LOGGER.info("Reading location data")
        FORMATTER.push()
        dframe = GetDataFrame(worksheet)

        # Dupe headers likely cause serious issues, so stop
        if 'duplicated' in dframe.bad_headers:
            LOGGER.error('Cannot parse locations with duplicated headers')
            return

        # Reduce to lower case 
        # TODO - not trapping dupes that are only case
        #        differences. Do that in GetDataFrame?
        headers = IsLower(dframe.headers).values

        # Check location names are available
        if 'location name' not in headers:
            LOGGER.error('Location name column not found')
            return

        # Get dictionaries of values for each row
        locs = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]

        # Split up old and new if there are there any new ones?
        if 'new' in headers:
            
            # Check the New column is just yes, no
            new_vals = IsLower([rw['new'] for rw in locs])
            new_vals = IsNotBlank(new_vals)
            if not new_vals:
                LOGGER.error('New locations field contains blank rows.')

            # check only yes or no entries
            valid_new = {'yes', 'no'}
            bad_new = set(new_vals) - valid_new
            if bad_new:
                LOGGER.error('New locations field contains values other than yes and no: ',
                             extra={'join': bad_new})
        
            known_locs = [rw for rw in locs 
                          if isinstance(rw['new'], str) and rw['new'].lower() == 'no']
            new_locs = [rw for rw in locs 
                          if isinstance(rw['new'], str) and rw['new'].lower() == 'yes']
        else:
            new_locs = []
            known_locs = locs

        if new_locs:
            self.validate_new_locations(new_locs)
        
        if known_locs:
            known_loc_names = [lc['location name'] for lc in known_locs]
            self.validate_known_locations(known_loc_names)
        
        # summary of processing
        n_errors = CH.counters['ERROR'] - start_errors

        if n_errors > 0:
            LOGGER.info('Locations contains {} errors'.format(n_errors))
        else:
            LOGGER.info('{} locations loaded correctly'.format(len(self.locations)))

    @loggerinfo_push_pop('Checking new locations')
    def add_new_locations(self, locs):
        """This method takes list of dictionaries giving the details of
        new locations to be added to the instance. These dictionaries
        should contain keys `location name` and `type` and then at least 
        one of  both `latitude` _and_ `longitude` or `wkt`.
        """

        # Validation - TODO check locs is a list of dicts 
        # - Do all the dicts have the same keys
        loc_keys = set([tuple(k.keys()) for k in locs])
        if len(loc_keys) > 1:
            # TODO - should this be CRITICAL? Check error types.
            LOGGER.error('Inconsistent keys in add_new_locations')
            return
        
        # - Do they provide location names...
        loc_keys = loc_keys.pop()
        if 'location name' not in loc_keys:
            LOGGER.error("No location name entries in add_new_locations")
            return

        # - ... and are any of those names blank ...
        loc_names = [itm['location name'] for itm in locs]
        loc_names = IsNotBlank(loc_names, keep_failed=False)
        if not loc_names:
            LOGGER.error('Location names contains empty cells or whitespace text')

        # - ... or not strings. Only allow strings for new locations - known 
        #   locations get a pass for integer site codes but not here.
        loc_names = IsString(loc_names)
        if not loc_names:
            LOGGER.error('New location names include non-string values: ',
                         extra={'join': loc_names.failed})

        # Look for duplicated names in inputs - this includes all types and so
        # could give messy information but can't detect duplicates late on
        # cleaned names because we're dealing with sets by that point and there
        # can be no duplication in sets
        dupes = HasDuplicates(loc_names)
        if dupes:
            LOGGER.error('New location names contain duplicated values: ', 
                         extra={'join': dupes})

        # Look for new names that duplicate known names
        duplicates_existing = [rw['location name'] for rw in locs
                               if rw['location name'] in self.known_loc_names]

        if duplicates_existing:
            LOGGER.error('New location names duplicate known names and aliases: ',
                            extra={'join': duplicates_existing})

        # Type is required - used to indicate the kind of location in the absence
        # of any actual geodata
        if 'type' not in loc_keys:
            LOGGER.error('New locations do not provide the location type')
            for this_loc in locs:
                this_loc['type'] = "MISSING"
        else:
            # get lowercase types
            geo_types = set(IsLower([vl['type'] for vl in locs]))

            # Handle blanks
            geo_types = IsNotBlank(geo_types, keep_failed=False)
            if not geo_types:
                LOGGER.error('Types for new locations contains blank or whitespace entries.')

            # Handle unknown geo types
            valid_geotypes = {'point', 'linestring', 'polygon', 'transect', 'area'}
            bad_geo_types = set(geo_types) - valid_geotypes
            if bad_geo_types:
                LOGGER.error('New locations include unknown location types: ',
                                extra={'join': bad_geo_types})
            
        # Geometry information
        # Record which geom columns are present
        if ('latitude' in loc_keys) ^ ('longitude' in loc_keys):
            LOGGER.error('New locations should either latitude _and_ longitude or neither')

        lonlat_provided = 'latitude' in loc_keys and 'longitude' in loc_keys
        wkt_provided = 'wkt' in loc_keys

        # TODO - supplying both is not an error, and probably shouldn't be. WKT
        # takes priority in the index when both are present, but not testing for
        # congruence at present

        # Look for geographic data (even if it is just the explicit statement
        # that none is available using NA)
        if not (lonlat_provided or wkt_provided):
            LOGGER.error('New locations reported: you must provide Lat/Long or WKT,'
                         'using NA explicitly when this data is missing.')
        else:
            # Check Lat Long and WKT using check_field_geo to validate the values,
            # since this method automatically updates the dataset extents. 

            if lonlat_provided:

                LOGGER.info('Validating lat / long data')

                for axs in ['latitude', 'longitude']:
                    
                    # Allow NAs for unknown location points
                    axs_vals = [vl[axs] for vl in locs if vl[axs] != 'NA']
                    axs_vals = IsNotBlank(axs_vals, keep_failed = False)
                    if not axs_vals:
                        LOGGER.error(f'Blank {axs} values for new locations: use NA.')
                    
                    # Check for data types _here_ to keep interpretable errors
                    # TODO - maybe simplify Extent objects to _assume_ types.
                    #        Probably no but there is some duplication of effort here
                    axs_vals = IsNumber(axs_vals, keep_failed = False)
                    if not axs_vals:
                        LOGGER.error(f'Non-numeric {axs} values for new locations: ',
                                     extra={'join': axs_vals.failed})

                    # Update extents
                    if axs_vals.values:
                        self.extents[axs].update(axs_vals)

            if wkt_provided:

                LOGGER.info('Validating WKT data')
                blank_wkt = []
                non_string_wkt = []
                bad_wkt = []

                for this_new_loc in locs:
                    
                    if this_new_loc['wkt'] is None:
                        blank_wkt.append(this_new_loc['location name'])
                    elif not isinstance(this_new_loc['wkt'], str):
                        non_string_wkt.append(this_new_loc['location name'])
                    elif this_new_loc['wkt'].isspace():
                        blank_wkt.append(this_new_loc['location name'])
                    elif this_new_loc['wkt'] == 'NA':
                        pass
                    
                    else:

                        # Run the potential WKT through the parser
                        try:
                            this_new_geom = wkt.loads(this_new_loc['wkt'])
                        except WKTReadingError:
                            bad_wkt.append(this_new_loc['location name'])
                        else:
                            # Is it a valid 2D geom
                            if not this_new_geom.is_valid or this_new_geom.has_z:
                                bad_wkt.append(this_new_loc['location name'])
                            # Store the extents to check for sensible coordinates
                            bnds = this_new_geom.bounds
                            self.extents['latitude'].update([bnds[1], bnds[3]])
                            self.extents['longitude'].update([bnds[0], bnds[2]])

                if blank_wkt:
                    LOGGER.error('Blank WKT values for new locations: use NA.')

                if non_string_wkt:
                    LOGGER.error('WKT values for new location not a string: ',
                                 extra={'join': non_string_wkt})

                if bad_wkt:
                    LOGGER.error('WKT information badly formatted, not geometrically valid or 3D: ',
                                    extra={'join': bad_wkt})

        # new location names
        # - test for duplicated names to already added values
        dupes = [lc for lc in loc_names if lc in self.locations]
        if dupes:
            LOGGER.error('Location names already added to Location instance: ',
                         extra={'join': dupes})
        
        self.locations.update(loc_names)

        for this_new_loc in locs:
            
            if wkt_provided and this_new_loc['wkt'] != 'NA':
                geom = this_new_loc['wkt']
            elif lonlat_provided and ((this_new_loc['latitude'] != 'NA') and (this_new_loc['longitude'] != 'NA')):
                geom = 'Point({longitude} {latitude})'.format(**this_new_loc)
            else:
                geom = None

            self.location_index.append((this_new_loc['location name'], True,
                                        this_new_loc['type'], geom))

    @loggerinfo_push_pop('Checking known locations')
    def add_known_locations(self, loc_names: list):
        """This method takes a list of values and tries to validate those
        values against known locations from the loaded resources. The values
        are expected to be strings"""

        # Check for blanks
        loc_names = IsNotBlank(loc_names, keep_failed=False)
        if not loc_names:
            LOGGER.error('Location names contains empty cells or whitespace text')

        # Look for duplicated values in names - this includes all types and so
        # could give messy information but can't detect duplicates late on
        # cleaned names because we're dealing with sets by that point and there
        # can be no duplication
        dupes = HasDuplicates(loc_names)
        if dupes:
            LOGGER.error('Added names contain duplicated values: ', 
                         extra={'join': dupes})

        # Validate and standardise types - strings or integer codes.
        loc_names = IsLocName(loc_names, keep_failed=False)

        if not loc_names:
            LOGGER.error('Location names contains values that are not strings or integers: ',
                          extra={'join': loc_names.failed})

        # Enforce strings and check loc names exist
        loc_names = set([str(v) for v in loc_names])
        unknown = loc_names - self.known_loc_names
        if unknown:
            LOGGER.error('Unknown locations found: ', extra={'join': unknown, })

        # are aliases being used?
        aliased = loc_names & set(self.location_aliases.keys())
        if aliased:
            LOGGER.warning('Locations aliases used. Maybe change to primary location names: ',
                           extra={'join': aliased})

        # Get the bounding box of known locations and aliased locations
        bbox_keys = (loc_names - (unknown | aliased)) | {self.location_aliases[ky] for ky in aliased}

        # get the extents of known unaliased locations
        if bbox_keys:
            bbox = [vl for ky, vl in list(self.valid_locations.items()) if ky in bbox_keys]
            bbox = list(zip(*bbox))
            self.longitudinal_extent.update((min(bbox[0]), max(bbox[1])))
            self.latitudinal_extent.update((min(bbox[2]), max(bbox[3])))

        # Update location names and index
        # - test for duplicated names to already added values
        dupes = [lc for lc in loc_names if lc in self.locations]
        if dupes:
            LOGGER.error('Location names already added to Location instance: ',
                         extra={'join': dupes})

        self.locations.update(loc_names)
        index_entries = [(lc, False, None, None) for lc in loc_names]
        self.location_index.append(index_entries)
    
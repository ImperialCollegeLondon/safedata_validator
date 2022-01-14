
from dataclasses import dataclass
from safedata_validator.logger import LOGGER, FORMATTER, CH, loggerinfo_push_pop
from safedata_validator.validators import (GetDataFrame, IsLower, IsNotBlank,
                                           IsString, HasDuplicates)
from safedata_validator.extent import Extent

from shapely import wkt
from shapely.errors import WKTReadingError



class Locations:

    def __init__(self, resources) -> None:
        
        self.locations = set()
        self.location_index = []
        self.latitudinal_extent = Extent('latitude', (float, int))
        self.longitudinal_extent = Extent('latitude', (float, int))
        # self.new = False
        # self.lonlat = False
        # self.wkt = False
        self.locations_used = set()

        self.valid_locations = resources.valid_locations
        self.location_aliases = resources.location_aliases
        self.known_loc_names = set(list(resources.valid_locations.keys()) +
                                   list(resources.location_aliases.keys()))

    @loggerinfo_push_pop('Loading Locations worksheet')
    def load(self, worksheet, resources):

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
        # TODO - not trapping dupes that are only case differences. Do that in GetDataFrame?
        headers = IsLower(dframe.headers).values

        # Check location names are available
        if 'location name' not in headers:
            LOGGER.error('Location name column not found')
            return

        # Record which geom columns are present
        if ('latitude' in headers) ^ ('longitude' in headers):
            LOGGER.error('Provide both or neither of latitude and longitude')

        self.lonlat = 'latitude' in headers and 'longitude' in headers
        self.wkt = 'wkt' in headers
        self.new = 'new' in headers

        # Location name cleaning
        loc_names = dframe.data_columns[headers.index('location name')]

        loc_names = IsNotBlank(loc_names)

        # Check for blanks
        if not loc_names:
            LOGGER.error('Blank location names: ',
                         extra={'join': loc_names.failed})

        # Convert non-strings to strings (e.g. numeric site codes)
        loc_names = IsString(loc_names)

        # look for duplicates
        unique = HasDuplicates(loc_names)
        if unique:
            LOGGER.error('Duplicated location names: ',
                         extra={'join': unique.duplicated})

        dframe.data_columns[headers.index('location name')] = loc_names

        # VALIDATE LOCATIONS

        # Get dictionaries of values for each row
        locs = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]

        # Get existing location names and aliases - new location names must
        # not appear in it and existing locations must.
        existing_loc_names = set(list(resources.valid_locations.keys()) +
                                 list(resources.location_aliases.keys()))

        # Split up old and new if there are there any new ones?
        if self.new:
            # Check the New column is just yes, no
            is_new = set(IsLower([rw['new'] for rw in locs]))
            is_new = IsNotBlank(is_new)
            if not is_new:
                LOGGER.error('New locations field contains blank rows.')

            # check only yes or no entries
            valid_new = {'yes', 'no'}
            bad_new = is_new - valid_new
            if bad_new:
                LOGGER.error('New field contains values other than yes and no: ',
                             extra={'join': bad_new})

            # extract the new and old locations - dropping None
            new_locs = [rw for rw in locs if rw['new'].lower() == 'yes']
            known_locs = [rw for rw in locs if rw['new'].lower() == 'no']
        else:
            new_locs = []
            known_locs = locs

        if new_locs:
            self.validate_new_locations(new_locs)
        
        if known_locs:
            known_loc_names = {lc['location name'] for lc in known_locs}
            self.validate_known_locations(known_loc_names)
        
        # summary of processing
        n_errors = CH.counters['ERROR'] - start_errors

        if n_errors > 0:
            LOGGER.info('Locations contains {} errors'.format(n_errors))
        else:
            LOGGER.info('{} locations loaded correctly'.format(len(self.locations)))

    @loggerinfo_push_pop('Checking new locations')
    def validate_new_locations(self, locations):


        # Type is required - used to indicate the kind of location in the absence
        # of any actual geodata
        if 'type' not in headers:
            LOGGER.error('New locations reported but Type field missing')

        else:
            # get lowercase types
            geo_types = set(IsLower([vl['type'] for vl in new_locs]))

            # Handle blanks
            if None in geo_types:
                LOGGER.error('Types for new locations contains blank entries.')
                geo_types -= {None}

            # Handle unknown geo types
            valid_geotypes = {'point', 'linestring', 'polygon', 'transect', 'area'}
            bad_geo_types = geo_types - valid_geotypes
            if bad_geo_types:
                LOGGER.error('Unknown location types: ',
                                extra={'join': bad_geo_types})

        # Look for duplicate names
        duplicates_existing = [rw['location name'] for rw in new_locs
                                if rw['location name'] in existing_loc_names]

        if duplicates_existing:
            LOGGER.error('New location names duplicate existing names and aliases: ',
                            extra={'join': duplicates_existing})

        # Look for geographic data (even if it is just the explicit statement
        # that none is available using NA)
        if not self.lonlat or self.wkt:
            LOGGER.error('New locations reported: you must provide Lat/Long or WKT')

        else:
            # Check Lat Long and WKT using check_field_geo to validate the values,
            # since this method automatically updates the dataset extents.

            if self.lonlat:

                LOGGER.info('Validating lat / long data')

                for axs in ['latitude', 'longitude']:

                    axs_vals = [vl[axs] for vl in new_locs if vl[axs] != 'NA']
                    axs_vals = IsNotBlank(axs_vals)

                    if not axs_vals:
                        LOGGER.error(f'Blank {axs} values for new locations: use NA.')

                    # Create an extent object to check values and report extents
                    self.extents[axs].update(axs_vals)

            if self.wkt:

                LOGGER.info('Validating WKT data')

                bad_wkt = []

                # Get locs where WKT not explicitly missing
                wkt_data = [l for l in new_locs if l['wkt'] != 'NA']
                wkt_found = [l for l in wkt_data if l['wkt'] is not None]

                if len(wkt_found) < len(wkt_data):
                    LOGGER.error('WKT field contains blanks for new locations, use NA')

                for this_new_loc in wkt_found:

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

                if bad_wkt:
                    LOGGER.error('WKT information badly formatted, not geometrically valid or 3D: ',
                                    extra={'join': bad_wkt})

        # new location names
        new_loc_names = {rw['location name'] for rw in new_locs}

        self.locations.update(new_loc_names)

        for this_new_loc in new_locs:

            if self.wkt and this_new_loc['wkt'] != 'NA':
                geom = this_new_loc['wkt']
            elif self.lonlat and ((this_new_loc['latitude'] != 'NA') and (this_new_loc['longitude'] != 'NA')):
                geom = 'Point({longitude} {latitude})'.format(**this_new_loc)
            else:
                geom = None

            self.location_index.append((this_new_loc['location name'], True,
                                        this_new_loc['type'], geom))


    @loggerinfo_push_pop('Checking known locations')
    def validate_known_locations(self, loc_names):

        # check names exist
        unknown = loc_names - self.known_loc_names
        if unknown:
            LOGGER.error('Unknown locations found: ', extra={'join': unknown, })

        # are aliases being used?
        aliased = loc_names & set(self.location_aliases.keys())
        if aliased:
            LOGGER.warn('Locations aliases used. Maybe change to primary location names: ',
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
        self.locations.update(loc_names)
        index_entries = [(lc, False, None, None) for lc in loc_names]
        self.location_index.append(index_entries)
    
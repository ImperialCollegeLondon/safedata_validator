#!/usr/bin/python
# -*- coding: UTF8 -*-

"""
Module containing code to verify the format of a SAFE project Excel dataset.

The functions are written to extract as much information as possible from
checking a file so, rather than raising an error and halting at the first
fault, functions are written to check as much as they can on the inputs.

Metadata is stored in a Dataset object, which also keeps a record of any
information and warnings issued during processing. These are optionally
printed to screen during processing but can also be retrieved from the
Dataset object.

Taxon names and locations are both validated:
i)  By default, taxonomic names are validated against the GBIF backbone taxonomy
    using API queries over the internet. A user can download and build a local
    sqlite copy of the database to avoid the need for a web connection and to
    improve the speed of validation.

ii) Locations are validated against a list of valid locations. By default, this
    is downloaded from a web service provided from the SAFE website, but a local
    file can be provided for off line use.
"""

from __future__ import print_function
import os
import datetime
import argparse
import re
from collections import Counter
from StringIO import StringIO
import numbers
import openpyxl
from openpyxl import utils
import requests
import simplejson
import sqlite3

# define some regular expressions used to check validity
RE_ORCID = re.compile(r'[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{4}')
RE_EMAIL = re.compile(r'\S+@\S+\.\S+')
RE_NAME = re.compile(r'[^,]+,[ ]?[^,]+')
RE_WSPACE_ONLY = re.compile(r'^\s*$')
# note that logically the next one also catches whitespace at both ends
RE_WSPACE_AT_ENDS = re.compile(r'^\s+\w+|\w+\s+$')
RE_DMS = re.compile(r'[°\'"dms’”]+')

"""
Some simple helper functions
"""


def is_blank(value):
    """
    Helper function that checks if a value is None or contains only whitespace
    Args:
        value: A single value

    Returns:
        Boolean
    """

    return (value is None) or (RE_WSPACE_ONLY.match(unicode(value)) is not None)


def is_padded(value):
    """
    Helper function that checks if the value is padded with whitespace
    Args:
        value: Contents of a cell

    Returns:
        Boolean
    """

    return (value is not None) and (RE_WSPACE_AT_ENDS.match(value))


def duplication(data):
    """
    Looks for duplication in a list of strings, returning the duplicated elements.

    Args:
        data: A list of strings

    Returns:
        A list of duplicated elements
    """

    return [ky for ky, vl in Counter(data).iteritems() if vl > 1]


def is_integer_string(txt):
    """
    Checks if a string value can represent an integer.
    Args:
        txt: A string

    Returns:
        A boolean.
    """
    try:
        int(txt)
        return True
    except ValueError:
        return False


def all_numeric(data):
    """
    Tests if all values in an iterable are numeric.
    Args:
        data: An iterable claiming to contain numbers

    Returns:
        A list containing a boolean indicating whether all of the
        values were numbers and then a list of the genuine numeric values
    """

    nums = [vl for vl in data if isinstance(vl, numbers.Number)]

    return [len(nums) == len(data), nums]


class DataWorksheet(object):

    """
    This is just a container for the metadata on a data worksheet.
    It basically only exists to set defaults on a data worksheet and
    to provide dot notation rather than using dictionary labels.

    Args:
        meta: The dictionary of worksheet name, title and description
        loaded from the summary worksheet.
    """

    def __init__(self, meta):

        # set defaults
        self.name = meta['name']
        self.description = meta['description']
        self.title = meta['title']
        self.max_row = 0
        self.max_col = 0
        self.descriptors = None
        self.taxa_fields = None
        self.field_name_row = None
        self.fields = []


class Dataset(object):

    """
    This class provides methods to load and store the metadata associated
    with a dataset stored in a SAFE project formatted Excel file.

    Args:
        filename: The Excel file from which to read the Dataset
        verbose: Should the reporting print to the screen during runtime
        gbif_database: A local path to an sqlite database containing the GBIF database
        to be used instead of the web API.
    """

    def __init__(self, filename, verbose=True, gbif_database=None):

        try:
            self.workbook = openpyxl.load_workbook(filename=filename, data_only=True,
                                                   read_only=True)
        except IOError:
            raise IOError('Could not open file {}'.format(filename))

        self.filename = os.path.basename(filename)

        # report tracking
        self.message_header = "Checking file '{}'".format(filename)
        self.messages = []
        self.verbose = verbose
        self.kinds = {'info': '-', 'hint': '?', 'warn': '!'}
        self.n_warnings = 0
        if verbose:
            print(self.message_header)

        # basic info
        self.sheet_names = set(self.workbook.get_sheet_names())

        # initialise data contents
        self.project_id = None
        self.authors = []
        self.title = None
        self.description = None
        self.dataworksheet_summaries = []
        self.access = 'Open'
        self.embargo_date = None
        self.publication_doi = []
        self.dataworksheets = []
        self.keywords = []
        self.temporal_extent = None
        self.latitudinal_extent = None
        self.longitudinal_extent = None
        self.locations = set()
        self.locations_used = set()
        self.taxon_names = set()
        self.taxon_names_used = set()
        self.taxon_index = set()

        # Setup the taxonomy validation mechanism.
        if gbif_database is None:
            self.info('Using GBIF online API to validate taxonomy', 0)
            self.use_local_gbif = False
            self.gbif_conn = None
        else:
            # If a GBIF file is provided, does the file exist?
            if not os.path.exists(gbif_database):
                self.info('Local GBIF database not found, defaulting to online API', 0)
                self.use_local_gbif = False
                self.gbif_conn = None

            # And does the file behave like i) a sqlite database ii) containing a table 'backbone'
            try:
                # can connect sqlite to any existing path, but only attempts to query it
                # throw up exceptions
                conn = sqlite3.connect(gbif_database)
                _ = conn.execute('select count(*) from backbone;')
            except sqlite3.OperationalError:
                self.info('Local GBIF database does not contain the backbone table, defaulting '
                          'to online API', 0)
                self.use_local_gbif = False
                self.gbif_conn = None
            except sqlite3.DatabaseError:
                self.info('Local SQLite database not valid, defaulting to online API', 0)
                self.use_local_gbif = False
                self.gbif_conn = None
            else:
                self.info('Using local GBIF database to validate taxonomy', 0)
                self.use_local_gbif = True
                self.gbif_conn = conn
                self.gbif_conn.row_factory = sqlite3.Row

    # Logging methods
    def info(self, message, level=0):
        """
        Adds an information message to the instance.
        Args:
            message: The information message text.
            level: The message indent level
        """
        msg = ('info', message, level)
        self.messages.append(msg)
        if self.verbose:
            self.print_msg(*msg)

    def warn(self, message, level=0, join=None, as_repr=False):
        """
        Adds a warning message to the instance, optionally adding
        a joined set of values onto the end, possibly converting them
        to a string representation to capture formatting more clearly

        Args:
            message: The warning message text.
            level: The message indent level
            join: A list of values to join and append to the warning
            as_repr: Should join values be converted to a representation
        """
        if join and as_repr:
            message += ', '.join([repr(unicode(vl)) for vl in join])
        elif join:
            message += ', '.join(join)

        msg = ('warn', message, level)
        self.messages.append(msg)
        self.n_warnings += 1
        if self.verbose:
            self.print_msg(*msg)

    def hint(self, message, level=0, join=None):
        """
        Adds an hint message to the instance. This is used to indicate where
        something might need to be checked but isn't necessarily an error
        Args:
            message: The hint text.
            level: The hint indent level
            join: A list of values to join and append to the warning
        """
        if join:
            message += ', '.join(join)

        msg = ('hint', message, level)
        self.messages.append(msg)
        if self.verbose:
            self.print_msg(*msg)

    def print_msg(self, kind, msg, level):
        """
        Prints a message to the screen
        Args:
            kind: The message kind
            msg: The message text
            level: The message indent level
        """
        print('  ' * level + self.kinds[kind] + ' ' + msg)

    def report(self):
        """
        Formats all the messages into a report.
        Returns:
            A StringIO() instance containing the report text
        """
        report = StringIO()
        report.writelines([self.message_header + '\n'])
        msgs = ['  ' * lv + self.kinds[kn] + msg + '\n' for kn, msg, lv in self.messages]
        report.writelines(msgs)
        return report

    # Utility method to update extents
    def update_extent(self, extent, val_type, which):
        """
        Updates one of the three dataset extents
        Args:
            extent: A two tuple (min, max) for one of the extents
            val_type: The type that the values in the extent should be
            which: Which extent to update
        """

        if (not isinstance(extent, tuple) or
                len(extent) != 2 or
                not all([isinstance(vl, val_type) for vl in extent]) or
                extent[0] > extent[1]):
            raise ValueError('extent must be an increasing two-tuple of values of '
                             'type {}'.format(val_type))

        current_vals = getattr(self, which)
        if current_vals is None:
            setattr(self, which, extent)
        else:
            extent = zip(current_vals, extent)
            setattr(self, which, (min(extent[0]), max(extent[1])))

    # Main methods to populate class attributes
    def load_summary(self, validate_doi=False):

        """
        Checks the information in the summary worksheet and looks for the metadata and
        dataset worksheets. The function is intended to try and get as much information
        as possible from the worksheet: the dictionary of metadata returned will have
        None for any missing data, which should be handled by downstream code.

        Args:
            validate_doi: Check any publication DOIs, requiring a web connection.
        """

        # try and get the summary worksheet
        self.info("Checking Summary worksheet")
        start_warn = self.n_warnings
        try:
            worksheet = self.workbook['Summary']
        except KeyError:
            self.warn("Summary worksheet not found, moving on.")
            return None

        # load dictionary of summary information block, allowing for multiple
        # columns for fields (data compilers, dataset sheets).
        # Note that worksheet.columns can load some odd invisible empty cells:
        # it doesn't guarantee an equal number of cells in each column. So 
        # use the squared range.
        mxcl = worksheet.max_column
        mxrw = worksheet.max_row

        # populate block of rows of cell objects
        data = list(worksheet.get_squared_range(1, 1, mxcl, mxrw))

        # transpose and get values
        data = zip(*data)
        cols = tuple(tuple(cell.value for cell in col) for col in data)

        # convert into dictionary
        hdrs = cols[0]
        values = zip(*cols[1:])
        summary = {ky: vl for ky, vl in zip(hdrs, values)}

        # Check the minimal keys are expected
        required = {"SAFE Project ID", "Access status", "Title", "Description",
                    "Author name", "Author email", "Author affiliation", "Author ORCID",
                    "Worksheet name", "Worksheet title", "Worksheet description", 'Keywords'}

        found = set(summary.keys())

        # don't bail here - try and get as far as possible
        if not found.issuperset(required):
            self.warn('Missing metadata fields: ', 1, join=required - found)

        # # check for contents
        # if any([set(x) == {None} for x in summary.values()]):
        #     self.warn('Metadata fields with no information.', 1)

        # CHECK PROJECT ID
        if 'SAFE Project ID' not in summary:
            self.warn('SAFE Project ID missing', 1)
        elif not isinstance(summary['SAFE Project ID'][0], long):
            self.warn('SAFE Project ID is not an integer.', 1)
        else:
            self.project_id = summary['SAFE Project ID'][0]

        # CHECK DATASET TITLE
        if 'Title' not in summary:
            self.warn('Dataset title row missing', 1)
        elif is_blank(summary['Title'][0]):
            self.warn('Dataset title is blank', 1)
        else:
            self.title = summary['Title'][0]

        # CHECK DATASET DESCRIPTION
        if 'Description' not in summary:
            self.warn('Dataset description row missing', 1)
        elif is_blank(summary['Description'][0]):
            self.warn('Dataset description is blank', 1)
        else:
            self.description = summary['Description'][0]

        # CHECK ACCESS STATUS AND EMBARGO DETAILS
        if 'Access status' not in summary:
            self.warn('Access status missing', 1)
        elif summary['Access status'][0] not in ['Open', 'Embargo']:
            self.warn('Access status must be Open or Embargo '
                      'not {}'.format(summary['Access status'][0]), 1)
        elif summary['Access status'] == 'Embargo':
            self.access = 'Embargo'
            if 'Embargo date' not in summary:
                self.warn('Dataset embargoed but embargo date row missing.', 1)
            elif is_blank(summary['Embargo'][0]):
                self.warn('Dataset embargo date  is blank', 1)
            else:
                embargo_date = summary['Embargo date'][0]
                now = datetime.datetime.now()
                if not isinstance(embargo_date, datetime.datetime):
                    self.warn('Embargo date not formatted as date.', 1)
                elif embargo_date < now:
                    self.warn('Embargo date is in the past.', 1)
                elif embargo_date > now + datetime.timedelta(days=2 * 365):
                    self.warn('Embargo date more than two years in the future.', 1)
                else:
                    self.embargo_date = embargo_date.date().isoformat()

        # CHECK KEYWORDS
        if 'Keywords' not in summary:
            self.warn('Dataset keywords row missing', 1)
        elif all([is_blank(kywd) for kywd in summary['Keywords']]):
            self.warn('No keywords provided', 1)
        else:
            # drop any blanks
            self.keywords = [vl for vl in summary['Keywords'] if not is_blank(vl)]
            if len(self.keywords) == 1 and ',' in self.keywords[0]:
                self.warn('Put keywords in separate cells, not comma delimited in one cell', 1)

        # CHECK FOR PUBLICATION DOIs if any are provided
        if 'Publication DOI' in summary:
            pub_doi = [vl for vl in summary['Publication DOI'] if not is_blank(vl)]
            # check formatting - basically make sure they have a proxy URL
            doi_is_url = [vl.startswith('https://doi.org/') for vl in pub_doi]
            if not all(doi_is_url):
                self.warn('Please provide publication DOIs as a URL: https://doi.org/...', 1)

            if validate_doi:
                for doi, is_doi in zip(pub_doi, doi_is_url):
                    if is_doi:
                        api_call = 'https://doi.org/api/handles/{}'.format(doi[16:])
                        r = requests.get(api_call)
                        if r.json()['responseCode'] != 1:
                            self.warn('DOI not found: {}'.format(doi), 1)

            self.publication_doi = pub_doi

        # CHECK AUTHORS
        # Get the set of author fields and create blank entries for any missing fields
        # to simplify handling the error checking.
        author_keys = ['Author name', 'Author affiliation', 'Author email', 'Author ORCID']
        authors = {k: summary[k] if k in summary else [None] * (len(cols) - 1)
                   for k in author_keys}

        # - switch in zenodo keys
        zenodo_keys = ['name', 'affiliation', 'email', 'orcid']
        for old, new in zip(author_keys, zenodo_keys):
            authors[new] = authors[old]
            authors.pop(old)

        # rotate and remove completely blank columns
        authors_list = zip(*authors.values())
        authors_list = [x for x in authors_list if x != tuple([None] * 4)]
        # return to the original orientation
        authors = {k: v for k, v in zip(authors, zip(*authors_list))}

        # now check the contents
        # i) Names
        author_names = [unicode(vl) for vl in authors['name'] if vl is not None]
        if len(author_names) < len(authors['name']):
            self.warn('Missing author names', 1)
        if author_names:
            bad_names = [vl for vl in author_names if not RE_NAME.match(vl)]
            if bad_names:
                self.warn('Author name not formatted as last_name, first_names: ', 1,
                          join=bad_names, as_repr=True)

        # ii) Affiliations (no regex checking)
        blank_affil = [is_blank(vl) for vl in authors['affiliation']]
        if any(blank_affil):
            self.hint('Missing affiliations - please provide if available', 1)

        # iii) Email
        author_emails = [unicode(vl) for vl in authors['email'] if vl is not None]
        if len(author_emails) < len(authors['email']):
            self.hint('Missing author emails - please provide if available', 1)
        if author_emails:
            bad_emails = [vl for vl in author_emails if not RE_EMAIL.match(vl)]
            if bad_emails:
                self.warn('Email not properly formatted: ', 1, join=bad_emails, as_repr=True)

        # iii) ORCiD (not mandatory)
        author_orcid = [unicode(vl) for vl in authors['orcid'] if vl is not None]
        if len(author_orcid) < len(authors['orcid']):
            self.hint('Missing ORCiDs, consider adding them!', 1)
        if author_orcid:
            bad_orcid = [vl for vl in author_orcid if not RE_ORCID.match(vl)]
            if bad_orcid:
                self.warn('ORCID not properly formatted: ', 1, join=bad_orcid)

        # and finally store as a dictionary per author
        self.authors = [dict(zip(authors.keys(), vals)) for vals in zip(*authors.values())]

        # CHECK DATA WORKSHEETS
        ws_keys = ['Worksheet name', 'Worksheet title', 'Worksheet description']
        data_worksheets = {k: summary[k] if k in summary else [None] * (len(cols) - 1)
                           for k in ws_keys}

        # - switch in short keys
        short_ws_keys = ['name', 'title', 'description']
        for old, new in zip(ws_keys, short_ws_keys):
            data_worksheets[new] = data_worksheets[old]
            data_worksheets.pop(old)

        # rotate and remove completely blank columns
        data_worksheets_list = zip(*data_worksheets.values())
        data_worksheets_list = [x for x in data_worksheets_list if x != tuple([None] * 3)]
        # return to the original orientation
        data_worksheets = {k: v for k, v in zip(data_worksheets, zip(*data_worksheets_list))}

        # now check the contents
        # i) Names
        ws_names = [vl for vl in data_worksheets['name'] if vl is not None]
        if len(ws_names) < len(data_worksheets['name']):
            self.warn('Missing worksheet names', 1)
        if ws_names:
            bad_names = [vl for vl in ws_names if vl not in self.sheet_names]
            if bad_names:
                self.warn('Worksheet names not found in workbook: ', 1,
                          join=bad_names, as_repr=True)

        # ii) Titles
        blank_names = [is_blank(vl) for vl in data_worksheets['title']]
        if any(blank_names):
            self.warn('Missing worksheet title', 1)

        # ii) Descriptions
        blank_desc = [is_blank(vl) for vl in data_worksheets['description']]
        if any(blank_desc):
            self.warn('Missing worksheet description', 1)

        # and finally store a list of dictionaries of data worksheet summary details
        self.dataworksheet_summaries = [dict(zip(data_worksheets.keys(), vals))
                                        for vals in zip(*data_worksheets.values())]

        # check for extra undocumented spreadsheets
        if 'Worksheet name' in summary:
            expected_sheets = set(data_worksheets['name']) | {'Summary', 'Taxa', 'Locations'}
            if not self.sheet_names.issubset(expected_sheets):
                self.warn('Undocumented sheets found in  workbook: ', 1,
                          join=self.sheet_names - expected_sheets)

        # summary of processing
        if (self.n_warnings - start_warn) > 0:
            self.info('Summary contains {} errors'.format(self.n_warnings - start_warn), 1)
        else:
            self.info('Summary formatted correctly', 1)

    def load_locations(self, locations_json=None):

        """
        Attempts to load and check the contents of the Locations worksheet and
        compile the geographic extent of the locations used.

        Args:
            locations_json: A path to a JSON file containing a valid set of location names.
                With the default value of None, the function tries to get this from
                a SAFE project website service.
        Returns:
            Populates the locations field of the Dataset object
        """

        # try and get the locations worksheet
        self.info("Checking Locations worksheet")
        start_warn = self.n_warnings
        try:
            locs_wb = self.workbook['Locations']
        except KeyError:
            # No locations is pretty implausible, but still persevere as if
            # they aren't going to be required
            self.hint("No locations worksheet found - moving on", 1)
            return

        # GET THE GAZETTEER VALIDATION INFORMATION
        loc_payload = None
        if locations_json is None:
            # If no file is provided then try and get locations from the website service
            loc_get = requests.get('https://www.safeproject.net/call/json/get_locations_bbox')
            if loc_get.status_code != 200:
                self.warn('Could not download locations. Use a local json file.', 1)
            else:
                loc_payload = loc_get.json()
        else:
            # try and load the file
            try:
                loc_payload = simplejson.load(file(locations_json))
            except IOError:
                self.warn('Could not load location names from file.', 1)

        # process the payload
        if loc_payload is not None:
            valid_locations = loc_payload['locations']
            aliases = loc_payload['aliases']
        else:
            # create empty dictionaries to continue with validation
            valid_locations = {}
            aliases = {}

        # ITERATE OVER THE WORKSHEET ROWS
        loc_rows = locs_wb.rows

        # Get the field headers:
        # Duplicated headers are a problem because the values
        # in the locations dictionaries get overwritten.
        hdrs = [cl.value for cl in loc_rows.next()]
        dupes = duplication([h for h in hdrs if not is_blank(h)])
        if dupes:
            self.warn('Duplicated location sheet headers: ', 1, join=dupes)

        # Check location names are available
        if 'Location name' not in hdrs:
            self.warn('Location name column not found', 1)
            self.locations = set()
            return

        # Convert remaining rows into a list of location dictionaries
        locs = [{ky: cl.value for ky, cl in zip(hdrs, rw)} for rw in loc_rows]

        # strip out any rows that consist of nothing but empty cells
        locs = [row for row in locs if not all([is_blank(vl) for vl in row.values()])]

        # Location name cleaning - get names as strings
        for rw in locs:
            rw['Location name'] = unicode(rw['Location name'])

        # check for rogue whitespace
        ws_padded = [rw['Location name'] for rw in locs if is_padded(rw['Location name'])]
        if ws_padded:
            self.warn('Locations names with whitespace padding: ', 1,
                      join=ws_padded, as_repr=True)
            # clean whitespace padding
            for row in locs:
                row['Location name'] = row['Location name'].strip()

        # look for duplicates
        dupes = duplication([rw['Location name'] for rw in locs])
        if dupes:
            self.warn('Duplicated location names: ', 1, join=dupes, as_repr=True)

        # VALIDATE LOCATIONS
        # Get existing location names and aliases -new location names must not appear
        # in it and existing locations must.
        existing_loc_names = set(valid_locations.keys() + aliases.keys())

        # Split up old and new if there are there any new ones?
        if 'New' in hdrs:

            # Check the New column is just yes, no
            is_new = {rw['New'] for rw in locs}
            valid_new = {'yes', 'no', 'Yes', 'No'}
            if not is_new.issubset(valid_new):
                self.warn('New field contains values other than Yes and No: ', 1,
                          join=is_new - valid_new, as_repr=True)

            # extract the new and old locations
            new_locs = [rw for rw in locs if rw['New'].lower() == 'yes']
            locs = [rw for rw in locs if rw['New'].lower() == 'no']
        else:
            new_locs = None

        # Process new locations if there are any
        if new_locs:
            self.hint('{} new locations reported'.format(len(new_locs)), 1)

            # check Lat Long and Type, which automatically updates the extents.
            # Unlike a data worksheet field, here we don't have any metadata or want
            # to keep it, so field checker gets passed an empty dictionary, which is discarded.
            if 'Latitude' in hdrs:
                lats = [vl['Latitude'] for vl in new_locs if vl['Latitude'] != u'NA']
                non_blank_lats = [vl for vl in lats if not is_blank(vl)]
                if len(non_blank_lats) < len(lats):
                    self.warn('Blank latitude values for new locations: use NA.', 2)
                self.check_field_geo({}, non_blank_lats, which='latitude')
            else:
                self.warn('New locations reported but Latitude field missing', 2)

            if 'Longitude' in hdrs:
                longs = [vl['Longitude'] for vl in new_locs if vl['Longitude'] != u'NA']
                non_blank_longs = [vl for vl in longs if not is_blank(vl)]
                if len(non_blank_longs) < len(longs):
                    self.warn('Blank longitude values for new locations: use NA.', 2)
                self.check_field_geo({}, non_blank_longs, which='longitude')
            else:
                self.warn('New locations reported but Longitude field missing', 2)

            if 'Type' in hdrs:
                geo_types = {vl['Type'] for vl in new_locs}
                bad_geo_types = geo_types - {'POINT', 'LINESTRING', 'POLYGON'}
                if bad_geo_types:
                    self.warn('Unknown location types: ', 2, join=bad_geo_types, as_repr=True)
            else:
                self.warn('New locations reported but Type field missing', 2)

            duplicates_existing = [rw['Location name'] for rw in new_locs
                                   if rw['Location name'] in existing_loc_names]

            if duplicates_existing:
                self.warn('New location names duplicate existing names and aliases: ', 2,
                          join=duplicates_existing)

            # new location names
            new_loc_names = {rw['Location name'] for rw in new_locs}
        else:
            new_loc_names = set()

        # Process existing locations if there are any
        if locs:
            # check names exist
            loc_names = {rw['Location name'] for rw in locs}
            unknown = loc_names - existing_loc_names
            if unknown:
                self.warn('Unknown locations found: ', 1, join=unknown, as_repr=True)

            # are aliases being used?
            aliased = loc_names & set(aliases.keys())
            if aliased:
                self.hint('Locations aliases used. Maybe change to primary location names: ',
                          1, join=aliased)

            # Get the bounding box of known locations and aliased locations
            bbox_keys = (loc_names - (unknown | aliased)) | {aliases[ky] for ky in aliased}

            # get the extents of known unaliased locations
            if bbox_keys:
                bbox = [vl for ky, vl in valid_locations.iteritems() if ky in bbox_keys]
                bbox = zip(*bbox)
                self.update_extent((min(bbox[0]), max(bbox[1])), float, 'longitudinal_extent')
                self.update_extent((min(bbox[2]), max(bbox[3])), float, 'latitudinal_extent')
        else:
            loc_names = set()

        # reporting
        if (self.n_warnings - start_warn) > 0:
            self.info('Locations contains {} errors'.format(self.n_warnings - start_warn), 1)
        else:
            self.info('{} locations loaded correctly'.format(len(loc_names)), 1)

        self.locations = loc_names | new_loc_names

    @staticmethod
    def gbif_local_to_index_tuple(row):
        """
        Turns a row from the local GBIF database into a tuple for the taxon index
        Args:
            row: A row from the GBIF backbone database as a dictionary

        Returns:
            A 5 tuple: taxon ID, parent ID, name, rank, status
        """

        if row['parentNameUsageID'] == '':
            parent = 0
        else:
            parent = int(row['parentNameUsageID'])

        return [(int(row['taxonID']), parent, row['canonicalName'],
                row['taxonRank'], row['taxonomicStatus'])]

    def web_gbif_validate(self, tax, rnk, backbone_types):

        """
        Validates a taxon name and rank against the GBIF web API.

        Note that this function only really uses self as to get a reference to the logging,
        which could be reimplemented using the standard logging module, and this could then
        become a simple static method.

        Args:
            tax: The (case insensitive) taxon name
            rnk: The taxon rank
            backbone_types: A list of higher taxon rank names to be added
                to the index along with the target taxon

        Returns:
            If successful, a 2-list of lists of index tuples, separating the
            canonical match for this taxon from the rest of the taxonomic index for it
                [[(index tuple for canonical taxon)],
                 [(index tuples for higher taxa + synonym)]]
            Otherwise, None.
        """

        # use the GBIF API
        url_gbif = ("http://api.gbif.org/v1/species/match"
                    "?name={}&rank={}&strict=true".format(tax, rnk))
        tax_gbif = requests.get(url_gbif)
        if tax_gbif.status_code != 200:
            self.warn('API validation issue for {}'.format(tax))
            return
        else:
            resp = tax_gbif.json()

        # get a set of pairs of rank, parent rank and add a root key
        backbone_pairs = zip(backbone_types, ['root'] + backbone_types[:-1])

        if resp['matchType'] == u'NONE':
            # No match found: catch any notes from GBIF for odd cases (such as two equal
            # matches as in the case of Morus, for example).
            if 'note' in resp.keys():
                rnk += ': ' + resp['note']
            self.warn('No unique match found for {} at rank {}'.format(tax, rnk), 2)
            return
        elif resp['status'].lower() == 'doubtful':
            # We don't accept doubtful names
            self.warn('Taxon {} is considered doubtful'.format(tax), 2)
            return
        elif resp['status'].lower() == 'accepted':
            self.info('Match found for {} at rank {}'.format(tax, rnk), 2)
            # populate the index - from the single response we can populate not
            # only the taxon but also the higher taxa. Before returning, pop the
            # canonical response off and return that separately from the rest of
            # the index entries.
            resp['rootKey'] = 0
            to_add = {rk: (resp[rk + 'Key'], resp[pr + 'Key'], resp[rk], rk, 'accepted')
                      for rk, pr in backbone_pairs if rk in resp}
            canon = [to_add.pop(rnk)]
            return [canon, to_add.values()]
        else:
            # A non-accepted match, so look up the accepted usage to notify the user.
            # Need to get the taxon level key for this rank, which will be the accepted
            # usage key
            acc_key = resp[rnk + 'Key']
            acc_url = "http://api.gbif.org/v1/species/{}".format(acc_key)
            tax_acc = requests.get(acc_url)
            if tax_acc.status_code != 200:
                self.warn('API validation issue for ID {}'.format(acc_key))
                return
            else:
                accpt = tax_acc.json()
                tx_tup = (tax, resp['status'].lower(), accpt['canonicalName'])
                self.hint('{} is considered a {} of {}'.format(*tx_tup), 2)

                # extend the index with full taxonomy from the accepted usage
                accpt['rootKey'] = 0
                to_add = {rk: (accpt[rk + 'Key'], accpt[pr + 'Key'], accpt[rk], rk, 'accepted')
                          for rk, pr in backbone_pairs if rk in accpt}
                canon = [to_add.pop(rnk)]
                # and then add the taxon as used in the dataset.
                parent_key = canon[0][1]
                to_add['user'] = (resp['usageKey'], parent_key, resp['canonicalName'],
                                  resp['rank'].lower(), resp['status'].lower())
                return [canon, to_add.values()]

    def local_gbif_validate(self, tax, rnk, backbone_types):
        """
        Validates a taxon name and rank against the local GBIF database. Because the local
        GBIF database does not contain higher taxon keys in each row, the function returns
        a list of higher taxa to be looked up later. It could call itself recursively but
        many of these entries are likely to duplicated in other taxa, so it is more efficient
        to collate a local set and add the unique entries later.

        Note that this function only really uses self as to get a reference to the logging,
        which could be reimplemented using the standard logging module, and this would then
        become a simple method.

        Args:
            tax: The (case sensitive) taxon name
            rnk: The taxon rank
            backbone_types: A list of higher taxon rank names to be added
                to the index along with the target taxon

        Returns:
            If successful, a 3-list of lists of tuple indices:
                 [[(index tuple for canonical taxon)],
                  [(index tuple non-canon taxon)],
                  [(tuples for higher taxon lookup)]]
            Otherwise, None.
        """

        tax_sql = ("select * from backbone where canonicalName='{}' and "
                   "taxonRank='{}';".format(tax, rnk.lower()))
        tax_gbif = self.gbif_conn.execute(tax_sql).fetchall()
        tax_gbif = [dict(tx) for tx in tax_gbif]

        # Process the hits
        if len(tax_gbif) == 0:
            self.warn('No match found for {} at rank {}'.format(tax, rnk), 2)
            return

        elif len(tax_gbif) == 1:
            # A single matching canonical name - check status
            tax_gbif = tax_gbif[0]

            if tax_gbif['taxonomicStatus'] == 'doubtful':
                # We don't accept doubtful names
                self.warn('Taxon {} is considered doubtful'.format(tax), 2)
                return

            elif tax_gbif['taxonomicStatus'] == 'accepted':
                # accepted taxon
                self.info('Match found for {} at rank {}'.format(tax, rnk), 2)
                # store the higher taxa to lookup fully later
                add_to_higher = [(tax_gbif[hi], hi) for hi in backbone_types
                                 if hi in tax_gbif and tax_gbif[hi] != u'' and hi != rnk]
                # return the accepted usage, an empty list since there isn't an unaccepted
                # usage and the list of higher taxa
                return [self.gbif_local_to_index_tuple(tax_gbif), [], add_to_higher]
            else:
                # synonym or misapplied, so find the accepted usage
                acc_id = tax_gbif['acceptedNameUsageID']
                acc_sql = 'select * from backbone where taxonID = {};'.format(acc_id)
                acc_gbif = self.gbif_conn.execute(acc_sql).fetchall()
                acc_gbif = dict(acc_gbif[0])
                tx_tup = (tax_gbif['canonicalName'], tax_gbif['taxonomicStatus'],
                          acc_gbif['canonicalName'])
                self.hint('{} is considered a {} of {}'.format(*tx_tup), 2)
                # store the higher taxa to lookup fully later
                add_to_higher = [(acc_gbif[hi], hi) for hi in backbone_types
                                 if hi in acc_gbif and acc_gbif[hi] != u'']

                # return the accepted usage, the actual taxon and higher taxa
                return [self.gbif_local_to_index_tuple(acc_gbif),
                        self.gbif_local_to_index_tuple(tax_gbif),
                        add_to_higher]

        else:
            # Multiple hits: reduce the set of hits to give something like the preferred
            # hits reported by GBIF.
            # First, get the taxon statuses
            tx_status = [tx['taxonomicStatus'] for tx in tax_gbif]
            tx_counts = Counter(tx_status)

            if 'accepted' in tx_counts.keys() and tx_counts['accepted'] == 1:
                # only one accepted match alongside other usages
                self.info('Match found for {} at rank {}'.format(tax, rnk), 2)
                # extract that hit
                tax_gbif = tax_gbif[tx_status.index('accepted')]
                # store the higher taxa to lookup fully later
                add_to_higher = [(tax_gbif[hi], hi) for hi in backbone_types
                                 if hi in tax_gbif and tax_gbif[hi] != u'' and hi != rnk]
                # return the accepted usage, an empty list since there isn't an unaccepted
                # usage and the list of higher taxa
                return [self.gbif_local_to_index_tuple(tax_gbif), [], add_to_higher]

            elif 'accepted' in tx_counts.keys() and tx_counts['accepted'] > 1:
                # more than one accepted match
                self.warn('No unique match found for {} at '
                          'rank {}'.format(tax, rnk), 2)
                return

            else:
                # Rows can now contain doubtful (no accepted usage) and then synonyms
                # (of varying kinds) and misapplied (both of which have accepted usage)
                tx_acc = {tx['acceptedNameUsageID'] for tx in tax_gbif
                          if tx['acceptedNameUsageID'] != u''}
                if set(tx_counts.keys()) == {u'doubtful'}:
                    # only doubtful records?
                    self.warn('Taxon {} is considered doubtful'.format(tax), 2)
                    return
                elif len(tx_acc) > 1:
                    # More than one accepted usage
                    self.warn('Taxon {} is not accepted and can refer to more than one '
                              'accepted usage.'.format(tax), 2)
                    return
                else:
                    # A single accepted usage - pick the first row to index
                    tax_gbif = tax_gbif[0]
                    acc_sql = ('select * from backbone where taxonID '
                               '= {};'.format(tx_acc.pop()))
                    acc_gbif = self.gbif_conn.execute(acc_sql).fetchall()
                    acc_gbif = dict(acc_gbif[0])
                    tx_tup = (tax, ', '.join(set(tx_status)), acc_gbif['canonicalName'])
                    self.hint('{} is considered a {} of {}'.format(*tx_tup), 2)
                    # store the higher taxa to lookup fully later
                    add_to_higher = [(acc_gbif[hi], hi) for hi in backbone_types
                                     if hi in acc_gbif and acc_gbif[hi] != u'']

                    # return the accepted usage, the actual taxon and higher taxa
                    return [self.gbif_local_to_index_tuple(tax_gbif),
                            self.gbif_local_to_index_tuple(acc_gbif),
                            add_to_higher]

    def load_taxa(self):

        """
        Attempts to load and check the content of the Taxa worksheet. The
        method checks that:
        i)  all taxa have a taxon name, ascientific  name and a taxon type,
            and a parent name and parent type if the taxon requires another
            taxon to hook it into the hierarchy.
        ii) validate scientific and parent names against the GBIF backbone.

        Returns:
            Updates the class instance to add a list of taxon names to be used
            for matching taxa in data worksheets and a set of taxa to build a
            taxon index. The taxon index set contains tuples of:
                (GBIF ID, Parent ID, canonical name, taxonomic rank, status)
        """

        # try and get the taxon worksheet
        self.info("Checking Taxa worksheet")
        start_warn = self.n_warnings
        try:
            sheet = self.workbook['Taxa']
        except KeyError:
            # This might mean that the study doesn't have any taxa, so return an empty
            # set. If the datasets then contain taxonomic names, it'll fail gracefully.
            self.hint("No taxa worksheet found - assuming no taxa in data for now!", 1)
            return

        # A) SANITIZE INPUTS
        # get and check the headers
        tx_rows = sheet.rows
        hdrs = [cl.value for cl in tx_rows.next()]

        # duplicated headers are a problem in that it will cause values in
        # the taxon dictionaries to be overwritten.
        dupes = duplication([h for h in hdrs if not is_blank(h)])
        if dupes:
            self.warn('Duplicated column headers in Taxa worksheet: ', 1, join=dupes)

        # Load dictionaries of the taxa
        taxa = [{ky: cl.value for ky, cl in zip(hdrs, rw)} for rw in tx_rows]

        # strip out any rows that consist of nothing but empty cells
        taxa = [row for row in taxa if not all([is_blank(vl) for vl in row.values()])]

        # check number of taxa found and standardise names
        if len(taxa) == 0:
            self.info('No taxon rows found'.format(len(taxa)), 1)
            return
        else:
            # Remove values keyed to None (blank columns),
            # convert keys to lower case and update the list of headers
            self.info('Checking {} taxa'.format(len(taxa)), 1)
            _ = [tx.pop(None) for tx in taxa if None in tx]
            taxa = [{ky.lower(): vl for ky, vl in tx.iteritems()} for tx in taxa]
            hdrs = taxa[0].keys()

        # clean the contents of whitespace padding
        for idx, tx in enumerate(taxa):
            for ky, vl in tx.iteritems():
                if is_padded(vl):
                    self.warn('Whitespace padding for {} in row {}:'
                              ' {}'.format(ky, idx + 2, repr(vl)), 2)
                    tx[ky] = vl.strip()

        # check which fields are found
        if 'taxon name' not in hdrs:
            # taxon names are not found so can't continue
            self.warn('No taxon name column found - no further checking', 1)
            self.taxon_names = set()
            return
        else:
            # check the taxon names
            taxon_names = [tx['taxon name'] for tx in taxa if not is_blank(tx['taxon name'])]

            # Any duplication in cleaned names
            dupes = duplication(taxon_names)
            if dupes:
                self.warn('Duplicated taxon names found: ', 2, join=dupes)

            # set the unique taxon names for the dataset instance
            self.taxon_names = set(taxon_names)

            # check to see if the validation headers are present
            if 'scientific name' not in hdrs or 'taxon type' not in hdrs:
                self.warn('Both taxon type and scientific name columns required '
                          'for taxon verification', 2)
                return

        # B) VALIDATE TAXONOMY
        # Initialise a set to record failed validation and the set of taxonomic levels
        # that are available to validate in the backbone. See README.md for details of
        # the validation process

        # taxon and parent types that can be checked against gbif
        backbone_types = ['kingdom', 'phylum', 'order', 'class',
                          'family', 'genus', 'species', 'subspecies']
        # alternative types that are allowed as taxon type, but which need a parent
        alt_types = ['morphospecies', 'functional group']

        # The GBIF API returns all the information needed to add higher taxa into
        # the index but the local DB doesn't provide higher taxon ids. For local
        # validation, higher taxon names and rank are stored locally and then unique
        # entries are added after checking all rows.
        local_higher_taxa = set()
        unvalidated = set()

        for idx, tx in enumerate(taxa):

            # use lowercase for matching
            if tx['taxon type'] is not None:
                tx['taxon type'] = tx['taxon type'].lower()
            if tx['parent type'] is not None:
                tx['parent type'] = tx['parent type'].lower()

            # Sanitize the taxon dictionary contents - three cases of partial
            # or missing information that stop further validation
            if is_blank(tx['taxon name']) or is_blank(tx['scientific name']) \
                    or is_blank(tx['taxon type']):
                # - mandatory three entries cannot be blank
                self.warn('Row {}: blank taxon name, scientific name or '
                          'taxon type'.format(idx + 2), 2)
                continue
            elif (not is_blank(tx['parent type'])) ^ (not is_blank(tx['parent name'])):
                # Logical XOR: incomplete parent taxon information
                self.warn('Incomplete parent taxon information provided for {taxon type}: '
                          '{taxon name}'.format(**tx), 2)
                continue
            elif tx['taxon type'] in alt_types and \
                    (is_blank(tx['parent name']) or is_blank(tx['parent type'])):
                # No parent taxon info is unacceptable for morphospecies and functional groups
                self.warn('Parent taxon information must be provided for {taxon type}: '
                          '{taxon name}'.format(**tx), 2)
                continue

            # Are the taxon and parent ranks to be validated one of the big seven
            # that the GBIF backbone provided or a valid alternative in the case of taxon type?
            if tx['parent type'] is None and tx['taxon type'] not in backbone_types + alt_types:
                self.warn('Taxon type {} is not recognized'.format(tx['taxon type']), 2)
                continue

            if tx['parent type'] is not None and tx['parent type'] not in backbone_types:
                self.warn('Parent type {} is not recognized'.format(tx['parent type']), 2)
                continue

            # The taxon input is now sanitized of obvious problems so can be validated against
            # the GBIF backbone database. We now have a taxon to look up or a taxon with a parent
            # to look up. Either way, we need to add the taxon to the index, along with the
            # taxonomic hierarchy.

            if is_blank(tx['parent type']):
                tx_name = tx['scientific name']
                tx_rank = tx['taxon type']
                add_taxon = False
            else:
                tx_name = tx['parent name']
                tx_rank = tx['parent type']
                add_taxon = True

            if self.gbif_conn is None:
                new_index = self.web_gbif_validate(tx_name, tx_rank, backbone_types)
                if new_index is not None:
                    # add the canonical taxon and other ranks to the index
                    self.taxon_index.update(new_index[0] + new_index[1])
            else:
                new_index = self.local_gbif_validate(tx_name, tx_rank, backbone_types)
                if new_index is not None:
                    self.taxon_index.update(new_index[0] + new_index[1])
                    local_higher_taxa.update(new_index[2])

            if add_taxon:
                # For taxa that have been validated via a parent, add a name used in the file
                # and insert the validated parent taxon ID.
                # For non taxon groups (alt_types), insert the taxon name, otherwise use the
                # provided (unvalidated!) scientific name.
                if tx['taxon type'] in alt_types:
                    alt_tx_name = tx['taxon name']
                else:
                    alt_tx_name = tx['scientific name']
                self.taxon_index.add((-1, new_index[0][0][0], alt_tx_name,
                                     tx['taxon type'], 'user provided'))

        # Look up anything added to the local higher taxon set
        if self.gbif_conn:
            self.info('Adding higher taxa from local database', 1)
            for tax, rnk in local_higher_taxa:
                new_index = self.local_gbif_validate(tax, rnk, backbone_types)
                self.taxon_index.update(new_index[0])

        # report on unvalidated taxa
        if unvalidated:
            self.warn('Entrez taxon validation failed for: ', 1, join=unvalidated)

        if (self.n_warnings - start_warn) > 0:
            self.info('Taxa contains {} errors'.format(self.n_warnings - start_warn), 1)
        else:
            self.info('{} taxa loaded correctly'.format(len(taxa)), 1)

    def load_data_worksheet(self, meta):

        """
        A method to load and check the formatting and content of a
        data worksheet, turning the dictionary for this worksheet
        in the dataset object into a DataWorksheet object.

        Args:
            meta: A reference to the metadata dictionary for this
            worksheet in the datasets attribute, containing the
            worksheet name, title and description from the summary.
        """

        # now start populating with data
        if meta['name'] not in self.sheet_names:
            self.warn('Data worksheet {} not found'.format(meta['name']))
            return
        else:
            self.info('Checking data worksheet {}'.format(meta['name']))
            # Create a dataworksheet to store details: basically
            # just a dictionary with dot notation and defaults.
            dwsh = DataWorksheet(meta)
            start_warn = self.n_warnings

        # get the worksheet and data dimensions
        worksheet = self.workbook[dwsh.name]
        dwsh.max_col = worksheet.max_column
        dwsh.max_row = worksheet.max_row

        # trap completely empty worksheets
        if dwsh.max_row == 1:
            self.warn('Worksheet is empty', 1)
            return

        # get the metadata field names
        # - first search at most the first 10 rows for the 'field_name' descriptor
        #   which shows the end of the metadata and the start of the data
        descriptors = [worksheet.cell(column=1, row=rw).value
                       for rw in range(1, min(10, dwsh.max_row) + 1)]
        if 'field_name' in descriptors:
            dwsh.field_name_row = descriptors.index('field_name') + 1
            dwsh.descriptors = descriptors[:dwsh.field_name_row]
        else:
            self.warn('Cannot parse data: field_name row not found', 1)
            return

        # Neither of the row and col maxima are particularly reliable as Excel can hang on to
        # cell references for previously used cells. We can ignore blank columns easily but
        # we do need to know where to actually stop for finding blank data in rows.
        # So, we explicitly check the row numbers, making them a mandatory part of the setup.

        # - get the values
        row_number_cells = worksheet.get_squared_range(1, dwsh.field_name_row + 1, 1, dwsh.max_row)
        row_numbers = [cl[0].value for cl in row_number_cells]

        # - trim blank or whitespace values from the end and update the max row
        while is_blank(row_numbers[-1]):
            row_numbers.pop()

        dwsh.max_row = len(row_numbers) + dwsh.field_name_row

        # now check the row numbers are numbers and if they are
        # do they start at one and go up by one
        if not all([isinstance(vl, numbers.Number) for vl in row_numbers]):
            self.warn('Non-numeric data found in row numbering', 1)
        else:
            if row_numbers[0] != 1:
                self.warn('Row numbering does not start at 1', 1)

            one_increment = [(vl1 - vl2) == 1 for vl1, vl2 in
                             zip(row_numbers[1:], row_numbers[:-1])]
            if not all(one_increment):
                self.warn('Row numbering does not consistently increment by 1', 1)

        # report on detected size
        self.info('Worksheet contains {} rows and {} columns'.format(dwsh.max_row, dwsh.max_col), 1)

        # get the rows of metadata
        metadata = worksheet.get_squared_range(2, 1, dwsh.max_col, dwsh.field_name_row)
        # extract data
        metadata = [[cell.value for cell in row] for row in metadata]
        # turn into columns for each field
        metadata = zip(*metadata)
        # and hence dictionaries for each field
        metadata = [dict(zip(dwsh.descriptors, fld)) for fld in metadata]

        # insert column index and letter
        for idx, fld in enumerate(metadata):
            fld[u'col_idx'] = idx + 2
            fld[u'column'] = utils.get_column_letter(idx + 2)

        # check field names unique (drop None)
        field_names = [fld['field_name'] for fld in metadata if fld['field_name'] is not None]
        dupes = duplication(field_names)
        if dupes:
            self.warn('Field names duplicated: ', 1, join=dupes)

        # get taxa field names for cross checking observation and trait data
        dwsh.taxa_fields = [fld['field_name'] for fld in metadata if fld['field_type'] == 'Taxa']

        # check the data in each field against the metadata
        for meta in metadata:

            # read the values and check them against the metadata
            data = worksheet.get_squared_range(meta['col_idx'], dwsh.field_name_row + 1,
                                               meta['col_idx'], dwsh.max_row)
            data = [cl[0].value for cl in data]
            self.check_field(dwsh, meta, data)

        # add the new DataWorksheet into the Dataset
        self.dataworksheets.append(dwsh)

        # reporting
        if (self.n_warnings - start_warn) > 0:
            self.info('Dataframe contains {} errors'.format(self.n_warnings - start_warn), 1)
        else:
            self.info('Dataframe formatted correctly', 1)

    def check_field(self, dwsh, meta, data):
        """
        Method to test whether the contents of a field are compatible
        with the metadata for the field. Largely used for checking
        data worksheets, but also by other methods. Acts as a wrapper
        around data type specific subfunctions.

        Args:
            dwsh: A DataWorksheet, used to pass worksheet level information,
             which at the moment is just taxa fields.
            meta: A dictionary of metadata
            data: A list of values
        """

        # Print out column checking header
        if is_blank(meta['field_name']):
            self.info('Checking Column {}'.format(meta['column']), 1)
        else:
            self.info('Checking field {field_name}'.format(**meta), 1)

        # Skip any field with no user provided metadata or data
        blank_data = [is_blank(vl) for vl in data]
        blank_meta = [is_blank(vl) for ky, vl in meta.iteritems()
                      if ky not in ['col_idx', 'column']]
        if all(blank_data) and all(blank_meta):
            return
        elif all(blank_meta):
            self.warn('Field contains no descriptor information but does contain values', 2)
            return

        # try and figure out what else is available
        if is_blank(meta['field_name']):
            self.warn('Field name is blank', 2)

        # check the description
        if is_blank(meta['description']):
            self.warn('Description is missing', 2)

        # filter out missing and blank data, except for comments fields, where
        # blanks are not an error
        if meta['field_type'] != 'Comments':

            # Only NA is acceptable
            na_vals = [vl == u'NA' for vl in data]
            if any(na_vals):
                self.hint('{} / {} values missing'.format(sum(na_vals), len(na_vals)), 2)

            # We won't tolerate blank data
            if any(blank_data):
                self.warn('{} cells are blank or contain only whitespace '
                          'text'.format(sum(blank_data)), 2)

            # Return the values that aren't NA, blank or whitespace only
            na_or_blank = [any(tst) for tst in zip(na_vals, blank_data)]
            data = [dt for dt, nb in zip(data, na_or_blank) if not nb]

        # run consistency checks where needed and trap unknown field types
        if meta['field_type'] == 'Date':
            self.check_field_datetime(meta, data, which='date')
        elif meta['field_type'] == 'Datetime':
            self.check_field_datetime(meta, data, which='datetime')
        elif meta['field_type'] == 'Time':
            self.check_field_time(meta, data)
        elif meta['field_type'] == 'Taxa':
            self.check_field_taxa(data)
        elif meta['field_type'] == 'Location':
            self.check_field_locations(data)
        elif meta['field_type'] in ['Categorical', 'Ordered Categorical']:
            self.check_field_categorical(meta, data)
        elif meta['field_type'] == 'Numeric':
            self.check_field_numeric(meta, data)
        elif meta['field_type'] == 'Abundance':
            self.check_field_abundance(meta, data, dwsh.taxa_fields)
        elif meta['field_type'] in ['Categorical Trait', 'Ordered Categorical Trait']:
            self.check_field_trait(meta, data, dwsh.taxa_fields, which='categorical')
        elif meta['field_type'] == 'Numeric Trait':
            self.check_field_trait(meta, data, dwsh.taxa_fields, which='numeric')
        elif meta['field_type'] in ['Categorical Interaction', 'Ordered Categorical Interaction']:
            self.check_field_interaction(meta, data, dwsh.taxa_fields, which='categorical')
        elif meta['field_type'] == 'Numeric Interaction':
            self.check_field_interaction(meta, data, dwsh.taxa_fields, which='numeric')
        elif meta['field_type'] == 'Latitude':
            self.check_field_geo(meta, data, which='latitude')
        elif meta['field_type'] == 'Longitude':
            self.check_field_geo(meta, data, which='longitude')
        elif meta['field_type'] in ['Replicate', 'ID']:
            # We've looked for missing data, no other constraints.
            pass
        elif meta['field_type'] == 'Comments':
            pass
        elif meta['field_type'] is None:
            self.warn('Field type is empty', 2)
        else:
            self.warn('Unknown field type {field_type}'.format(**meta), 1)

        # extend the dataworksheet fields
        dwsh.fields.append(meta)

    # Helper functions for checking data fields
    def _check_meta(self, meta, descriptor):
        """
        A standardised check to see if a required descriptor is present for
        a field and that it isn't simply empty or whitespace. The function
        reports problems to the Messages instance and returns a boolean
        showing if the checks passed successfully.

        Args:
            meta: A dictionary of field metadata descriptors
            descriptor: The name of the descriptor to check.

        Returns:
            A boolean, with True showing no problems and False showing
            that warnings occurred.
        """

        if descriptor not in meta:
            self.warn('{} descriptor missing'.format(descriptor), 2)
            return False
        elif is_blank(meta[descriptor]):
            self.warn('{} descriptor is blank'.format(descriptor), 2)
            return False
        else:
            return True

    def _check_taxon_meta(self, meta, taxa_fields):

        """
        Checks the taxonomic metadata of abundance and trait fields.
        This is more involved that the simple _check_meta(), because
        of the option to provide taxon_name or taxon_field descriptors

        Args:
            meta: A dictionary of metadata descriptors for the field
            taxa_fields: A list of Taxa fields in this worksheet.

        Returns:
            A boolean, with True showing no problems and False showing
            that warnings occurred.
        """

        # Are taxon_name and taxon_field provided and not blank: note use of
        # 'and' rather than '&' to allow missing descriptors to short cut

        tx_nm_prov = ('taxon_name' in meta) and (not is_blank(meta['taxon_name']))
        tx_fd_prov = ('taxon_field' in meta) and (not is_blank(meta['taxon_field']))

        if tx_nm_prov and tx_fd_prov:
            self.warn('Taxon name and taxon field both provided, use one only', 2)
            return False
        elif tx_nm_prov and meta['taxon_name'] not in self.taxon_names:
            self.warn('Taxon name not found in the Taxa worksheet', 2)
            return False
        elif tx_fd_prov and meta['taxon_field'] not in taxa_fields:
            self.warn("Taxon field not found in this worksheet", 2)
            return False
        elif not tx_nm_prov and not tx_fd_prov:
            self.warn("One of taxon name or taxon field must be provided", 2)
            return False
        else:
            if tx_nm_prov:
                self.taxon_names_used.update([meta['taxon_name']])
            return True

    def _check_interaction_meta(self, meta, taxa_fields):

        """
        Checks the taxonomic metadata of interaction fields.
        This is more involved that the simple _check_meta(), because
        of the option to provide taxon_name or taxon_field descriptors
        describing at least two taxonomic identifiers.

        Args:
            meta: A dictionary of metadata descriptors for the field
            taxa_fields: A list of Taxa fields in this worksheet.

        Returns:
            A boolean, with True showing no problems and False showing
            that warnings occurred.
        """

        # Are interaction_name and/or interaction_field provided and not blank:
        #  note use of 'and' rather than '&' to allow missing descriptors to short cut

        iact_nm_prov = ('interaction_name' in meta) and (not is_blank(meta['interaction_name']))
        iact_fd_prov = ('interaction_field' in meta) and (not is_blank(meta['interaction_field']))

        if not iact_nm_prov and not iact_fd_prov:
            self.warn("At least one of interaction name or interaction field must be provided", 2)
            return False
        else:
            if iact_nm_prov:
                # get the taxon names and descriptions from interaction name providers
                int_nm_lab, int_nm_desc = self._parse_levels(meta['interaction_name'], 2)
                # add names to used taxa
                self.taxon_names_used.update(int_nm_lab)
                # check they are found
                if not all([lab in self.taxon_names for lab in int_nm_lab]):
                    self.warn('Unknown taxa in interaction_name descriptor', 2)
                    nm_check = False
                else:

                    nm_check = True
            else:
                int_nm_lab, int_nm_desc = [(), ()]
                nm_check = True

            if iact_fd_prov:
                # check any field labels match to known taxon fields
                int_fd_lab, int_fd_desc = self._parse_levels(meta['interaction_field'], 2)
                if not all([lab in taxa_fields for lab in int_fd_lab]):
                    self.warn('Unknown taxon fields in interaction_field descriptor', 2)
                    fd_check = False
                else:
                    fd_check = True
            else:
                int_fd_lab, int_fd_desc = [(), ()]
                fd_check = True

            if len(int_nm_lab + int_fd_lab) < 2:
                self.warn('At least two interacting taxon labels or fields must be identified', 2)
                num_check = False
            else:
                num_check = True

            if nm_check and fd_check and num_check:
                return True
            else:
                return False

    def _parse_levels(self, txt, warn_level=1):
        """
        Splits up category information formatted as label:desc;label:desc, which
        is used in both levels for categorical data and interaction descriptors.
        Args:
            txt: The text string to parse
            warn_level: The indent level for any warnings issued

        Returns:
            A list of two tuples of label and descriptions.
        """

        # remove terminal semi-colon, if used.
        if txt.endswith(';'):
            txt = txt[:-1]
        # - split the text up by semi-colon
        parts = txt.split(';')
        # - split descriptions
        parts = [pt.split(':') for pt in parts]
        n_parts = [len(pt) for pt in parts]

        # simple formatting checks
        if any([pt > 2 for pt in n_parts]):
            self.warn('Extra colons in level description.', warn_level)

        # standardise descriptions
        if all([pt == 1 for pt in n_parts]):
            parts = [[pt[0], None] for pt in parts]
        elif all([pt == 2 for pt in n_parts]):
            pass
        else:
            self.warn('Provide descriptions for either all or none of the categories', warn_level)
            parts = [pt[0:2] if len(pt) >= 2 else [pt[0], None] for pt in parts]

        return zip(*parts)

    def check_field_datetime(self, meta, data, which='datetime'):

        """
        Checks for data consistency in date and datetime fields. Excel
        doesn't distinguish, both are loaded as datetime.datetime objects
        so we check that the values are compatible with the user provided
        field type.

        Args:
            meta: The field metadata, to be updated with the range
            data: A list of data values, allegedly containing datetimes or dates
            which: The datetime type to be checked for.
        """

        # Check type (excluding NA values)
        is_datetime = [isinstance(dt, datetime.datetime) for dt in data]
        if not all(is_datetime):
            self.warn('Non-date data in field.', 2)
            is_time = [isinstance(dt, datetime.time) for dt in data]
            if any(is_time):
                self.hint('Some values _only_  contain time components', 2)

        # Check no time component in actual dates
        if which == 'date':
            no_time = [vl.time() == datetime.time(0, 0) for vl in data]
            if not all(no_time):
                self.warn('Some values also contain time components', 2)

        # update the field metadata and the dataset extent
        real_dates = [dt for dt in data if isinstance(dt, datetime.datetime)]
        extent = (min(real_dates), max(real_dates))
        meta['range'] = extent
        self.update_extent(extent, datetime.datetime, 'temporal_extent')

    def check_field_time(self, meta, data):

        """
        Checks for data consistency in time fields and reports to the
        Messages instance.

        Args:
            meta: The field metadata, to be updated with the range
            data: A list of data values, allegedly of type datetime.time
        """

        # Check type (excluding NA values)
        type_check = set([type(vl) for vl in data])
        if type_check != {datetime.time}:
            self.warn('Non-time formatted data found.', 2)

        # update the field metadata
        real_times = [dt for dt in data if isinstance(dt, datetime.time)]
        meta['range'] = (min(real_times), max(real_times))

    def check_field_taxa(self, data):

        """
        Checks if all the values provided in a Taxon field are found
        in the Taxa worksheet.

        Args:
            data: A list of data values, allegedly taxon names
        """

        found = set(data)
        if self.taxon_names == set():
            self.warn('No taxa loaded', 2)
        if not found.issubset(self.taxon_names):
            self.warn('Includes taxa missing from Taxa worksheet: ', 2,
                      join=found - self.taxon_names, as_repr=True)

        # add the found taxa to the list of taxa used
        self.taxon_names_used.update(found)

    def check_field_locations(self, data):

        """
        Checks if all the values provided in a Locations field are
        found in the Locations worksheet, reporting to the Messages instance.

        Args:
            data: A list of data values, allegedly taxon names
        """

        # location names should be strings
        data = [unicode(dt) for dt in data]

        # check if locations are all provided
        found = set(data)
        if self.locations == set():
            self.warn('No locations loaded', 2)
        elif not found.issubset(self.locations):
            self.warn('Includes locations missing from Locations worksheet:', 2,
                      join=found - self.locations, as_repr=True)

        # add the locations to the set of locations used
        self.locations_used.update(found)

    def check_field_abundance(self, meta, data, taxa_fields):

        """
        Checks abundance type data, reporting to the Messages instance.

        Args:
            meta: A dictionary of metadata descriptors for the field
            data: A list of data values, allegedly numeric
            taxa_fields: A list of Taxa fields in this worksheet.
        """

        # check the required descriptors
        self._check_meta(meta, 'method')
        self._check_taxon_meta(meta, taxa_fields)

        # Can still check values are numeric, whatever happens above.
        # We're not going to insist on integers here - could be mean counts.
        all_nums, nums = all_numeric(data)
        if not all_nums:
            self.warn('Field contains non-numeric data', 2)

        # update the field metadata
        meta['range'] = (min(nums), max(nums))

    def check_field_categorical(self, meta, data):

        """
        Checks factor data, reporting to the Messages instance.

        Args:
            meta: A dictionary of metadata descriptors for the field
            data: A list of data values, allegedly strings from a provided set of levels
        """

        # Are levels found?
        ct_ok = self._check_meta(meta, 'levels')

        if not ct_ok:
            # Can't check further if no levels descriptor
            pass
        elif ct_ok and not isinstance(meta['levels'], unicode):
            # Can't really check anything here either
            self.warn('Category description does not seem to be text', 2)
        else:
            # Now we can test if the labels match up
            level_labels, level_desc = self._parse_levels(meta['levels'], 2)

            # - repeated labels?
            if len(set(level_labels)) < len(level_labels):
                self.warn('Repeated level labels', 2)

            # - check for integer level names
            integer_codes = [is_integer_string(vl) for vl in level_labels]
            if any(integer_codes):
                self.warn('Integer level names not permitted', 2)

            # Now look for consistency: get the unique values reported in the
            # data, convert to unicode to handle checking of integer labels and
            # then check the reported levels are a subset of the descriptors.
            reported = set(data)
            reported = {unicode(lv) for lv in reported}

            if not reported.issubset(level_labels):
                self.warn('Categories found in data missing from '
                          'description: ', 2, join=reported - set(level_labels))

    def check_field_numeric(self, meta, data):

        """
        Checks numeric type data, reporting to the Messages instance.

        Args:
            meta: A dictionary of metadata descriptors for the field
            data: A list of data values, allegedly numeric
        """

        # Check required descriptors
        self._check_meta(meta, 'units')
        self._check_meta(meta, 'method')

        # Regardless of the outcome of the meta checks, can still check the
        # data is all numeric, as it claims to be.
        all_nums, nums = all_numeric(data)
        if not all_nums:
            self.warn('Field contains non-numeric data', 2)

        # update the field metadata
        meta['range'] = (min(nums), max(nums))

    def check_field_trait(self, meta, data, taxa_fields, which='categorical'):

        """
        Checks trait type data and reports to the Messages instance. Just a wrapper
        to check that a valid taxon has been provided before handing off to
        check_field_categorical or check_field_numeric

        Args:
            meta: A dictionary of metadata descriptors for the field
            data: A list of data values, allegedly numeric
            taxa_fields: A list of Taxa fields in this worksheet.
            which: The type of trait data in the field
        """

        # Check required descriptors
        self._check_taxon_meta(meta, taxa_fields)

        # Regardless of the outcome of the meta checks, check the contents:
        if which == 'categorical':
            self.check_field_categorical(meta, data)
        elif which == 'numeric':
            self.check_field_numeric(meta, data)

    def check_field_interaction(self, meta, data, taxa_fields, which='categorical'):

        """
        Checks interaction type data and reports to the Messages instance.
        Just a wrapper to check that intercating taxa have been identified
        before handing off to check_field_categorical or check_field_numeric

        Args:
            meta: A dictionary of metadata descriptors for the field
            data: A list of data values, allegedly numeric
            taxa_fields: A list of Taxa fields in this worksheet.
            which: The type of trait data in the field
        """

        # Check required descriptors
        self._check_interaction_meta(meta, taxa_fields)

        # Regardless of the outcome of the meta checks, check the contents:
        if which == 'categorical':
            self.check_field_categorical(meta, data)
        elif which == 'numeric':
            self.check_field_numeric(meta, data)

    def check_field_geo(self, meta, data, which='latitude'):

        """
        Checks geographic coordinates. It also automatically updates
        the geographic extent of the dataset.

        Args:
            meta: A dictionary of metadata descriptors for the field
            data: A list of data values
            which: One of latitude or longitude
        """

        # Are the values represented as decimal degrees - numeric.
        all_nums, nums = all_numeric(data)
        if not all_nums:
            self.warn('Non numeric data found', 2)
            if any([RE_DMS.search(unicode(vl)) for vl in data]):
                self.hint('Possible degrees minutes and seconds formatting? Use decimal degrees', 2)

        # Check the locations
        if len(nums):
            min_geo = float(min(nums))
            max_geo = float(max(nums))
            extent = (min_geo, max_geo)

            if which == 'latitude':
                bnds = [-90, -4, 8, 90]
            elif which == 'longitude':
                bnds = [-180, 108, 120, 180]

            out_of_bounds = min_geo < bnds[0] or max_geo > bnds[3]
            out_of_borneo = bnds[0] <= min_geo < bnds[1] or bnds[2] < max_geo <= bnds[3]

            if out_of_bounds:
                self.warn('{0} values not in valid range[{1[0]}, {1[3]}]: '
                          '{2}'.format(which.capitalize(), bnds, extent), 2)
            elif out_of_borneo:
                self.hint('{0} values not in Borneo [{1[1]}, {1[2]}]: '
                          '{2}'.format(which.capitalize(), bnds, extent), 2)

            # update the field metadata and the dataset extent
            meta['range'] = extent
            # Look up the extent name to update and then update it
            which_extent = {'latitude': 'latitudinal_extent', 'longitude': 'longitudinal_extent'}
            self.update_extent(extent, float, which_extent[which])

    def check_locations_and_taxa_used(self):
        """
        A method to check that the locations and taxa provided have actually been used
        in the data worksheets scanned.
        """
        self.info('Checking provided locations and taxa all used in data worksheets')

        # check locations
        if self.locations_used != self.locations:
            self.warn('Provided locations not used: ', 1,
                      join=self.locations - self.locations_used)
        else:
            self.info('Provided locations all used in datasets', 1)

        # check taxa
        if self.taxon_names_used != self.taxon_names:
            self.warn('Provided taxa  not used: ', 1,
                      join=self.taxon_names - self.taxon_names_used)
        else:
            self.info('Provided taxa all used in datasets', 1)

    def export_metadata_dict(self):
        """
        Function to return a simple dictionary of the metadata content, combining
        dataworksheets, into a single representation of the contents

        Returns:
            A dictionary of metadata
        """

        # get the required components
        component_keys = ['access', 'authors', 'description', 'embargo_date', 'filename',
                          'keywords', 'latitudinal_extent', 'longitudinal_extent',
                          'project_id', 'temporal_extent', 'title']

        components = {ky: vl for ky, vl in self.__dict__.iteritems() if ky in component_keys}

        # add in the data worksheets
        components['dataworksheets'] = [dwsh.__dict__ for dwsh in self.dataworksheets]

        return components


# High level functions
def check_file(fname, verbose=True, gbif_database=None, locations_json=None, validate_doi=False):

    """
    Runs the format checking across an Excel workbook.

    Parameters:
        fname: Path to an Excel file
        verbose: Boolean to indicate whether to print messages as the program runs?
        gbif_database: Path to a local copy of the GBIF taxonomy backbone if that is to be used
        instead of the GBIF web API.
        locations_json: The path to a json file of valid location names
        validate_doi: Check any publication DOIs resolve, requiring a web connection.

    Returns:
        A Dataset object
    """

    # initialise the dataset object
    dataset = Dataset(fname, verbose=verbose, gbif_database=gbif_database)

    # load the metadata sheets
    dataset.load_summary(validate_doi=validate_doi)
    dataset.load_taxa()
    dataset.load_locations(locations_json=locations_json)

    # check the datasets
    if dataset.dataworksheet_summaries:
        for ws in dataset.dataworksheet_summaries:
            dataset.load_data_worksheet(ws)
    else:
        dataset.warn('No data worksheets found')

    # finally check if the provide locations and taxa are all used
    dataset.check_locations_and_taxa_used()

    if dataset.n_warnings:
        dataset.info('FAIL: file contained {} errors'.format(dataset.n_warnings))
    else:
        dataset.info('PASS: file formatted correctly')

    return dataset


def main():

    """
    This program validates an Excel file formatted as a SAFE dataset. As it runs, it outputs
    a report that highlights any problems with the formatting.

    Much of the validation is to check that the data meets our metadata standards and is
    internally consistent. However, it uses external sources to perform validation in three
    areas.

    1. Taxon validation. The program validates taxonomic names against the GBIF taxonomy
    backbone. By default, it uses the GBIF web API to validate names, but can also use a
    local copy of the backbone provided in a sqlite database: this will work offline and
    is much faster but requires some simple setup.

    2. Location names. The program also validate sampling location names against the SAFE gazeteer.
    By default, this is loaded automatically from the SAFE website so requires an internet
    connection, but a local copy can be provided for offline use.

    3. DOI checking. Optionally, the program will validate any DOIs provided as having used
    the database. This requires a web connection and cannot be performed offline.
    """

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('fname', help="Path to the Excel file to be validated.")
    parser.add_argument('-l', '--locations_json', default=None,
                        help='Path to a locally stored json file of valid location names')
    parser.add_argument('-g', '--gbif_database', default=None,
                        help=('The path to a local sqlite database containing the GBIF '
                              'taxonomy backbone.'))
    parser.add_argument('--validate_doi', action="store_true", default=False,
                        help=('Check the validity of any publication DOIs, '
                              'provided by the user. Requires a web connection.'))

    args = parser.parse_args()

    check_file(fname=args.fname, verbose=True, locations_json=args.locations_json,
               gbif_database=args.gbif_database, validate_doi=args.validate_doi)

if __name__ == "__main__":
    main()

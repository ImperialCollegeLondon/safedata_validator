#!/usr/bin/python

from __future__ import print_function
import openpyxl
from openpyxl import utils
import sys
import datetime
import re
from collections import Counter
from StringIO import StringIO
import numbers

"""
Module containing code to verify the format of a SAFE project Excel dataset.

The functions are written to extract as much information as possible from
checking a file so, rather than raising an error and halting at the first
fault, functions are written to check as much as they can on the inputs.
Any information and warnings are added on to a Messages instance, passed
to each function, which is used to keep a report of issues throughout the
file check.
"""

# define some regular expressions used to check validity
re_orcid = re.compile('[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{4}')
re_email = re.compile('\S+@\S+\.\S+')
re_name = re.compile('[^,]+,[ ]?[^,]+')
re_whitespace_only = re.compile('^\s*$')


class Messages(object):

    """
    This class keeps a record of messages to be conveyed to the
    user and keeps a tally of warnings.

    Attributes:
        header: A header for the messages report.
        verbose: If True, messages are printed to the screen as they are added.
        kinds: A dictionary keyed by available message kinds, with values as string
            prefixes for messages of that kind.
        messages: A list of tuples of messages giving (kind, text, level).
        new_warnings: A count of the warning messages added since the instance was
            created or since the get_new_warnings_count() method was last called.
    """

    def __init__(self, header, verbose):
        self.header = header
        self.messages = []
        self.verbose = verbose
        self.kinds = {'info': '-', 'hint': '?', 'warn': '!'}
        self.new_warnings = 0
        if verbose:
            print(self.header)

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

    def warn(self, message, level=0):
        """
        Adds a warning message to the instance.
        Args:
            message: The warning message text.
            level: The message indent level
        """
        msg = ('warn', message, level)
        self.messages.append(msg)
        self.new_warnings += 1
        if self.verbose:
            self.print_msg(*msg)

    def hint(self, message, level=0):
        """
        Adds an hint message to the instance. This is used to indicate where
        something might need to be checked but isn't necessarily an error
        Args:
            message: The hint text.
            level: The hint indent level
        """
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
        report.writelines([self.header])
        msgs = ['  ' * lv + self.kinds[kn] + msg for lv, kn, msg in self.messages]
        report.writelines(msgs)
        return report

    # Functions to check warnings
    def get_new_warning_count(self):
        """
        Reports the number of warnings since the instance was initialised or
        since this function was last called.
        Returns:
            An integer
        """
        if self.new_warnings:
            new = self.new_warnings
            self.new_warnings = 0
            return new
        else:
            return 0

    def count_warnings(self):
        """
        Calculates the total number of warnings in the instances
        Returns:
            An integer
        """
        return sum([msg[0] == 'warn' for msg in self.messages])


def get_summary(wb, m):

    """
    Checks the information in the summary worksheet and looks for the metadata and
    dataset worksheets. The function is intended to try and get as much information
    as possible from the worksheet: the dictionary of metadata returned will have
    None for any missing data, which should be handled by downstream code.

    Args:
        wb: An openpyxl Workbook instance
        m: A Messages instance
    Returns:
        A dictionary of the available summary metadata
    """

    # try and get the summary worksheet
    m.info("Checking Summary worksheet")
    try:
        ws = wb['Summary']
    except KeyError:
        m.warn("Summary worksheet not found, moving on.")
        return None

    # load dictionary of summary information block, allowing for multiple
    # columns for fields (data compilers, dataset sheets).
    summary_dict = {}
    col_range = range(2, ws.max_column + 1)
    max_row = ws.max_row

    for rw in range(1, max_row + 1):
        vals = [ws.cell(row=rw, column=cl).value for cl in col_range]
        summary_dict[ws.cell(row=rw, column=1).value] = vals

    # build return dictionary
    ret_dict = {}

    # Check the minimal keys are expected
    required = {"SAFE Project ID", "Access status", "Title", "Description",
                "Author name", "Author email", "Author affiliation", "Author ORCID",
                "Worksheet name", "Worksheet title", "Worksheet description"}

    found = set(summary_dict.keys())

    # don't bail here - try and get as far as possible
    if not found.issuperset(required):
        m.warn('Missing metadata fields: {}'.format(','.join(required - found)), 1)

    # check for contents
    if any([set(x) == {None} for x in summary_dict.values()]):
        m.warn('Metadata fields with no information.', 1)

    # CHECK PROJECT ID
    if 'SAFE Project ID' not in summary_dict:
        m.warn('SAFE Project ID missing', 1)
        ret_dict['project_id'] = None
    else:
        if type(summary_dict['SAFE Project ID'][0]) != long:
            m.warn('SAFE Project ID is not an integer.', 1)
        ret_dict['project_id'] = summary_dict['SAFE Project ID'][0]

    # CHECK DATASET TITLE
    if 'Title' not in summary_dict:
        m.warn('Dataset title missing', 1)
        ret_dict['title'] = None
    else:
        ret_dict['title'] = summary_dict['Title'][0]

    # CHECK DATASET DESCRIPTION
    if 'Description' not in summary_dict:
        m.warn('Dataset title missing', 1)
        ret_dict['description'] = None
    else:
        ret_dict['description'] = summary_dict['Title'][0]

    # CHECK ACCESS STATUS AND EMBARGO DETAILS
    if 'Access status' in summary_dict:
        access_status = summary_dict['Access status'][0]
        ret_dict['access_status'] = summary_dict['Access status'][0]
        ret_dict['embargo_date'] = None
        if access_status == 'Embargo':
            if 'Embargo date' not in summary_dict:
                m.warn('Dataset embargoed but no date provided.', 1)
            embargo_date = summary_dict['Embargo date'][0]
            now = datetime.datetime.now()
            if type(embargo_date) != datetime.datetime:
                m.warn('Embargo date not formatted as date.', 1)
            elif embargo_date < now:
                m.warn('Embargo date is in the past.', 1)
            elif embargo_date > now + datetime.timedelta(days=2*365):
                m.warn('Embargo date more than two years in the future.', 1)
            else:
                ret_dict['embargo_date'] = embargo_date
        elif access_status == 'Open':
            pass
        else:
            m.warn('Access status must be Open or Embargo not {}'.format(access_status), 1)

    # CHECK AUTHORS
    author_keys = ['Author name', 'Author affiliation', 'Author email', 'Author ORCID']
    if set(summary_dict.keys()).issuperset(author_keys):
        authors = zip(*[summary_dict[x] for x in author_keys])

        # remove any completely blank entries
        authors = [x for x in authors if x != tuple([None] * 4)]

        # convert to dict in Zenodo style and check completeness and validity
        for ind, au in enumerate(authors):
            au = {k: v for k, v in zip(['name', 'affiliation', 'email', 'orcid'], au)}
            # look for missing values
            for k, v in au.iteritems():
                if v is None:
                    m.warn('Author {} missing'.format(k), 1)
            # check validity
            if not re_name.match(au['name']):
                m.warn('Author name not formated as last_name, '
                       'first_names: {}'.format(au['name']), 1)
            if not re_orcid.match(str(au['orcid'])) or au['orcid'] is None:
                m.warn('ORCID not properly formatted: {}'.format(au['orcid']), 1)
            if not re_email.match(au['email']):
                m.warn('Email not properly formatted: {}'.format(au['email']), 1)
            authors[ind] = au

        ret_dict['authors'] = authors
    else:
        m.warn('Author metadata block incomplete.', 1)
        ret_dict['authors'] = None

    # CHECK DATASETS
    ds_keys = ['Worksheet name', 'Worksheet title', 'Worksheet description']
    valid_sheets = set(wb.get_sheet_names())

    if set(summary_dict.keys()).issuperset(ds_keys):
        datasets = zip(*[summary_dict[x] for x in ds_keys])

        # remove any completely blank entries
        datasets = [x for x in datasets if x != tuple([None] * 4)]

        # validate
        for ind, dt in enumerate(datasets):
            dt = {k: v for k, v in zip(['worksheet', 'title', 'description'], dt)}
            # look for missing values
            for k, v in dt.iteritems():
                if v is None:
                    m.warn('Dataset {} missing'.format(k), 1)
            # validate
            if dt['worksheet'] not in valid_sheets:
                m.warn('Dataset worksheet not found: {}.'.format(dt['worksheet']), 1)
            datasets[ind] = dt

        ret_dict['datasets'] = datasets
    else:
        m.warn('Dataset metadata block incomplete.', 1)
        ret_dict['datasets'] = None

    # check for extra undocumented spreadsheets
    if 'Worksheet name' in summary_dict:
        expected_sheets = set(summary_dict['Worksheet name'] +
                              ['Summary', 'Taxa', 'Locations'])
        if valid_sheets != expected_sheets:
            m.warn('Extra sheets found in workbook: '
                   '{}'.format(','.join(valid_sheets - expected_sheets)), 1)

    new_warn = m.get_new_warning_count()
    if new_warn:
        m.info('Summary contains {} errors'.format(new_warn), 1)
    else:
        m.info('Summary formatted correctly', 1)

    return ret_dict


def get_locations(wb, m):

    """
    Attempts to load and check the contents of the Locations worksheet.
    Args:
        wb: An openpyxl Workbook instance
        m: A Messages instance
    Returns:
        A set of provided location names or None if the worksheet cannot
        be found or no locations can be loaded.
    """

    # try and get the locations worksheet
    m.info("Checking Locations worksheet")
    try:
        locs = wb['Locations']
    except KeyError:
        # No locations is pretty implausible, but still persevere as if
        # they aren't going to be required
        m.warn("No locations worksheet found - moving on", 1)
        return None

    # Check the headers
    fields = []
    max_col = locs.max_column
    max_row = locs.max_row

    for cl in range(1, max_col):
        fields.append(locs.cell(row=1, column=cl).value)

    # check the key fields are there
    if not set(fields).issuperset({'Location name'}):
        m.warn('Location name column not found', 1)
        loc_names = None
    else:
        # get the location names
        names_col = fields.index('Location name') + 1
        rows = range(2, max_row + 1)
        loc_names = [locs.cell(row=rw, column=names_col).value for rw in rows]

        # check for duplicate names
        if len(set(loc_names)) != len(loc_names):
            m.warn('Duplicated location names', 1)

        loc_names = set(loc_names)

    new_warn = m.get_new_warning_count()
    if new_warn:
        m.info('Locations contains {} errors'.format(new_warn), 1)
    else:
        m.info('{} locations loaded correctly'.format(len(loc_names)), 1)

    return loc_names


def get_taxa(wb, m):

    """
    Attempts to load and check the content of the Taxa worksheet.
    Args:
        wb: An openpyxl Workbook instance.
        m: A Messages instance.
    Returns:
        A set of provided taxon names or None if the worksheet cannot
        be found or no taxa can be loaded.
    """

    # try and get the taxon worksheet
    m.info("Checking Taxa worksheet")
    try:
        taxa = wb['Taxa']
    except KeyError:
        # This might mean that the study doesn't have any taxa, so return an empty
        # set. If the datasets then contain taxonomic names, it'll fail gracefully.
        m.hint("No taxa worksheet found - assuming no taxa in data for now!", 1)
        return set()

    # Check the headers
    fields = []
    max_col = taxa.max_column
    max_row = taxa.max_row

    for cl in range(1, max_col):
        fields.append(taxa.cell(row=1, column=cl).value)

    # check the two key fields are there
    if not set(fields).issuperset({'Taxon name', 'Taxon type'}):
        m.warn('One or both of Taxon name and Taxon type columns not found', 1)
        return None

    # get the taxon names and types
    names_col = fields.index('Taxon name') + 1
    types_col = fields.index('Taxon type') + 1
    rows = range(2, max_row + 1)
    taxon_names = [taxa.cell(row=rw, column=names_col).value for rw in rows]
    taxon_types = [taxa.cell(row=rw, column=types_col).value for rw in rows]

    # check for duplicate names
    if len(set(taxon_names)) != len(taxon_names):
        m.warn('Duplicated taxon names', 1)

    # Check for some common taxonomic levels. We want as much taxonomic
    # context as possible, but if all taxa are morphospecies / functional
    # groups to family level then genus and species are just empty columns.
    taxon_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    if len(set(taxon_types) & set(taxon_levels)) == 0:
        m.warn('Please provide fields giving at least some of the core taxonomic levels', 1)

    # now cross check the given types against the field headers
    types_found = set(taxon_types)
    types_valid = set(fields) - {'Taxon name', 'Taxon type'}
    if not types_found.issubset(types_valid):
        m.warn('Some rows have taxon types that do not match a column name: '
               '{}'.format(','.join(types_found - types_valid)), 1)

    new_warn = m.get_new_warning_count()
    if new_warn:
        m.info('Taxa contains {} errors'.format(new_warn), 1)
    else:
        m.info('{} taxa loaded correctly'.format(len(taxon_names)), 1)

    return set(taxon_names)


def check_data_worksheet(wb, ws_name, taxa, locations, m):

    """
    Attempt to load and checks the formatting and content of a data worksheet,
    updating the Messages instance `m` with the results.

    Args:
        wb: An openpyxl Workbook instance
        ws_name: The worksheet name
        taxa: A list of valid taxa
        locations: A list of valid locations
        m: A Messages instance
    """

    if ws_name not in wb.get_sheet_names():
        m.warn('Data worksheet {} not found'.format(ws_name))
        return None
    else:
        m.info('Checking data worksheet {}'.format(ws_name))

    # get the worksheet and data dimensions
    ws = wb[ws_name]
    max_col = ws.max_column
    max_row = ws.max_row

    # trap completely empty worksheets
    if max_row == 1:
        m.warn('Worksheet is empty', 2)
        return None

    # get the metadata field names
    # - first search at most the first 20 rows for the 'field_name' descriptor
    #   which shows the end of the metadata and the start of the data
    descriptors = [ws.cell(column=1, row=rw).value for rw in range(1, min(20, max_row) + 1)]
    if 'field_name' in descriptors:
        field_name_row = descriptors.index('field_name') + 1
        descriptors = descriptors[:field_name_row]
    else:
        m.warn('Cannot parse data: field_name row not found', 1)
        return None

    # get the metadata for each field
    field_metadata = []
    for cl in range(2, max_col + 1):
        this_field = {k: ws.cell(column=cl, row=i + 1).value for i, k in enumerate(descriptors)}
        field_metadata.append(this_field)

    # check required descriptors are present and if locations and taxa
    # turn out to be needed after all!
    ft_found = [fld['field_type'] for fld in field_metadata]
    if 'Categorical' in ft_found and 'levels' not in descriptors:
        m.warn('Categorical data fields found but no levels descriptor provided.', 1)

    if 'Numeric' in ft_found:
        if 'units' not in descriptors:
            m.warn('Numeric data fields found but no units descriptor provided.', 1)
        if 'method' not in descriptors:
            m.warn('Numeric data fields found but no units descriptor provided.', 1)

    if 'Location' in ft_found and locations is None:
        m.warn('Location field found but no Location worksheet provided.', 1)

    if 'Abundance' in ft_found:
        if 'taxon_name' not in descriptors:
            m.warn('Abundance field found but no taxon name descriptor provided.', 1)
        if 'method' not in descriptors:
            m.warn('Abundance field found but no sampling method descriptor provided.', 1)

    # check field names unique (drop None)
    fn_found = Counter([fld['field_name'] for fld in field_metadata
                        if fld['field_name'] is not None])
    duplicated = [k for k, v in fn_found.iteritems() if v > 1]
    if len(duplicated) > 0:
        m.warn('Field names duplicated: {}'.format(', '.join(duplicated)), 1)

    # get taxa field names for cross checking observation and trait data
    taxa_fields = [fld['field_name'] for fld in field_metadata if fld['field_type'] == 'Taxa']

    # TODO - check mandatory fields

    # check each field
    for idx, meta in enumerate(field_metadata):

        # skip any field with no metadata
        if set(meta.values()) == {None}:
            break

        # prep the messages instance to pass to functions
        if meta['field_name'] is None or re_whitespace_only.match(meta['field_name']):
            m.info('Checking Column {}'.format(openpyxl.utils.get_column_letter(idx + 2)), 1)
            m.warn('Field name is blank', 2)
        else:
            m.info('Checking field {field_name}'.format(**meta), 1)

        # check the description
        if meta['description'] is None or re_whitespace_only.match(meta['description']):
            m.warn('Description is missing', 2)

        # read the values
        data_block = ws.get_squared_range(idx + 2, field_name_row + 1, idx + 2, max_row)
        data = [cl[0].value for cl in data_block]

        # filter out missing and blank data, except for comments fields, where
        # blanks are not an error
        if meta['field_type'] != 'Comments':
            data = filter_missing_or_blank_data(data, m)

        # run consistency checks where needed and trap unknown field types
        if meta['field_type'] == 'Date':
            check_field_date(data, m)
        elif meta['field_type'] == 'Datetime':
            check_field_datetime(data, m)
        elif meta['field_type'] == 'Time':
            check_field_time(data, m)
        elif meta['field_type'] == 'Taxa':
            check_field_taxa(data, taxa, m)
        elif meta['field_type'] == 'Location':
            check_field_locations(data, locations, m)
        elif meta['field_type'] == 'Categorical':
            check_field_categorical(meta, data, m)
        elif meta['field_type'] == 'Numeric':
            check_field_numeric(meta, data, m)
        elif meta['field_type'] == 'Abundance':
            check_field_abundance(meta, data, taxa, taxa_fields, m)
        elif meta['field_type'] == 'Trait':
            check_field_trait(meta, data, taxa, taxa_fields, m)
        elif meta['field_type'] == 'Comments':
            pass
        elif meta['field_type'] is None:
            m.warn('Field type is empty', 2)
        else:
            m.warn('Unknown field type {field_type}'.format(**meta), 2)

    new_warn = m.get_new_warning_count()
    if new_warn:
        m.info('Dataframe contains {} errors'.format(new_warn), 1)
    else:
        m.info('Dataframe formatted correctly', 1)

    return None


def filter_missing_or_blank_data(data, m):

    """
    Takes a list of data and filters out any missing or blank data,
    reporting to the Messages instance as it goes. The filtered list
    is returned for feeding to field type checker functions so that
    they can assume no NA or blank data in their checks.

    Args:
        data: A list of values read from the Worksheet.
        m: A Messages instance
    Returns:
        A list of data values filtered to remove NA and blanks.
    """

    # Only NA is acceptable
    na_vals = [vl == u'NA' for vl in data]
    if any(na_vals):
        m.hint('{} / {} values missing'.format(sum(na_vals), len(na_vals)), 2)

    # We won't tolerate:
    # 1) empty cells (just to avoid ambiguity - e.g. in abundance data)
    is_empty = [vl is None for vl in data]
    if sum(is_empty):
        m.warn('{} cells are blank'.format(sum(is_empty)), 2)
    # 2) non-empty cells containing only whitespace strings
    ws_only = [re_whitespace_only.match(unicode(vl)) is not None for vl in data]
    if any(ws_only):
        m.warn('{} cells contain whitespace only text'.format(sum(ws_only)), 2)

    # Return the values that aren't NA, blank or whitespace only
    na_or_blank = [any(tst) for tst in zip(na_vals, is_empty, ws_only)]
    data = [dt for dt, nb in zip(data, na_or_blank) if not nb]

    return data


"""
Field checker functions and a couple of helper functions below.
"""


def _check_meta(meta, descriptor, m):
    """
    A standardised check to see if a required descriptor is present for
    a field and that it isn't simply empty or whitespace. The function
    reports problems to the Messages instance and returns a boolean
    showing if the checks passed successfully.

    Args:
        meta: A dictionary of field metadata descriptors
        descriptor: The name of the descriptor to check.
        m: An instance of class Messages

    Returns:
        A boolean, with True showing no problems and False showing
        that warnings occurred.
    """

    if descriptor not in meta:
        m.warn('{} descriptor missing'.format(descriptor), 2)
        return False
    elif meta[descriptor] is None or re_whitespace_only.match(meta[descriptor]):
        m.warn('{} descriptor is blank'.format(descriptor), 2)
        return False
    else:
        return True


def is_integer_string(s):
    """
    Checks if a string value can represent an integer.
    Args:
        s: A string

    Returns:
        A boolean.
    """
    try:
        int(s)
        return True
    except ValueError:
        return False


def check_field_date(data, m):

    """
    Checks for data consistency in date fields and reports to the
    Messages instance.

    Args:
        data: A list of data values, allegedly of type datetime.date
        m: A Messages instance
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data])
    if type_check != {datetime.datetime}:
        m.warn('Non-date data in field.', 2)
        if datetime.time in type_check:
            m.hint('Some values _only_  contain time components', 2)

    # Check no time component in actual dates
    no_time = [vl.time() == datetime.time(0, 0) for vl in data]
    if not all(no_time):
        m.warn('Some values also contain time components', 2)


def check_field_datetime(data, m):

    """
    Checks for data consistency in datetime fields and reports to the
    Messages instance.

    Args:
        data: A list of data values, allegedly of type datetime.datetime
        m: A Messages instance
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data])
    if type_check != {datetime.datetime}:
        m.warn('Non-date data in field.', 2)
        if datetime.time in type_check:
            m.hint('Some values _only_  contain time components', 2)


def check_field_time(data, m):

    """
    Checks for data consistency in time fields and reports to the
    Messages instance.

    Args:
        data: A list of data values, allegedly of type datetime.time
        m: A Messages instance
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data])
    if type_check != {datetime.time}:
        m.warn('Non-time formatted data found.', 2)


def check_field_taxa(data, taxa, m):

    """
    Checks if all the values provided in a Taxon field are found
    in the Taxa worksheet, reporting to the Messages instance.

    Args:
        data: A list of data values, allegedly taxon names
        taxa: A set containing taxon names from the Taxa worksheet
        m: A Messages instance
    """

    found = set(data)
    if taxa is None:
        m.warn('Taxa worksheet not provided or no taxon names were found', 2)
    if not found.issubset(taxa):
        m.warn('Includes taxa missing from Taxa worksheet:'
        ' {}'.format(', '.join(found - taxa)), 2)


def check_field_locations(data, locations, m):

    """
    Checks if all the values provided in a Locations field are
    found in the Locations worksheet, reporting to the Messages instance.

    Args:
        data: A list of data values, allegedly taxon names
        locations: A set containing locations from the Locations worksheet
        m: A Messages instance
    """

    # check if taxa are all provided
    found = set(data)
    if locations is None:
        m.warn('No Locations worksheet provided', 2)
    elif not found.issubset(locations):
        m.warn('Includes locations missing from Locations worksheet:'
        ' {}'.format(', '.join(found - locations)), 2)


def check_field_abundance(meta, data, taxa, taxa_fields,  m):

    """
    Checks abundance type data, reporting to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly numeric
        taxa: A set containing taxon names from the Taxa worksheet
        taxa_fields: A list of Taxa fields in this worksheet.
        m: A Messages instance
    """

    # check the required descriptors
    tx_ok = _check_meta(meta, 'taxon_name', m)
    mt_ok = _check_meta(meta, 'method', m)

    if tx_ok and meta['taxon_name'] not in taxa and meta['taxon_name'] not in taxa_fields:
        m.warn('Taxon name neither in the Taxa worksheet nor the name of a Taxa field', 2)

    # Can still check values are numeric, whatever happens above.
    # We're not going to insist on integers here - could be mean counts.
    is_numeric = [isinstance(vl, numbers.Number) for vl in data]
    if not all(is_numeric):
        m.warn('Field contains non-numeric data', 2)


def check_field_categorical(meta, data, m):

    """
    Checks factor data, reporting to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly strings from a provided set of levels
        m: A Messages instance
    """

    # this has already been tested but make it robust
    ct_ok = _check_meta(meta, 'levels', m)

    if not ct_ok:
        # Can't check further if no levels descriptor
        pass
    elif ct_ok and type(meta['levels']) is not unicode:
        # Can't really check anything here either
        m.warn('Category description does not seem to be text', 2)
    else:
        # Now we can test if the labels match up
        # - strip terminal semicolon if present
        if meta['levels'][-1] == ';':
            meta['levels'] = meta['levels'][:-1]
        # - split the text up by semi-colon
        levels = meta['levels'].split(';')
        # - strip off any description after a colon
        level_labels = set([ct.split(':')[0] for ct in levels])

        # - check for integer level names
        integer_codes = [is_integer_string(vl) for vl in level_labels]
        if any(integer_codes):
            m.warn('Integer level names not permitted', 2)

        # Now look for consistency: get the unique values reported in the
        # data, convert to unicode to handle checking of integer labels and
        # then check the reported levels are a subset of the descriptors.
        reported = set(data)
        reported = {unicode(lv) for lv in reported}

        if not reported.issubset(level_labels):
            m.warn('Categories found in data missing from description: '
                   '{}'.format(', '.join(reported - level_labels)), 2)


def check_field_numeric(meta, data, m):

    """
    Checks numeric type data, reporting to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly numeric
        m: A Messages instance
    """
    un_ok = _check_meta(meta, 'units', m)
    mt_ok = _check_meta(meta, 'method', m)

    # Regardless of the outcome of the meta checks, can still check the
    # data is all numeric, as it claims to be.
    is_numeric = [isinstance(vl, numbers.Number) for vl in data]
    if not all(is_numeric):
        m.warn('Non numeric data found', 2)


def check_field_trait(meta, data, taxa, taxa_fields, m):

    """
    Checks trait type data - things measured on an organism - and
    reports to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly numeric
        taxa: A set containing taxon names from the Taxa worksheet
        taxa_fields: A list of Taxa fields in this worksheet.
        m: A Messages instance
    """

    un_ok = _check_meta(meta, 'units', m)
    mt_ok = _check_meta(meta, 'method', m)
    tx_ok = _check_meta(meta, 'taxon_name', m)

    # check we can find the taxon that the trait refers to
    if tx_ok and meta['taxon_name'] not in taxa and meta['taxon_name'] not in taxa_fields:
        m.warn('Taxon name neither in the Taxa worksheet nor the name of a Taxa field', 2)

    # Regardless of the outcome of the meta checks, can still check the
    # data is all numeric, as it claims to be.
    is_numeric = [isinstance(vl, numbers.Number) for vl in data]
    if not all(is_numeric):
        m.warn('Non numeric data found', 2)


# High level functions

def check_file(fname, verbose=True):

    """
    Runs the format checking across an Excel workbook.

    Parameters:
        fname: Path to an Excel file
        verbose: Boolean to indicate whether to print messages as the program runs?

    Returns:
        A dictionary containing descriptive information about the workbook:
        - messages: The Messages instance
        - summary: The contents of the summary Worksheet metadata
    """

    try:
        wb = openpyxl.load_workbook(filename=fname, data_only=True, read_only=True)
    except:
        raise IOError('Could not open file {}'.format(fname))

    # now that we have a file, initialise the message tracker
    m = Messages("Checking file '{}'".format(fname), verbose)

    # check the metadata sheets
    summary = get_summary(wb, m)
    locations = get_locations(wb, m)
    taxa = get_taxa(wb, m)

    if 'datasets' in summary and len(summary['datasets']):
        for ds in summary['datasets']:
            check_data_worksheet(wb, ds['worksheet'], taxa, locations, m)
    else:
        m.info('No data worksheets found')

    if verbose:
        if m.count_warnings():
            m.info('FAIL: file contained {} errors'.format(m.count_warnings()))
        else:
            m.info('PASS: file formatted correctly')

    return {'messages': m, 'summary':summary}


def main():

    check_file(sys.argv[1])


if __name__ == "__main__":
    main()

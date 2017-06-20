#!/usr/bin/python

from __future__ import print_function
import openpyxl
import sys
import datetime
import re
from collections import Counter
from StringIO import StringIO

"""
Functions to check the contents of a SAFE project formatted Excel dataset.
Where possible, the functions will try to keep on running to give a list

"""

# Global definition of set of acceptable NA values: 'NA'
NA = {u'NA'}

# define some regular expressions used to check validity
re_orcid = re.compile('[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{4}')
re_email = re.compile('\S+@\S+\.\S+')
re_name = re.compile('[^,]+,[ ]?[^,]+')
re_whitespace_only = re.compile('^\s*$')


class Messages(object):

    """
    This class keeps a record of messages to be conveyed to the
    user and keeps a tally of warnings.
    """

    def __init__(self, header, verbose):
        self.header = header
        self.messages = []
        self.verbose = verbose
        self.kinds = {'info': '-', 'hint': '?', 'warn': '!'}
        self.count_new_warnings = 0
        if verbose:
            print(self.header)

    # several different levels of message:
    # - info: just reports progress and when all is ok.
    # - hint: things we might like differently, but we'll let pass
    # - warn: this is not ok.

    def info(self, message, level=0):
        msg = ('info', message, level)
        self.messages.append(msg)
        if self.verbose:
            self.print_msg(*msg)

    def warn(self, message, level=0):
        msg = ('warn', message, level)
        self.messages.append(msg)
        self.count_new_warnings += 1
        if self.verbose:
            self.print_msg(*msg)

    def hint(self, message, level=0):
        msg = ('hint', message, level)
        self.messages.append(msg)
        if self.verbose:
            self.print_msg(*msg)

    # Functions to print a message to the screen and to
    # return the text report as a StringIO object

    def print_msg(self, kind, msg, level):
        print('  ' * level + self.kinds[kind] + ' ' + msg)

    def report(self):
        report = StringIO()
        report.writelines([self.header])
        msgs = ['  ' * lv + self.kinds[kn] + msg for lv, kn, msg in self.messages]
        report.writelines(msgs)
        return report

    # Functions to check warnings
    def new_warnings(self):
        if self.count_new_warnings:
            new = self.count_new_warnings
            self.count_new_warnings = 0
            return new
        else:
            return 0

    def count_warnings(self):
        return sum([msg[0] == 'warn' for msg in self.messages])

    def no_warnings(self):
        if self.count_warnings() > 0:
            return False
        else:
            return True


def get_summary(wb, m):

    """
    Checks the information in the summary worksheet and looks for the metadata and
    dataset worksheets. The function is intended to try and get as much information
    as possible from the worksheet: the dictionary of metadata returned will have
    None for any missing data, which should be handled by downstream code.

    Parameters:
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
                "Worksheet name", "Worksheet title", "Worksheet description", "Worksheet type", }

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
    ds_keys = ['Worksheet name', 'Worksheet title', 'Worksheet description', 'Worksheet format']
    valid_sheets = set(wb.get_sheet_names())

    if set(summary_dict.keys()).issuperset(ds_keys):
        datasets = zip(*[summary_dict[x] for x in ds_keys])

        # remove any completely blank entries
        datasets = [x for x in datasets if x != tuple([None] * 4)]

        # validate
        for ind, dt in enumerate(datasets):
            dt = {k: v for k, v in zip(['worksheet', 'title', 'description', 'format'], dt)}
            # look for missing values
            for k, v in dt.iteritems():
                if v is None:
                    m.warn('Dataset {} missing'.format(k), 1)
            # validate
            if dt['worksheet'] not in valid_sheets:
                m.warn('Dataset worksheet not found: {}.'.format(dt['worksheet']), 1)
            if dt['format'] not in ['Dataframe', 'Matrix']:
                m.warn('Dataset format not valid: {}.'.format(dt['format']), 1)
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

    new_warn = m.new_warnings()
    if new_warn:
        m.info('Summary contains {} errors'.format(new_warn), 1)
    else:
        m.info('Summary formatted correctly', 1)

    return ret_dict


def get_locations(wb, m):

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

    new_warn = m.new_warnings()
    if new_warn:
        m.info('Locations contains {} errors'.format(new_warn), 1)
    else:
        m.info('{} locations loaded correctly'.format(len(loc_names)), 1)

    return loc_names


def get_taxa(wb, m):

    """
    Checks and reads the content of the taxa worksheet.

    Parameters:
        wb: An openpyxl Workbook instance.
        m: A Messages instance.
    Returns:
        A set of the taxon names to be used in the datasets.
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

    new_warn = m.new_warnings()
    if new_warn:
        m.info('Taxa contains {} errors'.format(new_warn), 1)
    else:
        m.info('{} taxa loaded correctly'.format(len(taxon_names)), 1)

    return set(taxon_names)


def check_dataframe(wb, ws_name, taxa, locations, m):

    """
    This function carries out basic checks on a 'dataframe' worksheet -
    that is, tabular data with fields of different types of variable

    Parameters:
        wb: An openpyxl Workbook instance
        ws: The worksheet name
        taxa: A list of valid taxa
        locations: A list of valid locations
        m: A Messages instance

    :return:
    """

    m.info('Checking dataframe {}'.format(ws_name))
    # get the worksheet and data dimensions
    ws = wb[ws_name]
    max_col = ws.max_column
    max_row = ws.max_row

    # get the metadata field names
    # - first search the first 20 rows for the 'field_name' descriptor which
    #   shows the end of the metadata and the start of the data
    descriptors = [ws.cell(column=1, row=rw).value for rw in range(1, 21)]
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

    # check field types - allow for empty field metadata, could be:
    # i) spacer column within fields
    # ii) the user has some comments or noise out to the right (or did,
    #     and Excel just can't let it go).
    # We don't much care about this.
    ft_found = set(fld['field_type'] for fld in field_metadata)
    ft_known = {'Location', 'Date', 'Time', 'Datetime', 'Categorical',
                'Continuous', 'Abundance', 'Comments', 'Taxa', None}

    # check extra descriptors
    if 'Categorical' in ft_found and 'categories' not in descriptors:
        m.warn('Categorical data fields found but no category descriptions provided.', 1)

    if 'Continuous' in ft_found and 'units' not in descriptors:
        m.warn('Continuous data fields found but no units provided.', 1)

    if 'Location' in ft_found and locations is None:
        m.warn('Location field found but no Location worksheet provided.', 1)

    # check field names unique (drop None)
    fn_found = Counter([fld['field_name'] for fld in field_metadata
                        if fld['field_name'] is not None])
    duplicated = [k for k, v in fn_found.iteritems() if v > 1]
    if len(duplicated) > 0:
        m.warn('Field names duplicated: {}'.format(', '.join(duplicated)), 1)

    # check mandatory fields

    # check each field
    for idx, meta in enumerate(field_metadata):

        # skip any field with no metadata
        if set(meta.values()) == {None}:
            break

        # prep the messages instance to pass to functions
        m.info('Checking field {field_name}'.format(**meta), 1)

        # read the values
        data_block = ws.get_squared_range(idx + 2, field_name_row + 1, idx + 2, max_row)
        data = [cl[0].value for cl in data_block]

        # check for missing data
        check_missing_data(data, m)

        # run consistency checks where needed
        if meta['field_type'] == 'Categorical':
            check_field_categorical(meta, data, m)
        elif meta['field_type'] in ['Date', 'Datetime']:
            check_field_datelike(meta, data, m)
        elif meta['field_type'] == 'Time':
            check_field_time(meta, data, m)
        elif meta['field_type'] == 'Taxa':
            check_field_taxa(meta, data, taxa, m)
        elif meta['field_type'] == 'Location':
            check_field_locations(meta, data, locations, m)
        elif meta['field_type'] is None:
            m.warn('Field type is empty', 2)
        else:
            m.warn('Unknown field type {field_type}'.format(**meta), 2)

    if m.no_warnings():
        m.info('Dataframe formatted correctly', 1)
    else:
        m.info('Dataframe contains {} errors'.format(m.count_warnings()), 1)

    return None




def check_missing_data(data, m):

    # Only NA is acceptable
    na_vals = [vl in NA for vl in data]
    if any(na_vals):
        m.hint('{} / {} values missing'.format(sum(na_vals), len(na_vals)), 2)

    # We won't tolerate:
    # 1) empty cells (just to avoid ambiguity - e.g. in abundance data)
    is_empty = [vl is None for vl in data]
    if sum(is_empty):
        m.warn('cells are empty'.format(sum(is_empty)), 2)
    # 2) non-empty cells containing only whitespace strings
    ws_only = [re_whitespace_only.match(vl) for vl in data if type(vl) in [str, unicode]]
    if any(ws_only):
        m.warn('{} cells contain whitespace only text'.format(sum(ws_only)), 2)


"""
Field checker functions below. Must be capable of handling the NA values elegantly
"""


def check_field_datelike(meta, data, m):

    """
    Checks for consistency in date and datetime fields

    :param meta:
    :param data:
    :return:
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data if vl not in NA])
    if type_check != {datetime.datetime}:
        if datetime.time in type_check:
            m.warn('Some data only contains time information, not date.', 2)
        else:
            m.warn('Non-date formatted data found.', 2)

    # Check no time component in date fields
    if meta['field_type'] == 'Date':
        no_time = [vl.time() == datetime.time(0, 0) for vl in data if vl not in NA]
        if not all(no_time):
            m.warn('Datetimes found in field of type Date.', 2)


def check_field_time(meta, data, m):

    """
    Checks for consistency in time fields

    :param meta:
    :param data:
    :return:
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data if vl not in NA])
    if type_check != {datetime.time}:
        m.warn('Non-time formatted data found.', 2)


def check_field_categorical(meta, data, m):

    # this has already been tested but make it robust
    if 'categories' not in meta:
        m.warn('Category levels metadata not provided', 2)
    elif meta['categories'] is None:
        m.warn('Category descriptions empty', 2)
    elif type(meta['categories']) is not unicode:
        m.warn('Category description does not seem to be text', 2)
    else:
        # split the text up by semi-colon
        categories = meta['categories'].split(';')
        # strip off any description after a colon
        category_labels = set([ct.split(':')[0] for ct in categories])

        # get the unique values reported in the data, removing missing values and
        # then convert to unicode in case anyone has used numeric level names.
        reported = set(data) - NA
        reported = {unicode(lv) for lv in reported}

        # check for completeness
        if not reported.issubset(category_labels):
            m.warn('Categories found in data missing from description'
                   '{}: '.format(','.join(reported - category_labels)), 2)


def check_field_continuous(meta, data, m):

    """
    Checks for consistency in time fields

    :param meta:
    :param data:
    :return:
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data if vl not in NA])
    if not type_check.issubset({long, float}):
        m.warn('Non numeric data found in {}.'.format(meta['field_name']), 2)


def check_field_taxa(meta, data, taxa, m):

    # check if taxa are all provided
    found = set(data)
    if not found.issubset(taxa):
        m.warn('Includes taxa missing from Taxa worksheet:'
               ' {}'.format(', '.join(found - taxa)), 2)

def check_field_locations(meta, data, locations, m):

    # check if taxa are all provided
    found = set(data)
    if not found.issubset(locations):
        m.warn('Includes locations missing from Locations worksheet:'
               ' {}'.format(', '.join(found - locations)), 2)


def check_file(fname, verbose=True):

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
            if ds['format'] == 'Dataframe':
                check_dataframe(wb, ds['worksheet'], taxa, locations, m)
            elif ds['format'] == 'Matrix':
                pass
                # check_matrix(wb, ds['name'], taxa, locations)
            else:
                m.warn("Data worksheet {worksheet} in "
                       "unknown format {format}".format(**ds))
    else:
        m.info('No data worksheets found')

    if verbose:
        if m.count_warnings():
            m.info('FAIL: file contained {} errors'.format(m.count_warnings()))
        else:
            m.info('PASS: file formatted correctly')

    return m


def main():

    check_file(sys.argv[1])


if __name__ == "__main__":
    main()

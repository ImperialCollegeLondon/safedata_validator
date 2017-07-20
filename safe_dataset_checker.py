#!/usr/bin/python

"""
Module containing code to verify the format of a SAFE project Excel dataset.

The functions are written to extract as much information as possible from
checking a file so, rather than raising an error and halting at the first
fault, functions are written to check as much as they can on the inputs.
Any information and warnings are added on to a Messages instance, passed
to each function, which is used to keep a report of issues throughout the
file check.

Taxon names and locations are both validated:
i) Taxonomic names are validated against the NCBI database using either Entrez
   queries or the ete2 python package. The ete2 package installation creates
   a 300 MB local SQLITE version of the NCBI data, but does allow offline use.
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

# if there is a local install of the ete2 package, with a built database
# then use that for speed. Otherwise the code will try and use Entrez
try:
    from ete2 import NCBITaxa
    ncbi_db = os.path.join(os.environ.get('HOME', '/'), '.etetoolkit', 'taxa.sqlite')
    if not os.path.exists(ncbi_db):
        raise RuntimeError()
    USE_ETE = True
except (ImportError, RuntimeError) as e:
    USE_ETE = False


# define some regular expressions used to check validity
RE_ORCID = re.compile(r'[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{4}')
RE_EMAIL = re.compile(r'\S+@\S+\.\S+')
RE_NAME = re.compile(r'[^,]+,[ ]?[^,]+')
RE_WSPACE_ONLY = re.compile(r'^\s*$')


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
        n_warnings: The number of warning messages added.
    """

    def __init__(self, header, verbose=True):
        self.header = header
        self.messages = []
        self.verbose = verbose
        self.kinds = {'info': '-', 'hint': '?', 'warn': '!'}
        self.n_warnings = 0
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
        self.n_warnings += 1
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
        report.writelines([self.header + '\n'])
        msgs = ['  ' * lv + self.kinds[kn] + msg + '\n' for kn, msg, lv in self.messages]
        report.writelines(msgs)
        return report


def get_summary(workbook, msg):

    """
    Checks the information in the summary worksheet and looks for the metadata and
    dataset worksheets. The function is intended to try and get as much information
    as possible from the worksheet: the dictionary of metadata returned will have
    None for any missing data, which should be handled by downstream code.

    Args:
        workbook: An openpyxl Workbook instance
        msg: A Messages instance
    Returns:
        A dictionary of the available summary metadata. The contents should be
        serialisable using simplejson, so that the summary can be stored in a
        database field.
    """

    # try and get the summary worksheet
    msg.info("Checking Summary worksheet")
    start_warn = msg.n_warnings
    try:
        worksheet = workbook['Summary']
    except KeyError:
        msg.warn("Summary worksheet not found, moving on.")
        return None

    # load dictionary of summary information block, allowing for multiple
    # columns for fields (data compilers, dataset sheets).
    summary_dict = {}
    col_range = range(2, worksheet.max_column + 1)
    max_row = worksheet.max_row

    for row in range(1, max_row + 1):
        vals = [worksheet.cell(row=row, column=cl).value for cl in col_range]
        summary_dict[worksheet.cell(row=row, column=1).value] = vals

    # build return dictionary
    ret_dict = {}

    # Check the minimal keys are expected
    required = {"SAFE Project ID", "Access status", "Title", "Description",
                "Author name", "Author email", "Author affiliation", "Author ORCID",
                "Worksheet name", "Worksheet title", "Worksheet description", 'Keywords'}

    found = set(summary_dict.keys())

    # don't bail here - try and get as far as possible
    if not found.issuperset(required):
        msg.warn('Missing metadata fields: {}'.format(', '.join(required - found)), 1)

    # check for contents
    if any([set(x) == {None} for x in summary_dict.values()]):
        msg.warn('Metadata fields with no information.', 1)

    # CHECK PROJECT ID
    if 'SAFE Project ID' not in summary_dict:
        msg.warn('SAFE Project ID missing', 1)
        ret_dict['project_id'] = None
    else:
        if not isinstance(summary_dict['SAFE Project ID'][0], long):
            msg.warn('SAFE Project ID is not an integer.', 1)
        ret_dict['project_id'] = summary_dict['SAFE Project ID'][0]

    # CHECK DATASET TITLE
    if 'Title' not in summary_dict:
        msg.warn('Dataset title missing', 1)
        ret_dict['title'] = None
    else:
        ret_dict['title'] = summary_dict['Title'][0]

    # CHECK DATASET DESCRIPTION
    if 'Description' not in summary_dict:
        msg.warn('Dataset title missing', 1)
        ret_dict['description'] = None
    else:
        ret_dict['description'] = summary_dict['Description'][0]

    # CHECK ACCESS STATUS AND EMBARGO DETAILS
    if 'Access status' in summary_dict:
        access_status = summary_dict['Access status'][0]
        ret_dict['access_status'] = summary_dict['Access status'][0]
        ret_dict['embargo_date'] = None
        if access_status == 'Embargo':
            if 'Embargo date' not in summary_dict:
                msg.warn('Dataset embargoed but no date provided.', 1)
            embargo_date = summary_dict['Embargo date'][0]
            now = datetime.datetime.now()
            if not isinstance(embargo_date, datetime.datetime):
                msg.warn('Embargo date not formatted as date.', 1)
            elif embargo_date < now:
                msg.warn('Embargo date is in the past.', 1)
            elif embargo_date > now + datetime.timedelta(days=2*365):
                msg.warn('Embargo date more than two years in the future.', 1)
            else:
                ret_dict['embargo_date'] = embargo_date.date().isoformat()
        elif access_status == 'Open':
            pass
        else:
            msg.warn('Access status must be Open or Embargo not {}'.format(access_status), 1)

    # CHECK KEYWORDS
    if 'Keywords' not in summary_dict:
        msg.warn('Dataset keywords missing', 1)
        ret_dict['keywords'] = None
    else:
        # drop any blanks
        keywords = [vl for vl in summary_dict['Keywords']
                    if vl is not None or not RE_WSPACE_ONLY.match(vl)]
        ret_dict['keywords'] = keywords

    # CHECK AUTHORS
    author_keys = ['Author name', 'Author affiliation', 'Author email', 'Author ORCID']
    if set(summary_dict.keys()).issuperset(author_keys):
        authors = zip(*[summary_dict[x] for x in author_keys])

        # remove any completely blank entries
        authors = [x for x in authors if x != tuple([None] * 4)]

        # convert to dict in Zenodo style and check completeness and validity
        # - ORCID isn't mandatory
        for ind, auth in enumerate(authors):
            auth = {k: v for k, v in zip(['name', 'affiliation', 'email', 'orcid'], auth)}
            # look for missing values
            for key, val in auth.iteritems():
                if val is None and key != 'orcid':
                    msg.warn('Author {} missing'.format(key), 1)
            # check validity
            if not RE_NAME.match(auth['name']):
                msg.warn('Author name not formated as last_name, '
                         'first_names: {}'.format(auth['name']), 1)

            if auth['orcid'] is None:
                msg.hint('Consider adding an ORCiD!', 1)
            elif not RE_ORCID.match(str(auth['orcid'])):
                msg.warn('ORCID not properly formatted: {}'.format(auth['orcid']), 1)

            if not RE_EMAIL.match(auth['email']):
                msg.warn('Email not properly formatted: {}'.format(auth['email']), 1)
            authors[ind] = auth

        ret_dict['authors'] = authors
    else:
        msg.warn('Author metadata block incomplete.', 1)
        ret_dict['authors'] = None

    # CHECK DATA WORKSHEETS
    ws_keys = ['Worksheet name', 'Worksheet title', 'Worksheet description']
    valid_sheets = set(workbook.get_sheet_names())

    if set(summary_dict.keys()).issuperset(ws_keys):
        data_worksheets = zip(*[summary_dict[x] for x in ws_keys])

        # remove any completely blank entries
        data_worksheets = [x for x in data_worksheets if x != tuple([None] * 3)]

        # validate
        for ind, data_ws in enumerate(data_worksheets):
            data_ws = {k: v for k, v in zip(['worksheet', 'title', 'description'], data_ws)}
            # look for missing values
            for key, val in data_ws.iteritems():
                if val is None:
                    msg.warn('Data worksheet {} missing'.format(key), 1)
            # validate
            if data_ws['worksheet'] not in valid_sheets:
                msg.warn('Dataset worksheet not found: {}.'.format(data_ws['worksheet']), 1)
            data_worksheets[ind] = data_ws

        ret_dict['data_worksheets'] = data_worksheets

        # check for extra undocumented spreadsheets
        if 'Worksheet name' in summary_dict:
            expected_sheets = set([vl['worksheet'] for vl in data_worksheets]
                                  + ['Summary', 'Taxa', 'Locations'])
            if valid_sheets != expected_sheets:
                msg.warn('Undocumented sheets found in workbook: '
                         '{}'.format(', '.join(valid_sheets - expected_sheets)), 1)
    else:
        msg.warn('Data worksheet metadata block incomplete.', 1)
        ret_dict['data_worksheets'] = None

    # summary of processing
    if (msg.n_warnings - start_warn) > 0:
        msg.info('Summary contains {} errors'.format(msg.n_warnings - start_warn), 1)
    else:
        msg.info('Summary formatted correctly', 1)

    return ret_dict


def get_locations(workbook, msg, locations_json=None):

    """
    Attempts to load and check the contents of the Locations worksheet.
    Args:
        workbook: An openpyxl Workbook instance
        msg: A Messages instance
        locations_json: A path to a JSON file containing a valid set of location names.
            With the default value of None, the function tries to get this from
            a SAFE project website service.
    Returns:
        A set of provided location names or None if the worksheet cannot
        be found or no locations can be loaded.
    """

    # try and get the locations worksheet
    msg.info("Checking Locations worksheet")
    start_warn = msg.n_warnings
    try:
        locs = workbook['Locations']
    except KeyError:
        # No locations is pretty implausible, but still persevere as if
        # they aren't going to be required
        msg.warn("No locations worksheet found - moving on", 1)
        return None

    # get the set of valid names
    if locations_json is None:
        # If no file is provided then try and get locations from the website service
        r = requests.get('https://www.safeproject.net/call/json/get_locations')
        if r.status_code != 200:
            msg.warn('Could not download valid location names. Use a local json file.', 1)
            valid_locations = []
        else:
            valid_locations = r.json()['locations']
    else:
        # try and load the file
        try:
            r = simplejson.load(file(locations_json))
            valid_locations = r['locations']
        except IOError:
            msg.warn('Could not load location names from file.', 1)
            valid_locations = []

    # Check the headers
    fields = []
    max_col = locs.max_column
    max_row = locs.max_row

    for col in range(1, max_col):
        fields.append(locs.cell(row=1, column=col).value)

    # check the key fields are there
    if not set(fields).issuperset({'Location name'}):
        msg.warn('Location name column not found', 1)
        loc_names = None
    else:
        # get the location names as strings
        names_col = fields.index('Location name') + 1
        rows = range(2, max_row + 1)
        loc_names = [str(locs.cell(row=rw, column=names_col).value) for rw in rows]

        # check for duplicate names
        if len(set(loc_names)) != len(loc_names):
            msg.warn('Duplicated location names', 1)

        loc_names = set(loc_names)

        # Are they valid?
        invalid = loc_names - set(valid_locations)
        if len(invalid) > 0:
            msg.warn('Invalid locations found: ' + ', '.join(invalid), 1)

    if (msg.n_warnings - start_warn) > 0:
        msg.info('Locations contains {} errors'.format(msg.n_warnings - start_warn), 1)
    else:
        msg.info('{} locations loaded correctly'.format(len(loc_names)), 1)

    return loc_names


def get_taxa(workbook, msg, force_entrez=False):

    """
    Attempts to load and check the content of the Taxa worksheet. The
    function checks that:
    i)   all taxa have a taxon name and a taxon type,
    ii)  that taxa have a complete taxonomic hierarchy up to the taxon type
    iii) that all taxonomic names are known to NCBI, unless they are
         explicitly marked as new species with an asterisk suffix.

    Args:
        workbook: An openpyxl Workbook instance.
        msg: A Messages instance.
        force_entrez: Use entrez queries even when a local NCBI query is available
    Returns:
        A set of provided taxon names or None if the worksheet cannot
        be found or no taxa can be loaded.
    """

    if force_entrez:
        use_ete = False
    else:
        use_ete = USE_ETE

    # try and get the taxon worksheet
    msg.info("Checking Taxa worksheet")
    start_warn = msg.n_warnings
    try:
        taxa = workbook['Taxa']
    except KeyError:
        # This might mean that the study doesn't have any taxa, so return an empty
        # set. If the datasets then contain taxonomic names, it'll fail gracefully.
        msg.hint("No taxa worksheet found - assuming no taxa in data for now!", 1)
        return set()

    # Check the headers
    max_col = taxa.max_column
    max_row = taxa.max_row

    # load the data
    cells = taxa.get_squared_range(1, 1, max_col, max_row)
    data = [[cl.value for cl in row] for row in cells]

    # remove null columns and rows
    null_row = [set(rw) == {None} for rw in data]
    data = [rw for rw, nl in zip(data, null_row) if not nl]
    data = zip(*data)
    null_col = [set(cl) == {None} for cl in data]
    data = [rw for rw, nl in zip(data, null_col) if not nl]
    data = zip(*data)

    # split off the field headers from the data
    fields = data[0]
    data = data[1:]

    # Look for the names and types and check pairwise complete
    if not set(fields).issuperset({'Taxon name', 'Taxon type'}):
        msg.warn('One or both of Taxon name and Taxon type columns not found', 1)
        return None

    # turn the taxa into dictionaries
    taxa = [dict(zip(fields, row)) for row in data]

    # Look for rows with either or both of taxon type and taxon name missing.
    # We've cleared rows with all None.
    incomplete_pairs = [((rw['Taxon name'] is None) or (rw['Taxon type'] is None))
                        for rw in taxa]
    if any(incomplete_pairs):
        msg.warn('Rows with taxon name and/or taxon type missing', 1)

    # check for duplicate names
    taxon_names = [tx['Taxon name'] for tx in taxa]
    if len(set(taxon_names)) != len(taxon_names):
        msg.warn('Duplicated taxon names', 1)

    # Now check the taxonomy
    # - We will only check the names for the big seven taxonomic levels
    taxon_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    # - What types are set by the user?
    taxon_types = set([tx['Taxon type'] for tx in taxa])
    # - Allow for missing taxon types, since that will already have issued a warning
    available_types = set(fields) | {None}
    if not taxon_types.issubset(available_types):
        msg.warn('Some rows have taxon types that do not match a column name: '
                 '{}'.format(', '.join(taxon_types - available_types)), 1)

    # What levels are required given the taxon types?
    max_taxon = max([taxon_levels.index(fld) for fld in fields if fld in taxon_levels])
    required_levels = set(taxon_levels[:max_taxon + 1])
    if not set(fields).issuperset(required_levels):
        m = 'Required columns missing: {}'
        m = m.format(', '.join(required_levels - set(fields)))
        msg.warn(m, 1)

    # Now check the available taxon columns against the NCBI database
    #  - NCBI knows about the following ranks:
    #    ncbi.db.execute('select distinct rank from species;').fetchall()
    #  - The user can provide these but we currently only check and index the major seven

    if use_ete:
        ncbi = NCBITaxa()
        msg.info('Using local NCBI database to validate names', 1)
    else:
        msg.info('Using Entrez queries to validate names', 1)

    # get the levels that are both required and available, and retain order
    available_types = [tx for tx in taxon_levels if tx in (set(fields) & required_levels)]

    # collect information for a taxonomic index
    tax_index = []

    for rnk in available_types:
        msg.info('Checking validity of {} level taxa'.format(rnk), 1)

        # a) Read the unique entries in that column, skipping None
        tax_vals = set([tx[rnk] for tx in taxa if tx[rnk] is not None])

        # b) Check for whitespace only and drop them
        ws_only = [RE_WSPACE_ONLY.match(x) for x in tax_vals]
        if any(ws_only):
            msg.warn('Cells with only whitespace text', 2)
            tax_vals = [tx for tx, ws in zip(tax_vals, ws_only) if not ws]

        # c) look for new species flagged as new, so they don't get tested
        new = set([tx for tx in tax_vals if tx.endswith('*')])
        if new:
            msg.hint('New taxa reported: {}'.format(', '.join(new)), 2)

        # c) Validate the remaining names at the given rank in NCBI using either:
        #    i) a local database via ete2 or
        #    ii) entrez queries over the web, which will be slow.
        to_validate = set(tax_vals) - new

        if use_ete:
            # look up the names
            ids = ncbi.get_name_translator(to_validate)
            # isolate names that aren't found
            not_found = set(tax_vals) - set(ids.keys())
            # check the ranks
            ranks_found = {ky: ncbi.get_rank(vl).values() for ky, vl in ids.iteritems()}
            bad_rank = set([ky for ky, vl in ranks_found.iteritems() if rnk.lower() not in vl])
            invalid = not_found & bad_rank
            # all taxa are validated when using local query
            unvalidated = set()
        else:
            # loop the names through the entrez interface to the taxonomy database
            # - search query to see if the name exists at the given rank
            entrez = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
                      'db=taxonomy&term={}+AND+{}[rank]&retmode=json')
            # compile problems into sets
            unvalidated = set()
            invalid = set()
            # loop over the values
            for tx in tax_vals:
                tx_id = requests.get(entrez.format(tx, rnk))
                if tx_id.status_code != 200:
                    # internet failures just add individual taxa to the unvalidated list
                    unvalidated.add(tx)
                else:
                    if tx_id.json()['esearchresult']['count'] == 0:
                        invalid.add(tx)

        # d) now report invalid taxa and unvalidated
        if invalid:
            msg.warn('Taxa not valid at this rank: {}'.format(', '.join(invalid)), 2)

        if unvalidated:
            msg.warn('Entrez taxon validation failed for: {}'.format(', '.join(invalid)), 2)

        tax_index.extend([(tx, rnk.lower()) for tx in tax_vals])

    # Finally, check completeness of taxonomic hierarchy
    msg.info('Checking taxon hierarchies complete'.format(rnk), 1)
    for tx in taxa:
        tx_tp = tx['Taxon type']

        # set level to check to for non-taxonomic levels
        if tx_tp == 'Morphospecies':
            tx_tp = 'Order'
            wrn = 'Morphospecies taxonomy not complete to {} level: {}'
        elif tx_tp == 'Functional Group':
            tx_tp = 'Kingdom'
            wrn = 'Functional group taxonomy not complete to {} level: {}'
        else:
            wrn = 'Taxonomy not complete to {} level: {}'

        # If the type is available, are all levels up to and including
        # it complete and non blank.
        if tx_tp is not None and tx_tp in tx:
            to_check = available_types[0:available_types.index(tx_tp) + 1]
            filled = [False if ((tx[rnk] is None) or (RE_WSPACE_ONLY.match(tx[rnk])))
                      else True for rnk in to_check]
            if not all(filled):
                msg.warn(wrn.format(tx_tp, tx['Taxon name']), 2)
        elif tx_tp is not None and tx_tp not in tx:
            msg.warn('Stated taxon type column {} is missing: {}'.format(tx_tp, tx['Taxon name']), 2)
        else:
            msg.warn('Taxon type blank: {}'.format(tx['Taxon name']), 2)

    if (msg.n_warnings - start_warn) > 0:
        msg.info('Taxa contains {} errors'.format(msg.n_warnings - start_warn), 1)
    else:
        msg.info('{} taxa loaded correctly'.format(len(taxon_names)), 1)

    return set(taxon_names), tax_index


def check_data_worksheet(workbook, ws_meta, taxa, locations, msg):

    """
    Attempt to load and checks the formatting and content of a data worksheet,
    updating the Messages instance `m` with the results.

    The function returns the meta dictionary, updated with available information
    on the size of the worksheet and the fields in it.

    Args:
        workbook: An openpyxl Workbook instance
        ws_meta: The metadata for this sheet from the summary
        taxa: A list of valid taxa
        locations: A list of valid locations
        msg: A Messages instance

    Returns:
        A possibly updated version of ws_meta.
    """

    ws_name = ws_meta['worksheet']
    if ws_name not in workbook.get_sheet_names():
        msg.warn('Data worksheet {} not found'.format(ws_name))
        return ws_meta
    else:
        msg.info('Checking data worksheet {}'.format(ws_name))
        start_warn = msg.n_warnings

    # get the worksheet and data dimensions
    worksheet = workbook[ws_name]
    max_col = worksheet.max_column
    max_row = worksheet.max_row

    # trap completely empty worksheets
    if max_row == 1:
        msg.warn('Worksheet is empty', 1)
        return ws_meta

    # get the metadata field names
    # - first search at most the first 20 rows for the 'field_name' descriptor
    #   which shows the end of the metadata and the start of the data
    descriptors = [worksheet.cell(column=1, row=rw).value
                   for rw in range(1, min(20, max_row) + 1)]
    if 'field_name' in descriptors:
        field_name_row = descriptors.index('field_name') + 1
        descriptors = descriptors[:field_name_row]
    else:
        msg.warn('Cannot parse data: field_name row not found', 1)
        return ws_meta

    # Neither of the row and col maxima are particularly reliable as Excel can hang on to
    # cell references for previously used cells. We can ignore blank columns easily but
    # we do need to know where to actually stop for finding blank data in rows.
    # So, we explicitly check the row numbers, making them a mandatory part of the setup.

    # - get the values
    row_number_cells = worksheet.get_squared_range(1, field_name_row + 1, 1, max_row)
    row_numbers = [cl[0].value for cl in row_number_cells]

    # - trim blank or whitespace values from the end and update the max col
    while row_numbers[-1] is None or RE_WSPACE_ONLY.match(str(row_numbers[-1])):
        row_numbers.pop()

    max_row = len(row_numbers) + field_name_row

    # now check the row numbers are numbers and if they are
    # do they start at one and go up by one
    if not all([isinstance(vl, numbers.Number) for vl in row_numbers]):
        msg.warn('Non-numeric data found in row numbering', 1)
    else:
        if row_numbers[0] != 1:
            msg.warn('Row numbering does not start at 1', 1)

        one_increment = [(vl1 - vl2) == 1 for vl1, vl2 in zip(row_numbers[1:], row_numbers[:-1])]
        if not all(one_increment):
            msg.warn('Row numbering does not consistently increment by 1', 1)

    # report on detected size
    msg.info('Worksheet contains {} rows and {} columns'.format(max_row, max_col), 1)
    ws_meta['ncol'] = max_col
    ws_meta['nrow'] = max_row

    # get the metadata for each field
    field_metadata = []
    for col in range(2, max_col + 1):
        this_field = {k: worksheet.cell(column=col, row=i + 1).value
                      for i, k in enumerate(descriptors)}
        field_metadata.append(this_field)

    # check required descriptors are present and if locations and taxa
    # turn out to be needed after all!
    ft_found = [fld['field_type'] for fld in field_metadata]
    if 'Categorical' in ft_found and 'levels' not in descriptors:
        msg.warn('Categorical data fields found but no levels descriptor provided.', 1)

    if 'Numeric' in ft_found:
        if 'units' not in descriptors:
            msg.warn('Numeric data fields found but no units descriptor provided.', 1)
        if 'method' not in descriptors:
            msg.warn('Numeric data fields found but no method descriptor provided.', 1)

    if 'Location' in ft_found and locations is None:
        msg.warn('Location field found but no Location worksheet provided.', 1)

    if 'Abundance' in ft_found:
        if 'taxon_name' not in descriptors:
            msg.warn('Abundance field found but no taxon name descriptor provided.', 1)
        if 'method' not in descriptors:
            msg.warn('Abundance field found but no sampling method descriptor provided.', 1)

    # check field names unique (drop None)
    fn_found = Counter([fld['field_name'] for fld in field_metadata
                        if fld['field_name'] is not None])
    duplicated = [k for k, v in fn_found.iteritems() if v > 1]
    if len(duplicated) > 0:
        msg.warn('Field names duplicated: {}'.format(', '.join(duplicated)), 1)

    # get taxa field names for cross checking observation and trait data
    taxa_fields = [fld['field_name'] for fld in field_metadata if fld['field_type'] == 'Taxa']

    # TODO - check mandatory fields

    # check each field
    for idx, meta in enumerate(field_metadata):

        # skip any field with no metadata
        if set(meta.values()) == {None}:
            break

        # prep the messages instance to pass to functions
        if meta['field_name'] is None or RE_WSPACE_ONLY.match(meta['field_name']):
            msg.info('Checking Column {}'.format(utils.get_column_letter(idx + 2)), 1)
            msg.warn('Field name is blank', 2)
        else:
            msg.info('Checking field {field_name}'.format(**meta), 1)

        # check the description
        if meta['description'] is None or RE_WSPACE_ONLY.match(meta['description']):
            msg.warn('Description is missing', 2)

        # read the values
        data_block = worksheet.get_squared_range(idx + 2, field_name_row + 1, idx + 2, max_row)
        data = [cl[0].value for cl in data_block]

        # filter out missing and blank data, except for comments fields, where
        # blanks are not an error
        if meta['field_type'] != 'Comments':
            data = filter_missing_or_blank_data(data, msg)

        # run consistency checks where needed and trap unknown field types
        if meta['field_type'] == 'Date':
            check_field_date(data, msg)
        elif meta['field_type'] == 'Datetime':
            check_field_datetime(data, msg)
        elif meta['field_type'] == 'Time':
            check_field_time(data, msg)
        elif meta['field_type'] == 'Taxa':
            check_field_taxa(data, taxa, msg)
        elif meta['field_type'] == 'Location':
            # location names should be strings
            data = [str(dt) for dt in data]
            check_field_locations(data, locations, msg)
        elif meta['field_type'] == 'Categorical':
            check_field_categorical(meta, data, msg)
        elif meta['field_type'] == 'Numeric':
            check_field_numeric(meta, data, msg)
        elif meta['field_type'] == 'Abundance':
            check_field_abundance(meta, data, taxa, taxa_fields, msg)
        elif meta['field_type'] == 'Trait':
            check_field_trait(meta, data, taxa, taxa_fields, msg)
        elif meta['field_type'] == 'Comments':
            pass
        elif meta['field_type'] is None:
            msg.warn('Field type is empty', 2)
        else:
            msg.warn('Unknown field type {field_type}'.format(**meta), 2)

    if (msg.n_warnings - start_warn) > 0:
        msg.info('Dataframe contains {} errors'.format(msg.n_warnings - start_warn), 1)
    else:
        msg.info('Dataframe formatted correctly', 1)

    # trim out field descriptions for blank columns
    field_metadata = [fld for fld in field_metadata if set(fld.values()) != {None}]

    # update ws_meta with the field information
    ws_meta['fields'] = field_metadata

    return ws_meta


def filter_missing_or_blank_data(data, msg):

    """
    Takes a list of data and filters out any missing or blank data,
    reporting to the Messages instance as it goes. The filtered list
    is returned for feeding to field type checker functions so that
    they can assume no NA or blank data in their checks.

    Args:
        data: A list of values read from the Worksheet.
        msg: A Messages instance
    Returns:
        A list of data values filtered to remove NA and blanks.
    """

    # Only NA is acceptable
    na_vals = [vl == u'NA' for vl in data]
    if any(na_vals):
        msg.hint('{} / {} values missing'.format(sum(na_vals), len(na_vals)), 2)

    # We won't tolerate:
    # 1) empty cells (just to avoid ambiguity - e.g. in abundance data)
    is_empty = [vl is None for vl in data]
    if sum(is_empty):
        msg.warn('{} cells are blank'.format(sum(is_empty)), 2)
    # 2) non-empty cells containing only whitespace strings
    ws_only = [RE_WSPACE_ONLY.match(unicode(vl)) is not None for vl in data]
    if any(ws_only):
        msg.warn('{} cells contain whitespace only text'.format(sum(ws_only)), 2)

    # Return the values that aren't NA, blank or whitespace only
    na_or_blank = [any(tst) for tst in zip(na_vals, is_empty, ws_only)]
    data = [dt for dt, nb in zip(data, na_or_blank) if not nb]

    return data


def _check_meta(meta, descriptor, msg):
    """
    A standardised check to see if a required descriptor is present for
    a field and that it isn't simply empty or whitespace. The function
    reports problems to the Messages instance and returns a boolean
    showing if the checks passed successfully.

    Args:
        meta: A dictionary of field metadata descriptors
        descriptor: The name of the descriptor to check.
        msg: An instance of class Messages

    Returns:
        A boolean, with True showing no problems and False showing
        that warnings occurred.
    """

    if descriptor not in meta:
        msg.warn('{} descriptor missing'.format(descriptor), 2)
        return False
    elif meta[descriptor] is None or RE_WSPACE_ONLY.match(meta[descriptor]):
        msg.warn('{} descriptor is blank'.format(descriptor), 2)
        return False
    else:
        return True


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


def check_field_date(data, msg):

    """
    Checks for data consistency in date fields and reports to the
    Messages instance.

    Args:
        data: A list of data values, allegedly of type datetime.date
        msg: A Messages instance
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data])
    if type_check != {datetime.datetime}:
        msg.warn('Non-date data in field.', 2)
        if datetime.time in type_check:
            msg.hint('Some values _only_  contain time components', 2)

    # Check no time component in actual dates
    no_time = [vl.time() == datetime.time(0, 0) for vl in data]
    if not all(no_time):
        msg.warn('Some values also contain time components', 2)


def check_field_datetime(data, msg):

    """
    Checks for data consistency in datetime fields and reports to the
    Messages instance.

    Args:
        data: A list of data values, allegedly of type datetime.datetime
        msg: A Messages instance
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data])
    if type_check != {datetime.datetime}:
        msg.warn('Non-date data in field.', 2)
        if datetime.time in type_check:
            msg.hint('Some values _only_  contain time components', 2)


def check_field_time(data, msg):

    """
    Checks for data consistency in time fields and reports to the
    Messages instance.

    Args:
        data: A list of data values, allegedly of type datetime.time
        msg: A Messages instance
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data])
    if type_check != {datetime.time}:
        msg.warn('Non-time formatted data found.', 2)


def check_field_taxa(data, taxa, msg):

    """
    Checks if all the values provided in a Taxon field are found
    in the Taxa worksheet, reporting to the Messages instance.

    Args:
        data: A list of data values, allegedly taxon names
        taxa: A set containing taxon names from the Taxa worksheet
        msg: A Messages instance
    """

    found = set(data)
    if taxa is None:
        msg.warn('Taxa worksheet not provided or no taxon names were found', 2)
    if not found.issubset(taxa):
        msg.warn('Includes taxa missing from Taxa worksheet:'
                 ' {}'.format(', '.join(found - taxa)), 2)


def check_field_locations(data, locations, msg):

    """
    Checks if all the values provided in a Locations field are
    found in the Locations worksheet, reporting to the Messages instance.

    Args:
        data: A list of data values, allegedly taxon names
        locations: A set containing locations from the Locations worksheet
        msg: A Messages instance
    """

    # check if locations are all provided
    found = set(data)
    if locations is None:
        msg.warn('No Locations worksheet provided', 2)
    elif not found.issubset(locations):
        msg.warn('Includes locations missing from Locations worksheet:'
                 ' {}'.format(', '.join(found - locations)), 2)


def check_field_abundance(meta, data, taxa, taxa_fields, msg):

    """
    Checks abundance type data, reporting to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly numeric
        taxa: A set containing taxon names from the Taxa worksheet
        taxa_fields: A list of Taxa fields in this worksheet.
        msg: A Messages instance
    """

    # check the required descriptors
    _check_meta(meta, 'method', msg)
    tx_ok = _check_meta(meta, 'taxon_name', msg)

    if tx_ok and meta['taxon_name'] not in taxa and meta['taxon_name'] not in taxa_fields:
        msg.warn('Taxon name neither in the Taxa worksheet nor the name of a Taxa field', 2)

    # Can still check values are numeric, whatever happens above.
    # We're not going to insist on integers here - could be mean counts.
    is_numeric = [isinstance(vl, numbers.Number) for vl in data]
    if not all(is_numeric):
        msg.warn('Field contains non-numeric data', 2)


def check_field_categorical(meta, data, msg):

    """
    Checks factor data, reporting to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly strings from a provided set of levels
        msg: A Messages instance
    """

    # this has already been tested but make it robust
    ct_ok = _check_meta(meta, 'levels', msg)

    if not ct_ok:
        # Can't check further if no levels descriptor
        pass
    elif ct_ok and not isinstance(meta['levels'], unicode):
        # Can't really check anything here either
        msg.warn('Category description does not seem to be text', 2)
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
            msg.warn('Integer level names not permitted', 2)

        # Now look for consistency: get the unique values reported in the
        # data, convert to unicode to handle checking of integer labels and
        # then check the reported levels are a subset of the descriptors.
        reported = set(data)
        reported = {unicode(lv) for lv in reported}

        if not reported.issubset(level_labels):
            msg.warn('Categories found in data missing from description: '
                     '{}'.format(', '.join(reported - level_labels)), 2)


def check_field_numeric(meta, data, msg):

    """
    Checks numeric type data, reporting to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly numeric
        msg: A Messages instance
    """
    _check_meta(meta, 'units', msg)
    _check_meta(meta, 'method', msg)

    # Regardless of the outcome of the meta checks, can still check the
    # data is all numeric, as it claims to be.
    is_numeric = [isinstance(vl, numbers.Number) for vl in data]
    if not all(is_numeric):
        msg.warn('Non numeric data found', 2)


def check_field_trait(meta, data, taxa, taxa_fields, msg):

    """
    Checks trait type data - things measured on an organism - and
    reports to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly numeric
        taxa: A set containing taxon names from the Taxa worksheet
        taxa_fields: A list of Taxa fields in this worksheet.
        msg: A Messages instance
    """

    _check_meta(meta, 'units', msg)
    _check_meta(meta, 'method', msg)
    tx_ok = _check_meta(meta, 'taxon_name', msg)

    # check we can find the taxon that the trait refers to
    if tx_ok and meta['taxon_name'] not in taxa and meta['taxon_name'] not in taxa_fields:
        msg.warn('Taxon name neither in the Taxa worksheet nor the name of a Taxa field', 2)

    # Regardless of the outcome of the meta checks, can still check the
    # data is all numeric, as it claims to be.
    is_numeric = [isinstance(vl, numbers.Number) for vl in data]
    if not all(is_numeric):
        msg.warn('Non numeric data found', 2)


# High level functions

def check_file(fname, verbose=True, force_entrez=False, locations_json=None):

    """
    Runs the format checking across an Excel workbook.

    Parameters:
        fname: Path to an Excel file
        verbose: Boolean to indicate whether to print messages as the program runs?
        force_entrez: Should the taxon checking use Entrez even when a local database is available
        locations_json: The path to a json file of valid location names

    Returns:
        A dictionary containing descriptive information about the workbook:
        - messages: The Messages instance
        - summary: The contents of the summary Worksheet metadata
    """

    try:
        workbook = openpyxl.load_workbook(filename=fname, data_only=True, read_only=True)
    except IOError:
        raise IOError('Could not open file {}'.format(fname))

    # now that we have a file, initialise the message tracker
    msg = Messages("Checking file '{}'".format(fname), verbose)

    # check the metadata sheets
    summary = get_summary(workbook, msg)
    locations = get_locations(workbook, msg, locations_json=locations_json)
    taxa, tax_index = get_taxa(workbook, msg, force_entrez=force_entrez)

    if 'data_worksheets' in summary and len(summary['data_worksheets']):
        for idx, data_ws in enumerate(summary['data_worksheets']):
            summary['data_worksheets'][idx] = check_data_worksheet(workbook, data_ws, taxa,
                                                                   locations, msg)
    else:
        msg.info('No data worksheets found')

    if msg.n_warnings:
        msg.info('FAIL: file contained {} errors'.format(msg.n_warnings))
    else:
        msg.info('PASS: file formatted correctly')

    return {'messages': msg, 'summary': summary, 'taxonomy': tax_index}


def main():

    """
    This program validates an Excel file formatted as a SAFE dataset. As it runs, it outputs
    a report that highlights any problems with the formatting.

    The program validates taxonomic names against the NCBI taxonomy database. It will use
    the ete2 python package to use a local version if possible, but will otherwise attempt
    to use the Entrez web service to validate names. The program also validate sampling
    location names: by default, this is loaded automatically from the SAFE website so requires
    an internet connection, but a local copy can be provided for offline use.
    """

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('fname', help="Path to the Excel file to be validated.")
    parser.add_argument('-l', '--locations_json',  default=None,
                        help=('Path to a locally stored json file of valid location names'))
    parser.add_argument('--force_entrez', action="store_true", default=False,
                        help=('Use entrez queries for taxon validation, even '
                              'if ete2 is available.'))

    args = parser.parse_args()

    check_file(fname=args.fname, verbose=True, locations_json=args.locations_json,
               force_entrez=args.force_entrez)


if __name__ == "__main__":
    main()

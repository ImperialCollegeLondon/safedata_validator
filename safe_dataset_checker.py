#!/usr/bin/python
# -*- coding: UTF8 -*-

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
from collections import Counter, OrderedDict
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
    if not os.path.exists(os.path.join(os.environ.get('HOME', '/'), '.etetoolkit', 'taxa.sqlite')):
        raise RuntimeError('ETE database not found')
    USE_ETE = True
except (ImportError, RuntimeError) as err:
    USE_ETE = False


# define some regular expressions used to check validity
RE_ORCID = re.compile(r'[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{4}')
RE_EMAIL = re.compile(r'\S+@\S+\.\S+')
RE_NAME = re.compile(r'[^,]+,[ ]?[^,]+')
RE_WSPACE_ONLY = re.compile(r'^\s*$')
# note that logically the next one also catches whitespace at both ends
RE_WSPACE_AT_ENDS = re.compile(r'^\s+\w+|\w+\s+$')
RE_DMS = re.compile(r'[°\'"dms’”]+')


"""
Define a class that carries a report about file processing
from function to function
"""


class Messenger(object):

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
            message += ', '.join([repr(str(vl)) for vl in join])
        elif join:
            message += ', '.join(join)

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


"""
Some simple helper functions that collect repeatedly used code
"""


def is_blank(value):
    """
    Helper function that checks if a value is None or contains only whitespace
    Args:
        value: A single value

    Returns:
        Boolean
    """

    return (value is None) or (RE_WSPACE_ONLY.match(str(value)) is not None)


def is_padded(value):
    """
    Helper function that checks if the value is padded with whitespace
    Args:
        value: Contents of a cell

    Returns:
        Boolean
    """

    return (value is not None) and (RE_WSPACE_AT_ENDS.match(value))


def duplication(data, msngr, message, depth):
    """
    Looks for duplication in a list of strings and prints a warning if
    duplication is found, returning the duplicated elements.

    Args:
        data: A list of strings
        msngr: A Messages instance
        message: A string for the warning message
        depth: Messages print depth for any warning raised

    Returns:
        A list of duplicated elements
    """

    # check for duplicate names
    duplicated = [ky for ky, vl in Counter(data).iteritems() if vl > 1]
    if duplicated:
        msngr.warn(message, depth, join=duplicated, as_repr=True)

    return duplicated

"""
Main functions to check contents of the file
"""


def get_summary(workbook, msngr):

    """
    Checks the information in the summary worksheet and looks for the metadata and
    dataset worksheets. The function is intended to try and get as much information
    as possible from the worksheet: the dictionary of metadata returned will have
    None for any missing data, which should be handled by downstream code.

    Args:
        workbook: An openpyxl Workbook instance
        msngr: A Messages instance
    Returns:
        A dictionary of the available summary metadata. The contents should be
        serialisable using simplejson, so that the summary can be stored in a
        database field.
    """

    # try and get the summary worksheet
    msngr.info("Checking Summary worksheet")
    start_warn = msngr.n_warnings
    try:
        worksheet = workbook['Summary']
    except KeyError:
        msngr.warn("Summary worksheet not found, moving on.")
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
        msngr.warn('Missing metadata fields: ', 1, join=required - found)

    # check for contents
    if any([set(x) == {None} for x in summary_dict.values()]):
        msngr.warn('Metadata fields with no information.', 1)

    # CHECK PROJECT ID
    if 'SAFE Project ID' not in summary_dict:
        msngr.warn('SAFE Project ID missing', 1)
        ret_dict['project_id'] = None
    else:
        if not isinstance(summary_dict['SAFE Project ID'][0], long):
            msngr.warn('SAFE Project ID is not an integer.', 1)
        ret_dict['project_id'] = summary_dict['SAFE Project ID'][0]

    # CHECK DATASET TITLE
    if 'Title' not in summary_dict:
        msngr.warn('Dataset title missing', 1)
        ret_dict['title'] = None
    else:
        ret_dict['title'] = summary_dict['Title'][0]

    # CHECK DATASET DESCRIPTION
    if 'Description' not in summary_dict:
        msngr.warn('Dataset title missing', 1)
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
                msngr.warn('Dataset embargoed but no date provided.', 1)
            embargo_date = summary_dict['Embargo date'][0]
            now = datetime.datetime.now()
            if not isinstance(embargo_date, datetime.datetime):
                msngr.warn('Embargo date not formatted as date.', 1)
            elif embargo_date < now:
                msngr.warn('Embargo date is in the past.', 1)
            elif embargo_date > now + datetime.timedelta(days=2*365):
                msngr.warn('Embargo date more than two years in the future.', 1)
            else:
                ret_dict['embargo_date'] = embargo_date.date().isoformat()
        elif access_status == 'Open':
            pass
        else:
            msngr.warn('Access status must be Open or Embargo '
                       'not {}'.format(repr(str(access_status))), 1)

    # CHECK KEYWORDS
    if 'Keywords' not in summary_dict:
        msngr.warn('Dataset keywords missing', 1)
        ret_dict['keywords'] = None
    else:
        # drop any blanks
        keywords = [vl for vl in summary_dict['Keywords'] if not is_blank(vl)]
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
                    msngr.warn('Author {} missing'.format(key), 1)
            # check validity
            if not RE_NAME.match(auth['name']):
                msngr.warn('Author name not formated as last_name, '
                           'first_names: {}'.format(auth['name']), 1)

            if auth['orcid'] is None:
                msngr.hint('Consider adding an ORCiD!', 1)
            elif not RE_ORCID.match(str(auth['orcid'])):
                msngr.warn('ORCID not properly formatted: {}'.format(auth['orcid']), 1)

            if not RE_EMAIL.match(auth['email']):
                msngr.warn('Email not properly formatted: {}'.format(auth['email']), 1)
            authors[ind] = auth

        ret_dict['authors'] = authors
    else:
        msngr.warn('Author metadata block incomplete.', 1)
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
                    msngr.warn('Data worksheet {} missing'.format(key), 1)
            # validate
            if data_ws['worksheet'] not in valid_sheets:
                msngr.warn('Dataset worksheet not found: {}.'.format(data_ws['worksheet']), 1)
            data_worksheets[ind] = data_ws

        ret_dict['data_worksheets'] = data_worksheets

        # check for extra undocumented spreadsheets
        if 'Worksheet name' in summary_dict:
            expected_sheets = set([vl['worksheet'] for vl in data_worksheets]
                                  + ['Summary', 'Taxa', 'Locations'])
            if valid_sheets != expected_sheets:
                msngr.warn('Undocumented sheets found in '
                           'workbook: ', 1, join=valid_sheets - expected_sheets)
    else:
        msngr.warn('Data worksheet metadata block incomplete.', 1)
        ret_dict['data_worksheets'] = None

    # summary of processing
    if (msngr.n_warnings - start_warn) > 0:
        msngr.info('Summary contains {} errors'.format(msngr.n_warnings - start_warn), 1)
    else:
        msngr.info('Summary formatted correctly', 1)

    return ret_dict


def get_locations(workbook, msngr, locations_json=None):

    """
    Attempts to load and check the contents of the Locations worksheet.
    Args:
        workbook: An openpyxl Workbook instance
        msngr: A Messages instance
        locations_json: A path to a JSON file containing a valid set of location names.
            With the default value of None, the function tries to get this from
            a SAFE project website service.
    Returns:
        A set of provided location names or None if the worksheet cannot
        be found or no locations can be loaded.
    """

    # try and get the locations worksheet
    msngr.info("Checking Locations worksheet")
    start_warn = msngr.n_warnings
    try:
        locs = workbook['Locations']
    except KeyError:
        # No locations is pretty implausible, but still persevere as if
        # they aren't going to be required
        msngr.warn("No locations worksheet found - moving on", 1)
        return None

    # get the set of valid names
    if locations_json is None:
        # If no file is provided then try and get locations from the website service
        loc_json = requests.get('https://www.safeproject.net/call/json/get_locations')
        if loc_json.status_code != 200:
            msngr.warn('Could not download valid location names. Use a local json file.', 1)
            valid_locations = []
        else:
            valid_locations = loc_json.json()['locations']
    else:
        # try and load the file
        try:
            loc_json = simplejson.load(file(locations_json))
            valid_locations = loc_json['locations']
        except IOError:
            msngr.warn('Could not load location names from file.', 1)
            valid_locations = []

    # Check the headers
    fields = []
    max_col = locs.max_column
    max_row = locs.max_row

    for col in range(1, max_col):
        fields.append(locs.cell(row=1, column=col).value)

    # check the key fields are there
    if not set(fields).issuperset({'Location name'}):
        msngr.warn('Location name column not found', 1)
        loc_names = None
    else:
        # get the location names as strings
        names_col = fields.index('Location name') + 1
        rows = range(2, max_row + 1)
        loc_names = [str(locs.cell(row=rw, column=names_col).value) for rw in rows]

        # report number of locations found
        msngr.info('Checking {} taxa'.format(len(loc_names)), 1)

        # check for duplicate names
        _ = duplication(loc_names, msngr, 'Duplicated location names: ', 1)

        # check for rogue whitespace and get representation of bad names
        ws_padded = [vl for vl in loc_names if is_padded(vl)]
        if ws_padded:
            msngr.warn('Locations names with whitespace padding: ', 1, join=ws_padded, as_repr=True)

        # check they are known - white strip first.
        loc_names = set([lc.strip() for lc in loc_names])
        unknown= loc_names - set(valid_locations)
        if unknown:
            msngr.warn('Unknown locations found: ', 1, join=unknown, as_repr=True)

    if (msngr.n_warnings - start_warn) > 0:
        msngr.info('Locations contains {} errors'.format(msngr.n_warnings - start_warn), 1)
    else:
        msngr.info('{} locations loaded correctly'.format(len(loc_names)), 1)

    return loc_names


def get_taxa(workbook, msngr, use_entrez=False, check_all_ranks=False):

    """
    Attempts to load and check the content of the Taxa worksheet. The
    function checks that:
    i)   all taxa have a taxon name and a taxon type,
    ii)  that taxa have a complete taxonomic hierarchy up to the taxon type
    iii) that all taxonomic names are known to NCBI, unless they are
         explicitly marked as new species with an asterisk suffix.

    Args:
        workbook: An openpyxl Workbook instance.
        msngr: A Messages instance.
        use_entrez: Use entrez queries even when a local NCBI query is available
        check_all_ranks: Validate all taxonomic ranks provided, not just the required ones
    Returns:
        A set of provided taxon names or None if the worksheet cannot
        be found or no taxa can be loaded.
    """

    if use_entrez:
        use_ete = False
    else:
        use_ete = USE_ETE

    # try and get the taxon worksheet
    msngr.info("Checking Taxa worksheet")
    start_warn = msngr.n_warnings
    try:
        sheet = workbook['Taxa']
    except KeyError:
        # This might mean that the study doesn't have any taxa, so return an empty
        # set. If the datasets then contain taxonomic names, it'll fail gracefully.
        msngr.hint("No taxa worksheet found - assuming no taxa in data for now!", 1)
        return set()

    # get and check the headers
    tx_rows = sheet.rows
    hdrs = [cl.value for cl in tx_rows.next()]
    # duplicated headers are a problem in that it will cause values in
    # the taxon dictionaries to be overwritten. We don't want to keep the
    # returned values, so discard
    _ = duplication(hdrs, msngr, 'Duplicated taxon sheet headers', 1)

    # Load dictionaries of the taxa and check some taxa are found
    taxa = [{ky: cl.value for ky, cl in zip(hdrs, rw)} for rw in tx_rows]

    # report number of taxa found
    if len(taxa) == 0:
        msngr.info('No taxon rows found'.format(len(taxa)), 1)
        return None

    # i) remove any values keyed to None
    # ii) convert keys to lower case
    # ii) update the list of headers
    msngr.info('Checking {} taxa'.format(len(taxa)), 1)
    _ = [tx.pop(None) for tx in taxa if None in tx]
    taxa = [{ky.lower(): vl for ky, vl in tx.iteritems()} for tx in taxa]
    hdrs = taxa[0].keys()

    # basic checks on taxon names: do they exist, no duplication, no padding
    if 'taxon name' not in hdrs:
        msngr.warn('No taxon name column found - no further checking', 1)
        return None

    tx_names = [tx['taxon name'] for tx in taxa]
    nm_empty = [str(ix + 2) for ix, vl in enumerate(tx_names) if is_blank(vl)]
    if nm_empty:
        msngr.warn('Taxon names blank in row(s): ', 1, join=nm_empty)

    _ = duplication(tx_names, msngr, 'Duplicated taxon names found: ', 1)

    nm_padded = [vl for vl in tx_names if is_padded(vl)]
    if nm_padded:
        msngr.warn('Taxon names with whitespace padding: ', 1, join=nm_padded, as_repr=True)

    # basic checks on taxon types: they exist and no padding
    if 'taxon type' not in hdrs:
        msngr.warn('No taxon type column found - no taxon verification', 1)
        return set(tx_names)

    tx_types = [tx['taxon type'] for tx in taxa]
    tp_padded = [vl for vl in tx_types if is_padded(vl)]
    if tp_padded:
        msngr.warn('Taxon types with whitespace padding: ', 1, join=tp_padded, as_repr=True)

    # Now check the taxonomic levels
    # - We only require the names for the big seven taxonomic levels
    rk_required = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    # - Users could provide any of these, which NCBI does know about
    #   but which are likely to be incompletely provided. The names key to the index
    #   of the required fields for that taxonomic level.
    # - Functional group and morphospecies are inserted here at the level we require
    # - Using an ordered dict to maintain hierarchy
    rk_known = OrderedDict([('kingdom', 0), ('functional group', 0), ('subkingdom', 0),
                           ('superphylum', 0), ('phylum', 1), ('subphylum', 1),
                           ('superclass', 1), ('class', 2), ('subclass', 2),
                           ('infraclass', 2), ('cohort', 2), ('superorder', 2),
                           ('order', 3), ('suborder', 3), ('infraorder', 3),
                           ('parvorder', 3), ('superfamily', 3), ('morphospecies', 3),
                           ('family', 4), ('subfamily', 4), ('tribe', 4),
                           ('subtribe', 4), ('genus', 5), ('subgenus', 5),
                           ('species group', 5), ('species subgroup', 5),
                           ('species', 6), ('subspecies', 6), ('varietas', 6), ('forma', 6)])

    # get the provided ranks
    rk_provided = set(hdrs) & set(rk_known)
    msngr.info('Fields provided for ' + str(len(rk_provided)) +
               ' taxonomic ranks: ' + ', '.join(rk_provided), 1)

    # get the taxon types and those types with no matching column
    tx_types = set([tp.lower() for tp in tx_types if tp is not None])
    tx_types_missing = tx_types - rk_provided

    # Check hierarchy for each taxon
    msngr.info('Checking completeness of taxonomic hierarchies', 1)
    for ridx, this_taxon in enumerate(taxa):
        tx_tp = this_taxon['taxon type']
        tx_nm = this_taxon['taxon name']
        tx_nm = '' if is_blank(tx_nm) else ', ' + tx_nm

        if is_blank(tx_tp):
            # Null taxon type
            msngr.warn('[R{}{}] Taxon type blank'.format(ridx + 2, tx_nm), 2)
        elif tx_tp.lower() in tx_types_missing:
            # Provided but no matching column
            msngr.warn('[R{}{}] Taxon type does not match to a column '
                       'name'.format(ridx + 2, tx_nm), 2)
        else:
            # Check there is some information up to the required level
            idx_required = rk_known[tx_tp.lower()]
            vals_required = [this_taxon[rnk] for rnk in rk_required[0:idx_required]
                             if rnk not in tx_types_missing]
            vals_empty = [is_blank(vl) for vl in vals_required]
            if any(vals_empty):
                txt = '[R{}{}] Taxon information not complete to {} level'
                msngr.warn(txt.format(ridx + 2, tx_nm, rk_required[idx_required]), 2)

    # Now setup to check the taxa against the NCBI database
    if use_ete:
        ncbi = NCBITaxa()
        msngr.info('Using local NCBI database to validate names', 1)
    else:
        # - search query to see if the name exists at the given rank
        entrez = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
                  'db=taxonomy&term={}+AND+{}[rank]&retmode=json')
        msngr.info('Using Entrez queries to validate names', 1)

    # Initialise a set to record when an Entrez query fails and a list to
    # compile the taxonomic index
    unvalidated = set()
    tx_index = []

    # which ranks to check - all or just required?
    # - don't try and check ranks that haven't been provided
    # - don't check morphospecies or functional groups against NCBI
    if check_all_ranks:
        rk_to_check = [rk for rk in rk_known if rk in rk_provided
                       and rk not in ['morphospecies', 'functional group']]
    else:
        rk_to_check = [rk for rk in rk_required if rk in rk_provided
                       and rk not in ['morphospecies', 'functional group']]

    for rnk in rk_to_check:

        msngr.info('Checking taxa at ' + rnk + ' level', 2)
        # Get the set of values that aren't None or whitespace at this rank
        tx_to_check = {tx[rnk] for tx in taxa if not is_blank(tx[rnk])}

        # Check for rogue whitespace
        ws_padded = set([tx for tx in tx_to_check if is_padded(tx)])
        if ws_padded:
            msngr.warn('Taxon name with whitespace padding: ', 3, join=ws_padded, as_repr=True)

        # Look for new taxa flagged as new
        new = {vl for vl in tx_to_check if vl.endswith('*')}
        if new:
            msngr.hint('New taxon reported: ' + ', '.join(new), 3)

        # drop new taxa
        tx_to_check -= new
        # tidy padded taxa to check for them now
        tx_to_check = {tx.strip() for tx in tx_to_check}

        # This merges taxa that aren't found at all or  which aren't found at the
        # stated rank - this is marginally less clear for users but two entrez
        # queries are needed to discriminate, so keep the mechanism simple.
        if use_ete:
            # look up the names
            ids = ncbi.get_name_translator(tx_to_check)
            # isolate names that aren't found
            not_found = tx_to_check - set(ids.keys())
            # check the ranks
            ranks_found = {ky: ncbi.get_rank(vl).values() for ky, vl in ids.iteritems()}
            bad_rank = set([ky for ky, vl in ranks_found.iteritems() if rnk not in vl])
            invalid = not_found | bad_rank
        else:
            # loop the names through the entrez interface to the taxonomy database
            invalid = set()
            # loop over the values
            for each_tx in tx_to_check:
                entrez_response = requests.get(entrez.format(each_tx, rnk))
                if entrez_response.status_code != 200:
                    # internet failures just add individual taxa to the unvalidated list
                    unvalidated.add(each_tx)
                else:
                    if entrez_response.json()['esearchresult']['count'] == '0':
                        invalid.add(each_tx)

        # d) Report invalid taxa
        if invalid:
            msngr.warn('Taxa not found or not valid at this rank: ', 3, join=invalid)

        # build taxon index - currently no validation
        tx_index.extend([(tx, rnk) for tx in tx_to_check])

    # report on unvalidated taxa
    if unvalidated:
        msngr.warn('Entrez taxon validation failed for: ', 2, join=unvalidated)

    if (msngr.n_warnings - start_warn) > 0:
        msngr.info('Taxa contains {} errors'.format(msngr.n_warnings - start_warn), 1)
    else:
        msngr.info('{} taxa loaded correctly'.format(len(tx_names)), 1)

    return set(tx_names), tx_index


def check_data_worksheet(workbook, ws_meta, taxa, locations, msngr):

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
        msngr: A Messages instance

    Returns:
        A possibly updated version of ws_meta.
    """

    ws_name = ws_meta['worksheet']
    if ws_name not in workbook.get_sheet_names():
        msngr.warn('Data worksheet {} not found'.format(ws_name))
        return ws_meta
    else:
        msngr.info('Checking data worksheet {}'.format(ws_name))
        start_warn = msngr.n_warnings

    # get the worksheet and data dimensions
    worksheet = workbook[ws_name]
    max_col = worksheet.max_column
    max_row = worksheet.max_row

    # trap completely empty worksheets
    if max_row == 1:
        msngr.warn('Worksheet is empty', 1)
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
        msngr.warn('Cannot parse data: field_name row not found', 1)
        return ws_meta

    # Neither of the row and col maxima are particularly reliable as Excel can hang on to
    # cell references for previously used cells. We can ignore blank columns easily but
    # we do need to know where to actually stop for finding blank data in rows.
    # So, we explicitly check the row numbers, making them a mandatory part of the setup.

    # - get the values
    row_number_cells = worksheet.get_squared_range(1, field_name_row + 1, 1, max_row)
    row_numbers = [cl[0].value for cl in row_number_cells]

    # - trim blank or whitespace values from the end and update the max col
    while is_blank(row_numbers[-1]):
        row_numbers.pop()

    max_row = len(row_numbers) + field_name_row

    # now check the row numbers are numbers and if they are
    # do they start at one and go up by one
    if not all([isinstance(vl, numbers.Number) for vl in row_numbers]):
        msngr.warn('Non-numeric data found in row numbering', 1)
    else:
        if row_numbers[0] != 1:
            msngr.warn('Row numbering does not start at 1', 1)

        one_increment = [(vl1 - vl2) == 1 for vl1, vl2 in zip(row_numbers[1:], row_numbers[:-1])]
        if not all(one_increment):
            msngr.warn('Row numbering does not consistently increment by 1', 1)

    # report on detected size
    msngr.info('Worksheet contains {} rows and {} columns'.format(max_row, max_col), 1)
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
        msngr.warn('Categorical data fields found but no levels descriptor provided.', 1)

    if 'Numeric' in ft_found:
        if 'units' not in descriptors:
            msngr.warn('Numeric data fields found but no units descriptor provided.', 1)
        if 'method' not in descriptors:
            msngr.warn('Numeric data fields found but no method descriptor provided.', 1)

    if 'Location' in ft_found and locations is None:
        msngr.warn('Location field found but no Location worksheet provided.', 1)

    if 'Abundance' in ft_found:
        if 'taxon_name' not in descriptors:
            msngr.warn('Abundance field found but no taxon name descriptor provided.', 1)
        if 'method' not in descriptors:
            msngr.warn('Abundance field found but no sampling method descriptor provided.', 1)

    # check field names unique (drop None)
    field_names = [fld['field_name'] for fld in field_metadata if fld['field_name'] is not None]
    _ = duplication(field_names, msngr, 'Field names duplicated: ', 1)

    # get taxa field names for cross checking observation and trait data
    taxa_fields = [fld['field_name'] for fld in field_metadata if fld['field_type'] == 'Taxa']

    # TODO - check mandatory fields

    # check each field
    for idx, meta in enumerate(field_metadata):

        # skip any field with no metadata
        if set(meta.values()) == {None}:
            break

        # prep the messages instance to pass to functions
        if is_blank(meta['field_name']):
            msngr.info('Checking Column {}'.format(utils.get_column_letter(idx + 2)), 1)
            msngr.warn('Field name is blank', 2)
        else:
            msngr.info('Checking field {field_name}'.format(**meta), 1)

        # check the description
        if is_blank(meta['description']):
            msngr.warn('Description is missing', 2)

        # read the values
        data_block = worksheet.get_squared_range(idx + 2, field_name_row + 1, idx + 2, max_row)
        data = [cl[0].value for cl in data_block]

        # filter out missing and blank data, except for comments fields, where
        # blanks are not an error
        if meta['field_type'] != 'Comments':
            data = filter_missing_or_blank_data(data, msngr)

        # run consistency checks where needed and trap unknown field types
        if meta['field_type'] == 'Date':
            check_field_date(data, msngr)
        elif meta['field_type'] == 'Datetime':
            check_field_datetime(data, msngr)
        elif meta['field_type'] == 'Time':
            check_field_time(data, msngr)
        elif meta['field_type'] == 'Taxa':
            check_field_taxa(data, taxa, msngr)
        elif meta['field_type'] == 'Location':
            # location names should be strings
            data = [str(dt) for dt in data]
            check_field_locations(data, locations, msngr)
        elif meta['field_type'] == 'Categorical':
            check_field_categorical(meta, data, msngr)
        elif meta['field_type'] == 'Numeric':
            check_field_numeric(meta, data, msngr)
        elif meta['field_type'] == 'Abundance':
            check_field_abundance(meta, data, taxa, taxa_fields, msngr)
        elif meta['field_type'] == 'Trait':
            check_field_trait(meta, data, taxa, taxa_fields, msngr)
        elif meta['field_type'] == 'Latitude':
            check_field_geo(meta, data, msngr, which='latitude')
        elif meta['field_type'] == 'Longitude':
            check_field_geo(meta, data, msngr, which='longitude')
        elif meta['field_type'] == 'Comments':
            pass
        elif meta['field_type'] is None:
            msngr.warn('Field type is empty', 2)
        else:
            msngr.warn('Unknown field type {field_type}'.format(**meta), 2)

    if (msngr.n_warnings - start_warn) > 0:
        msngr.info('Dataframe contains {} errors'.format(msngr.n_warnings - start_warn), 1)
    else:
        msngr.info('Dataframe formatted correctly', 1)

    # trim out field descriptions for blank columns
    field_metadata = [fld for fld in field_metadata if set(fld.values()) != {None}]

    # update ws_meta with the field information
    ws_meta['fields'] = field_metadata

    return ws_meta


"""
Helper functions for checking data worksheets
"""


def filter_missing_or_blank_data(data, msngr):

    """
    Takes a list of data and filters out any missing or blank data,
    reporting to the Messages instance as it goes. The filtered list
    is returned for feeding to field type checker functions so that
    they can assume no NA or blank data in their checks.

    Args:
        data: A list of values read from the Worksheet.
        msngr: A Messages instance
    Returns:
        A list of data values filtered to remove NA and blanks.
    """

    # Only NA is acceptable
    na_vals = [vl == u'NA' for vl in data]
    if any(na_vals):
        msngr.hint('{} / {} values missing'.format(sum(na_vals), len(na_vals)), 2)

    # We won't tolerate:
    # 1) empty cells (just to avoid ambiguity - e.g. in abundance data)
    # 2) non-empty cells containing only whitespace strings
    blank = [is_blank(vl) for vl in data]
    if any(blank):
        msngr.warn('{} cells are blank or contain only whitespace text'.format(sum(blank)), 2)

    # Return the values that aren't NA, blank or whitespace only
    na_or_blank = [any(tst) for tst in zip(na_vals, blank)]
    data = [dt for dt, nb in zip(data, na_or_blank) if not nb]

    return data


def _check_meta(meta, descriptor, msngr):
    """
    A standardised check to see if a required descriptor is present for
    a field and that it isn't simply empty or whitespace. The function
    reports problems to the Messages instance and returns a boolean
    showing if the checks passed successfully.

    Args:
        meta: A dictionary of field metadata descriptors
        descriptor: The name of the descriptor to check.
        msngr: An instance of class Messages

    Returns:
        A boolean, with True showing no problems and False showing
        that warnings occurred.
    """

    if descriptor not in meta:
        msngr.warn('{} descriptor missing'.format(descriptor), 2)
        return False
    elif is_blank(meta[descriptor]):
        msngr.warn('{} descriptor is blank'.format(descriptor), 2)
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


def all_numeric(data):
    """
    Tests if all values in an iterable are numeric.
    Args:
        data: An iterable claiming to contain numbers

    Returns:
        Boolean
    """

    return all([isinstance(vl, numbers.Number) for vl in data])


def check_field_date(data, msngr):

    """
    Checks for data consistency in date fields and reports to the
    Messages instance.

    Args:
        data: A list of data values, allegedly of type datetime.date
        msngr: A Messages instance
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data])
    if type_check != {datetime.datetime}:
        msngr.warn('Non-date data in field.', 2)
        if datetime.time in type_check:
            msngr.hint('Some values _only_  contain time components', 2)

    # Check no time component in actual dates
    no_time = [vl.time() == datetime.time(0, 0) for vl in data]
    if not all(no_time):
        msngr.warn('Some values also contain time components', 2)


def check_field_datetime(data, msngr):

    """
    Checks for data consistency in datetime fields and reports to the
    Messages instance.

    Args:
        data: A list of data values, allegedly of type datetime.datetime
        msngr: A Messages instance
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data])
    if type_check != {datetime.datetime}:
        msngr.warn('Non-date data in field.', 2)
        if datetime.time in type_check:
            msngr.hint('Some values _only_  contain time components', 2)


def check_field_time(data, msngr):

    """
    Checks for data consistency in time fields and reports to the
    Messages instance.

    Args:
        data: A list of data values, allegedly of type datetime.time
        msngr: A Messages instance
    """

    # Check type (excluding NA values)
    type_check = set([type(vl) for vl in data])
    if type_check != {datetime.time}:
        msngr.warn('Non-time formatted data found.', 2)


def check_field_taxa(data, taxa, msngr):

    """
    Checks if all the values provided in a Taxon field are found
    in the Taxa worksheet, reporting to the Messages instance.

    Args:
        data: A list of data values, allegedly taxon names
        taxa: A set containing taxon names from the Taxa worksheet
        msngr: A Messages instance
    """

    found = set(data)
    if taxa is None:
        msngr.warn('Taxa worksheet not provided or no taxon names were found', 2)
    if not found.issubset(taxa):
        msngr.warn('Includes taxa missing from Taxa worksheet: ', 2,
                   join=found - taxa, as_repr=True)


def check_field_locations(data, locations, msngr):

    """
    Checks if all the values provided in a Locations field are
    found in the Locations worksheet, reporting to the Messages instance.

    Args:
        data: A list of data values, allegedly taxon names
        locations: A set containing locations from the Locations worksheet
        msngr: A Messages instance
    """

    # check if locations are all provided
    found = set(data)
    if locations is None:
        msngr.warn('No Locations worksheet provided', 2)
    elif not found.issubset(locations):
        msngr.warn('Includes locations missing from Locations worksheet:', 2,
                   join=found - locations, as_repr=True)


def check_field_abundance(meta, data, taxa, taxa_fields, msngr):

    """
    Checks abundance type data, reporting to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly numeric
        taxa: A set containing taxon names from the Taxa worksheet
        taxa_fields: A list of Taxa fields in this worksheet.
        msngr: A Messages instance
    """

    # check the required descriptors
    _check_meta(meta, 'method', msngr)
    tx_ok = _check_meta(meta, 'taxon_name', msngr)

    if tx_ok and meta['taxon_name'] not in taxa and meta['taxon_name'] not in taxa_fields:
        msngr.warn('Taxon name neither in the Taxa worksheet nor the name of a Taxa field', 2)

    # Can still check values are numeric, whatever happens above.
    # We're not going to insist on integers here - could be mean counts.
    if not all_numeric(data):
        msngr.warn('Field contains non-numeric data', 2)


def check_field_categorical(meta, data, msngr):

    """
    Checks factor data, reporting to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly strings from a provided set of levels
        msngr: A Messages instance
    """

    # this has already been tested but make it robust
    ct_ok = _check_meta(meta, 'levels', msngr)

    if not ct_ok:
        # Can't check further if no levels descriptor
        pass
    elif ct_ok and not isinstance(meta['levels'], unicode):
        # Can't really check anything here either
        msngr.warn('Category description does not seem to be text', 2)
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
            msngr.warn('Integer level names not permitted', 2)

        # Now look for consistency: get the unique values reported in the
        # data, convert to unicode to handle checking of integer labels and
        # then check the reported levels are a subset of the descriptors.
        reported = set(data)
        reported = {unicode(lv) for lv in reported}

        if not reported.issubset(level_labels):
            msngr.warn('Categories found in data missing from '
                       'description: ', 2, join= reported - level_labels)


def check_field_numeric(meta, data, msngr):

    """
    Checks numeric type data, reporting to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly numeric
        msngr: A Messages instance
    """
    _check_meta(meta, 'units', msngr)
    _check_meta(meta, 'method', msngr)

    # Regardless of the outcome of the meta checks, can still check the
    # data is all numeric, as it claims to be.
    if not all_numeric(data):
        msngr.warn('Non numeric data found', 2)


def check_field_trait(meta, data, taxa, taxa_fields, msngr):

    """
    Checks trait type data - things measured on an organism - and
    reports to the Messages instance.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values, allegedly numeric
        taxa: A set containing taxon names from the Taxa worksheet
        taxa_fields: A list of Taxa fields in this worksheet.
        msngr: A Messages instance
    """

    _check_meta(meta, 'units', msngr)
    _check_meta(meta, 'method', msngr)
    tx_ok = _check_meta(meta, 'taxon_name', msngr)

    # check we can find the taxon that the trait refers to
    if tx_ok and meta['taxon_name'] not in taxa and meta['taxon_name'] not in taxa_fields:
        msngr.warn('Taxon name neither in the Taxa worksheet nor the name of a Taxa field', 2)

    # Regardless of the outcome of the meta checks, can still check the
    # data is all numeric, as it claims to be.
    if not all_numeric(data):
        msngr.warn('Non numeric data found', 2)

def check_field_geo(meta, data, msngr, which='latitude'):

    """
    Checks geographic coordinates.

    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values
        msngr: A Messages instance
        which: One of latitude or longitude
    """

    # Are the values represented as decimal degrees - numeric.
    if not all_numeric(data):
        msngr.warn('Non numeric data found', 2)
        if any([RE_DMS.search(unicode(vl)) for vl in data]):
            msngr.hint('Possible degrees minutes and seconds formatting? Use decimal degrees', 2)

    # Check the bounds
    nums = [vl for vl in data if isinstance(vl, numbers.Number)]
    if len(nums):
        if which == 'latitude':
            if min(nums) < -90 or max(nums) > 90:
                msngr.warn('Latitude values outside [-90, 90]',2)
            if min(nums) < -4 or max(nums) > 8:
                msngr.hint('Latitude values not in Borneo [-4, 8]',2)
        elif which == 'longitude':
            if min(nums) < -180 or max(nums) > 180:
                msngr.warn('Longitude values outside [-180, 180]', 2)
            if min(nums) < 108 or max(nums) > 120:
                msngr.hint('Longitude values not in Borneo [-108, 120]', 2)


# High level functions

def check_file(fname, verbose=True, use_entrez=False,
               locations_json=None, check_all_ranks=False):

    """
    Runs the format checking across an Excel workbook.

    Parameters:
        fname: Path to an Excel file
        verbose: Boolean to indicate whether to print messages as the program runs?
        use_entrez: Should the taxon checking use Entrez even when a local database is available
        check_all_ranks: Should all provided taxonomic ranks be validated?
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
    msngr = Messenger("Checking file '{}'".format(fname), verbose)

    # check the metadata sheets
    summary = get_summary(workbook, msngr)
    locations = get_locations(workbook, msngr, locations_json=locations_json)
    taxa, tax_index = get_taxa(workbook, msngr, use_entrez=use_entrez,
                               check_all_ranks=check_all_ranks)

    if 'data_worksheets' in summary and len(summary['data_worksheets']):
        for idx, data_ws in enumerate(summary['data_worksheets']):
            summary['data_worksheets'][idx] = check_data_worksheet(workbook, data_ws, taxa,
                                                                   locations, msngr)
    else:
        msngr.info('No data worksheets found')

    if msngr.n_warnings:
        msngr.info('FAIL: file contained {} errors'.format(msngr.n_warnings))
    else:
        msngr.info('PASS: file formatted correctly')

    return {'messages': msngr, 'summary': summary, 'taxonomy': tax_index}


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
    parser.add_argument('-l', '--locations_json', default=None,
                        help='Path to a locally stored json file of valid location names')
    parser.add_argument('--use_entrez', action="store_true", default=False,
                        help=('Use entrez queries for taxon validation, even '
                              'if ete2 is available.'))
    parser.add_argument('--check_all_ranks', action="store_true", default=False,
                        help=('Check the validity of all taxonomic ranks included, '
                              'not just the standard required ranks.'))

    args = parser.parse_args()

    check_file(fname=args.fname, verbose=True, locations_json=args.locations_json,
               use_entrez=args.use_entrez, check_all_ranks=args.check_all_ranks)


if __name__ == "__main__":
    main()

import requests
from enforce_typing import enforce_types
import datetime
import re
from safedata_validator.logger import LOGGER, FORMATTER, CH, log_and_raise
from safedata_validator.validators import (GetDataFrame, IsLower, IsNotBlank,
                                           IsNotPadded, IsNotSpace, IsString,
                                           NoPunctuation, HasDuplicates)
from safedata_validator.extent import Extent

# Compile some global regex expressions
RE_DOI = re.compile(r'https?://(dx.)?doi.org/')
RE_ORCID = re.compile(r'[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]')
RE_EMAIL = re.compile(r'\S+@\S+\.\S+')
RE_NAME = re.compile(r'[^,]+,[ ]?[^,]+')


class Summary:

    # Class instance to define the metadata that may be present.
    # These are defined in blocks of related rows described in 4-tuples:
    #  0: a list of tuples giving the rows in the block:
    #     * row header key
    #     * mandatory within the block
    #     * 'internal' name - also maps to Zenodo fields
    #     * accepted types
    #  1: is the block mandatory (bool)
    #  2: the title of the block for the logger
    #  3: should there be only one record (bool)
    #
    # Each entry in the list of rows provides the row header that appears
    # in the summary sheet, if that field is mandatory within the set and
    # an internal field name, if needed in the code or by Zenodo.

    fields = dict(core=([('safe project id', True, 'pid', int),
                         ('access status', True, 'access', str),
                         ('embargo date', False, 'embargo_date', datetime.datetime),
                         ('access conditions', False, 'access_conditions', str),
                         ('title', True, None, str),
                         ('description', True, None, str)],
                        True, 'Core fields', True),
                  keywords=([('keywords', True, None, str)],
                            True, 'Keywords', False),
                  doi=([('publication doi', True, None, str)],
                       False, 'DOI', False),
                  date=([('start date', True, None, datetime.datetime),
                         ('end date', True, None, datetime.datetime)],
                        False, 'Date Extents', True),
                  geo=([('west', True, None, float),
                        ('east', True, None, float),
                        ('south', True, None, float),
                        ('north', True, None, float)],
                       False, 'Geographic Extents', True),
                  authors=([('author name', True, 'name', str),
                            ('author affiliation', False, 'affiliation', str),
                            ('author email', False, 'email', str),
                            ('author orcid', False, 'orcid', str)],
                           True, 'Authors', False),
                  funding=([('funding body', True, 'body', str),
                            ('funding type', True, 'type', str),
                            ('funding reference', False, 'ref', (str, int, float)),
                            ('funding link', False, 'url', str)],
                           False, 'Funding Bodies', False),
                  external=([('external file', True, 'file', str),
                             ('external file description', True, 'description', str)],
                            False, 'External Files', False),
                  worksheet=([('worksheet name', True, 'name', str),
                              ('worksheet title', True, 'title', str),
                              ('worksheet description', True, 'description', str),
                              ('worksheet external file', False, 'external', str)],
                             False, 'Worksheets', False),
                  permits=([('permit type', True, 'type', str),
                            ('permit authority', True, 'authority', str),
                            ('permit number', True, 'number', (str, int, float))],
                           False, 'Permits', False))

    def __init__(self, worksheet, validate_doi=False, valid_pid=None):

        """
        Checks the information in the summary worksheet and looks for the
        metadata and dataset worksheets. The function is intended to try and
        get as much information as possible from the worksheet: the dictionary
        of metadata returned will have None for any missing data, which should
        be handled by downstream code.

        Args:
            worksheet: An openpyxl worksheet instance.
            validate_doi: Check any publication DOIs, requiring a web connection.
            valid_pid: If provided, an integer or list of integer values that are
                permitted in the Project ID field (usually one for a new dataset
                but more if a published dataset is associated with multiple projects
                and any of those ids would be valid).
        """

        self.worksheet = worksheet
        self.validate_doi = validate_doi

        self.project_id = None
        self.title = None
        self.description = None
        self.access = None
        self.access_conditions = None
        self.embargo_date = None
        self.authors = None
        self.publication_doi = None
        self.keywords = None
        self.temporal_extent = Extent('temporal extent', datetime.date)
        self.latitudinal_extent = Extent('latitudinal extent', float)
        self.longitudinal_extent = Extent('longitudinal extent', float)
        self._rows = None
        self._ncols = None
        self._valid_pid = []

        # validate project_id is one of None, an integer or a list of integers
        if valid_pid is None:
            pass
        elif isinstance(valid_pid, int):
            self._valid_pid = [valid_pid]
        elif isinstance(valid_pid, list):
            if not all([isinstance(pid, int) for pid in valid_pid]):
                LOGGER.error("Invalid value in list of project_ids.")
                self._valid_pid = None
            else:
                self._valid_pid = valid_pid
        else:
            LOGGER.error("Provided project id must be an integer or list of integers")
            self._valid_pid = None

    def load(self, worksheet, validate_doi=False, valid_pid=None):
        """

        Args:
            worksheet:
            validate_doi:
            valid_pid:

        Returns:

        """
        # try and get the summary worksheet
        LOGGER.info("Checking Summary worksheet")
        FORMATTER.push()

        start_errors = CH.counters['ERROR']

    def _load_rows(self):

        # load worksheet rows
        rows = list(self.worksheet.iter_rows(values_only=True))
        self._ncols = self.worksheet.max_column

        # convert into dictionary using the lower cased first entry as the key
        row_headers = IsString([r[0] for r in rows])
        if not row_headers:
            LOGGER.error('Summary metadata fields column contains non text values :',
                         extra={'join': row_headers.failed})

        self._rows = {rw[0].lower(): rw[1:] for rw in rows}

        # Check the minimal keys are expected - mandatory fields in mandatory blocks
        required_blocks = (blk[0] for blk in list(self.fields.values()) if blk[1])
        required = {fld[0] for blk in required_blocks for fld in blk if fld[1]}

        found = set(self._rows.keys())

        if not found.issuperset(required):
            LOGGER.error('Missing mandatory metadata fields: ',
                         extra={'join': required - found})

        # Check only valid keys are found
        valid_blocks = (blk[0] for blk in list(self.fields.values()))
        valid_fields = {fld[0] for blk in valid_blocks for fld in blk}

        if found - valid_fields:
            LOGGER.error('Unknown metadata fields: ',
                         extra={'join': found - valid_fields})

    def _read_block(self, field_desc, mandatory, title, only_one):
        """
        Internal function that takes a given block definition from
        SUMMARY_FIELDS and returns a list of dictionary records. This
        function automatically does some common checking for missing
        data, bad input types etc, leaving the block specific functions
        to handle unique tests.

        Args:
            field_desc:  A list of tuples describing fields.
            mandatory: Is the block mandatory?
            title: The display title for the block
            only_one: Are multiple records for the field an error?
        Returns:
            None or a list of records.
        """

        mandatory_fields = [f[0] for f in field_desc if f[1]]
        optional_fields = [f[0] for f in field_desc if not f[1]]
        field_map = [(f[0], f[2]) for f in field_desc]
        field_types = {f[0]: f[3] for f in field_desc}

        # Get the full set of field names
        all_fields = mandatory_fields + optional_fields

        # Get the data, filling in completely missing rows
        block = {k: self._rows[k] if k in self._rows else [None] * (self._ncols - 1)
                 for k in all_fields}

        # Empty cells are already None, but also filter values to catch
        # pure whitespace content and replace with None
        for ky, vals in block.items():
            vals = IsNotSpace(vals)
            if not vals:
                LOGGER.error(f'Whitespace only cells in field {ky}')

            block[ky] = vals.values

        # Pivot to dictionary of records
        block = [dict(zip(block.keys(), vals)) for vals in zip(*block.values())]

        # Drop empty records
        block = [bl for bl in block if any(bl.values())]

        if not block:
            if mandatory:
                LOGGER.error(f'No {title} metadata found')
                return None
        else:
            LOGGER.info(f'Metadata for {title} found: {len(block)} records')
            FORMATTER.push()

            if len(block) > 1 and only_one:
                LOGGER.error('Only a single record should be present')

            # report on block fields
            for fld in mandatory_fields:
                fld_values = [rec[fld] for rec in block]
                if not all(fld_values):
                    LOGGER.error(f'Missing metadata in mandatory field {fld}')

            # report on actual data that is of the wrong type
            for fld in all_fields:
                bad_values = [rec[fld] for rec in block
                              if rec[fld] is not None and
                              not isinstance(rec[fld], field_types[fld])]

                if bad_values:
                    LOGGER.error(f'Field {fld} contains values of wrong type: ',
                                 extra={'join': bad_values})

            # remap names if provided
            to_rename = (mp for mp in field_map if mp[1] is not None)
            for (old, new) in to_rename:
                for rec in block:
                    rec[new] = rec[old]
                    rec.pop(old)

            return block

    def _load_authors(self):

        authors = self._read_block(*self.fields['authors'])

        # Author specific validation
        if authors is not None:

            # Badly formatted names
            bad_names = [rec['name'] for rec in authors
                         if isinstance(rec['name'], str) and not RE_NAME.match(rec['name'])]
            if bad_names:
                LOGGER.error('Author names not formatted as last_name, first_names: ',
                             extra={'join': bad_names})

            # Emails not formatted properly
            bad_emails = [rec['email'] for rec in authors
                          if isinstance(rec['email'], str) and not RE_EMAIL.match(rec['email'])]
            if bad_emails:
                LOGGER.error('Author emails not properly formatted: ',
                             extra={'join': bad_emails})

            # ORCIDs not strings
            bad_orcid = [rec['orcid'] for rec in authors
                         if isinstance(rec['orcid'], str) and not RE_ORCID.match(rec['orcid'])]
            if bad_orcid:
                LOGGER.error('Author ORCIDs not properly formatted: ',
                             extra={'join': bad_orcid})

        self.authors = authors

    def _load_keywords(self):

        keywords = self._read_block(*self.fields['keywords'])

        # extra data validation for keywords
        if keywords:
            keywords = NoPunctuation(keywords)
            if not keywords:
                LOGGER.error('Put each keyword in a separate cell, do not separate '
                             'keywords using commas or semi-colons')

            self.keywords = keywords.values

    def _load_permits(self):

        # LOOK FOR PERMIT DETAILS - users provide a permit authority, number and permit type

        permits = self._read_block(*self.fields['permits'])

        # Permit specific checking for allowed permit types
        if permits:
            permit_types = [rec['type'].lower() for rec in permits if isinstance(rec['type'], str)]
            valid_permit_types = {'research', 'export', 'ethics'}
            if not set(permit_types).issubset(valid_permit_types):
                LOGGER.error('Unknown permit types: ',
                             extra={'join': permit_types - valid_permit_types})

        self.permits = permits

    def _load_doi(self):

        # CHECK FOR PUBLICATION DOIs
        pub_doi = self._read_block(*self.fields['doi'])

        # Extra data validation for DOIs
        if pub_doi is not None:

            # Check DOI URLS _are_ urls
            pub_doi_re = [RE_DOI.search(v['publication doi']) for v in pub_doi
                          if isinstance(v['publication doi'], str)]
            if not all(pub_doi_re):
                LOGGER.error('Publication DOIs not all in format: https://doi.org/...')

            if self.validate_doi:
                for is_doi in pub_doi_re:
                    if is_doi:
                        api_call = f'https://doi.org/api/handles/{is_doi.string[is_doi.end():]}'
                        r = requests.get(api_call)
                        if r.json()['responseCode'] != 1:
                            LOGGER.error(f'DOI not found: {is_doi.string}')

        self.publication_doi = pub_doi

    def _load_funders(self):

        # LOOK FOR FUNDING DETAILS - users provide a funding body and a description
        # of the funding type and then optionally a reference number and a URL

        funders = self._read_block(*self.fields['funding'])

        # TODO - currently no check beyond _read_block but maybe actually check
        #        the URL is a URL and maybe even opens? Could use urllib.parse

        self.funders = funders

    def _load_temporal_extent(self):

        temp_extent = self._read_block(*self.fields['date'])

        # temporal extent validation and updating
        if temp_extent is not None:

            start_date = temp_extent[0]['start date']
            end_date = temp_extent[0]['end date']

            if (isinstance(start_date, datetime.datetime) and
                    isinstance(end_date, datetime.datetime) and
                    start_date > end_date):
                LOGGER.error('Start date is after end date')
            else:
                self.temporal_extent.update([start_date, end_date])

    def _load_geographic_extent(self):

        # Geographic extents
        geo_extent = self._read_block(*self.fields['geo'])
        bbox = geo_extent[0]

        if all([isinstance(v, float) for v in bbox]):
            self.latitudinal_extent.update([bbox['south'], bbox['north']])
            self.longitudinal_extent.update([bbox['west'], bbox['east']])


    def _unused_blocks(self):


        # Now check core rows
        singletons = read_block(*fields['core'])
        singletons = singletons[0]

        # Project ID specific validation
        pid = singletons['pid']

        # Check the value is in the provided list
        if project_id is not None and pid not in project_id:
            pid_str = ', '.join([str(p) for p in project_id])
            LOGGER.error(f'SAFE Project ID in file ({pid}) does not match any '
                         f'provided project ids ({pid_str})')
        else:
            self.project_id = pid

        # Title validation
        if singletons['title'] is None:
            LOGGER.error('Dataset title is blank')
        else:
            self.title = singletons['title']

        # Description validation
        if singletons['description'] is not None:
            LOGGER.error('Dataset description is blank')
        else:
            self.description = singletons['description'].value

        # Access status and embargo validation
        if singletons['access'] is not None:
            access = singletons['access']
            if not isinstance(access, str):
                LOGGER.error(f'Access status not a text value: {access}')
            elif access.lower() not in ['open', 'embargo', 'restricted']:
                LOGGER.error('Access status must be Open, Embargo or Restricted '
                             f'not {access}')
            else:
                self.access = access.lower()

                if self.access == 'embargo':
                    embargo_date = singletons['embargo_date']
                    if embargo_date is None:
                        LOGGER.error('Dataset embargoed but no embargo date provided')
                    else:

                        if not isinstance(embargo_date, datetime.datetime):
                            LOGGER.error(f'Embargo date not a date value: {embargo_date }')
                        else:
                            now = datetime.datetime.now()

                            if embargo_date < now:
                                LOGGER.error('Embargo date is in the past.')
                            elif embargo_date > now + datetime.timedelta(days=2 * 365):
                                LOGGER.error('Embargo date more than two years in the future.')
                            else:
                                LOGGER.info(f'Dataset access: embargoed until {embargo_date }')
                                self.embargo_date = embargo_date.date().isoformat()
                elif self.access == 'restricted':
                    access_conditions = singletons['access_conditions']

                    if access_conditions is None:
                        LOGGER.error('Dataset restricted but no access conditions specified')
                    else:
                        if not isinstance(access_conditions, str):
                            LOGGER.error(f'Access conditions are not text: {access_conditions}')
                        else:
                            LOGGER.info(f'Dataset access: restricted with conditions {access_conditions}')
                            self.access_conditions = access_conditions
                else:
                    LOGGER.info(f'Dataset access: {self.access}')



        # PROCESS BLOCKS OF CELLS
        # Load the AUTHORS block




        # CHECK EXTENTS - there are six fields to provide temporal and geographic
        # extents for the data. These two extents are a mandatory component of Gemini
        # Metadata so need to be provided if the worksheets and locations don't
        # populate the extents.



        # Geographic extents
        geo_extent = read_block(fields['geo'])

        if geo_extent is not None:

            bbox = geo_extent[0]

            bbox_field_types = {None if val is None else val.ctype for val in list(bbox.values())}

            if not bbox_field_types.issubset({xlrd.XL_CELL_NUMBER, None}):
                LOGGER.error('Geographic extents not numeric')
            elif all(bbox.values()):

                strip_ctypes(geo_extent)
                valid = self.validate_geo_extent((bbox['west'], bbox['east']), 'longitude')
                if valid:
                    self.update_extent((bbox['west'], bbox['east']), float, 'longitudinal_extent')

                valid = self.validate_geo_extent((bbox['south'], bbox['north']), 'latitude')
                if valid:
                    self.update_extent((bbox['south'], bbox['north']), float, 'latitudinal_extent')

        # Finally, check the information on data content
        LOGGER.info('Checking data: worksheets and external files',
                    extra={'indent_before': 1})

        # LOAD EXTERNAL FILES - small datasets will usually be contained
        # entirely in a single Excel file, but where formatting or size issues
        # require external files, then names and descriptions are included in
        # the summary information

        external_files = read_block(fields['external'])

        # external file specific validation
        if external_files is not None:
            # simplify to cell values not xlrd.Cells
            strip_ctypes(external_files)

            bad_names = [exf['file'] for exf in external_files if RE_CONTAINS_WSPACE.search(exf['file'])]
            if any(bad_names):
                LOGGER.error('External file names must not contain whitespace: ',
                             extra={'join': bad_names, 'quote': True})

        self.external_files = external_files

        # Load the WORKSHEETS block
        dataworksheet_summaries = read_block(fields['worksheet'])

        # Strip out faulty inclusion of Taxa and Location worksheets in data worksheets
        if dataworksheet_summaries is not None:
            # simplify to cell values not xlrd.Cells
            strip_ctypes(dataworksheet_summaries)

            # - catch inclusion of Taxa and Locations in data worksheets
            cited_sheets = [ws['name'] for ws in dataworksheet_summaries]

            if ('Locations' in cited_sheets) or ('Taxa' in cited_sheets):
                LOGGER.error('Do not include Taxa or Locations metadata sheets in '
                             'Data worksheet details')

                dataworksheet_summaries = [ws for ws in dataworksheet_summaries
                                           if ws['name'] not in ('Locations', 'Taxa')]

                dataworksheet_summaries = None if not dataworksheet_summaries else dataworksheet_summaries

        # Look to see what data is available - must be one or both of data worksheets
        # or external files and validate worksheets if present.
        cited_sheets = set()

        if dataworksheet_summaries is None and external_files is None:
            LOGGER.error("No data worksheets or external files provided - no data.")
        elif dataworksheet_summaries is None:
            LOGGER.info("Only external file descriptions provided")
        elif dataworksheet_summaries is not None:
            # Check sheet names
            cited_sheets = {ws['name'] for ws in dataworksheet_summaries}
            # names not in list of sheets
            bad_names = cited_sheets - self.sheet_names
            if bad_names:
                LOGGER.error('Worksheet names not found in workbook: ',
                             extra={'join': bad_names, 'quote': True})
            # bad external files
            external_in_sheet = {ws['external'] for ws in dataworksheet_summaries
                                 if ws['external'] is not None}
            if self.external_files is not None:
                external_names = {ex['file'] for ex in self.external_files}
            else:
                external_names = set()

            bad_externals = external_in_sheet - external_names
            if bad_externals:
                LOGGER.error('Worksheet descriptions refer to unreported external files.',
                             extra={'join': bad_externals, 'quote': True})

        # Check for existing sheets without description
        extra_names = self.sheet_names - {'Summary', 'Taxa', 'Locations'} - cited_sheets
        if extra_names:
            LOGGER.error('Undocumented sheets found in workbook: ',
                         extra={'join': extra_names, 'quote': True})

        self.dataworksheet_summaries = dataworksheet_summaries

        # summary of processing
        n_errors = CH.counters['ERROR'] - start_errors
        if n_errors > 0:
            LOGGER.info('Summary contains {} errors'.format(n_errors),
                        extra={'indent_before': 1, 'indent_after': 0})
        else:
            LOGGER.info('Summary formatted correctly',
                        extra={'indent_before': 1, 'indent_after': 0})
import requests
# from enforce_typing import enforce_types
import datetime
import re
from safedata_validator.logger import LOGGER, FORMATTER, COUNTER_HANDLER, loggerinfo_push_pop
from safedata_validator.validators import (IsNotSpace, IsString, NoPunctuation, blank_value)
from safedata_validator.extent import Extent

# Compile some regex expressions used in Summary checking
RE_DOI = re.compile(r'https?://(dx.)?doi.org/')
RE_ORCID = re.compile(r'[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]')
RE_EMAIL = re.compile(r'\S+@\S+\.\S+')
RE_NAME = re.compile(r'[^,]+,[ ]?[^,]+')
RE_CONTAINS_WSPACE = re.compile(r'\s')


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
                         ('title', True, None, str),
                         ('description', True, None, str)],
                        True, 'Core fields', True),
                  access=([('access status', True, 'access', str),
                           ('embargo date', False, 'embargo_date', datetime.datetime),
                           ('access conditions', False, 'access_conditions', str)],
                          True, 'Access details', True),
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


    def __init__(self, resources):

        """
        Checks the information in the summary worksheet and looks for the
        metadata and dataset worksheets. The methods are intended to try and
        get as much information as possible from the worksheet: the dictionary
        of metadata returned will have None for any missing data, which should
        be handled by downstream code.


        """

        self.project_id = None
        self.title = None
        self.description = None
        self.access = None
        self.authors = None
        self.permits = None
        self.publication_doi = None
        self.funders = None
        self.keywords = None
        self.temporal_extent = Extent('temporal extent', (datetime.date,),
                                      hard_bounds=resources.extents.temporal_hard_extent,
                                      soft_bounds=resources.extents.temporal_soft_extent)
        self.latitudinal_extent = Extent('latitudinal extent', (float, int),
                                         hard_bounds=resources.extents.latitudinal_hard_extent,
                                         soft_bounds=resources.extents.latitudinal_soft_extent)
        self.longitudinal_extent = Extent('longitudinal extent', (float, int),
                                          hard_bounds=resources.extents.longitudinal_hard_extent,
                                          soft_bounds=resources.extents.longitudinal_soft_extent)
        self.external_files = None
        self.data_worksheets = []

        self._rows = None
        self._ncols = None
        self.n_errors = None
        self.valid_pid = None
        self.validate_doi = False


    @loggerinfo_push_pop('Checking Summary worksheet')
    def load(self, worksheet, sheetnames, validate_doi=False, valid_pid=None):
        """
        Checks the information in a summary worksheet and looks for the
        metadata and dataset worksheets. The methods are intended to try and
        get as much information as possible from the worksheet: the dictionary
        of metadata returned will have None for any missing data, which should
        be handled by downstream code.

        Args:
            worksheet: An openpyxl worksheet instance.
            sheetnames: A set of sheet names found in the workbook.
            validate_doi: Check any publication DOIs, requiring a web connection.
            valid_pid: If provided, an integer or list of integer values that are
                permitted in the Project ID field (usually one for a new dataset
                but more if a published dataset is associated with multiple projects
                and any of those ids would be valid).
        """

        start_errors = COUNTER_HANDLER.counters['ERROR']

        # validate project_id is one of None, an integer or a list of integers
        if valid_pid is None:
            pass
        elif isinstance(valid_pid, int):
            self.valid_pid = [valid_pid]
        elif isinstance(valid_pid, list):
            if not all([isinstance(pid, int) for pid in valid_pid]):
                LOGGER.error("Invalid value in list of project_ids.")
                self.valid_pid = None
            else:
                self.valid_pid = valid_pid
        else:
            LOGGER.error("Provided project id must be an integer or list of integers")
            self.valid_pid = None

        self.validate_doi = validate_doi

        # load worksheet rows, removing blank rows
        # TODO - make 'internal' blank rows an error.
        rows = []
        for this_row in worksheet.iter_rows(values_only=True):
            if not all([blank_value(vl) for vl in this_row]):
                rows.append(this_row)

        self._ncols = worksheet.max_column

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

        # Now process the field blocks
        self._load_core()
        self._load_access_details()
        self._load_authors()
        self._load_keywords()
        self._load_doi()
        self._load_temporal_extent()
        self._load_geographic_extent()
        self._load_funders()
        self._load_permits()
        self._load_external_files()
        self._load_data_worksheets(sheetnames)

        # summary of processing
        self.n_errors = COUNTER_HANDLER.counters['ERROR'] - start_errors
        if self.n_errors > 0:
            LOGGER.info('Summary contains {} errors'.format(self.n_errors))
        else:
            LOGGER.info('Summary formatted correctly')

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
            else:
                LOGGER.info(f'No {title} metadata found')
            return None
        else:
            LOGGER.info(f'Metadata for {title} found: {len(block)} records')

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

    @loggerinfo_push_pop('Loading author metadata')
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

    @loggerinfo_push_pop('Loading keywords metadata')
    def _load_keywords(self):

        keywords = self._read_block(*self.fields['keywords'])

        # extra data validation for keywords
        if keywords:
            keywords = [rec['keywords'] for rec in keywords]
            keywords = NoPunctuation(keywords)
            if not keywords:
                LOGGER.error('Put each keyword in a separate cell, do not separate '
                             'keywords using commas or semi-colons')

            self.keywords = keywords.values

    @loggerinfo_push_pop('Loading permit metadata')
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

    @loggerinfo_push_pop('Loading DOI metadata')
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

    @loggerinfo_push_pop('Loading funding metadata')
    def _load_funders(self):

        # LOOK FOR FUNDING DETAILS - users provide a funding body and a description
        # of the funding type and then optionally a reference number and a URL

        funders = self._read_block(*self.fields['funding'])

        # TODO - currently no check beyond _read_block but maybe actually check
        #        the URL is a URL and maybe even opens? Could use urllib.parse

        self.funders = funders

    @loggerinfo_push_pop('Loading temporal extent metadata')
    def _load_temporal_extent(self):

        temp_extent = self._read_block(*self.fields['date'])

        # temporal extent validation and updating
        if temp_extent is not None:

            start_date = temp_extent[0]['start date']
            end_date = temp_extent[0]['end date']

            if not (isinstance(start_date, datetime.datetime) and
                    isinstance(end_date, datetime.datetime)):
                LOGGER.error('Temporal extents are not date values')
                return

            if not (start_date.time() == datetime.time(0,0) and
                    end_date.time() == datetime.time(0,0)):
                LOGGER.error('Temporal extents should be date not datetime values')
                return

            if start_date > end_date:
                LOGGER.error('Start date is after end date')
                return

            self.temporal_extent.update([start_date.date(), end_date.date()])

    @loggerinfo_push_pop('Loading geographic extent metadata')
    def _load_geographic_extent(self):

        # Geographic extents
        geo_extent = self._read_block(*self.fields['geo'])

        if geo_extent is not None:

            bbox = geo_extent[0]

            if all([isinstance(v, float) for v in bbox.values()]):
                if bbox['south'] > bbox['north']:
                    LOGGER.error('South limit is greater than north limit')
                else:
                    self.latitudinal_extent.update([bbox['south'], bbox['north']])

                if bbox['west'] > bbox['east']:
                    LOGGER.error('West limit is greater than east limit')
                else:
                    self.longitudinal_extent.update([bbox['west'], bbox['east']])

    @loggerinfo_push_pop('Loading external file metadata')
    def _load_external_files(self):

        # LOAD EXTERNAL FILES - small datasets will usually be contained
        # entirely in a single Excel file, but where formatting or size issues
        # require external files, then names and descriptions are included in
        # the summary information

        external_files = self._read_block(*self.fields['external'])

        # external file specific validation - no internal spaces.
        if external_files is not None:

            bad_names = [exf['file'] for exf in external_files
                         if isinstance(exf['file'], str) and RE_CONTAINS_WSPACE.search(exf['file'])]
            if any(bad_names):
                LOGGER.error('External file names must not contain whitespace: ',
                             extra={'join': bad_names})

        self.external_files = external_files

    @loggerinfo_push_pop('Loading data worksheet metadata')
    def _load_data_worksheets(self, sheetnames):

        # Load the WORKSHEETS block
        data_worksheets = self._read_block(*self.fields['worksheet'])

        # Strip out faulty inclusion of Taxa and Location worksheets in
        # data worksheets before considering combinations of WS and external files
        if data_worksheets is not None:

            cited_sheets = [ws['name'] for ws in data_worksheets]

            if ('Locations' in cited_sheets) or ('Taxa' in cited_sheets):
                LOGGER.error('Do not include Taxa or Locations metadata sheets in '
                             'Data worksheet details')

                data_worksheets = [ws for ws in data_worksheets
                                   if ws['name'] not in ('Locations', 'Taxa')] or None

        # Look to see what data is available - must be one or both of data worksheets
        # or external files and validate worksheets if present.
        cited_sheets = set()

        if data_worksheets is None and self.external_files is None:
            LOGGER.error("No data worksheets or external files provided - no data.")
            return
        elif data_worksheets is None:
            LOGGER.info("Only external file descriptions provided")
            return

        # Check sheet names in list of sheets
        cited_sheets = {ws['name'] for ws in data_worksheets}

        # Names not in list of sheets
        for each_ws in data_worksheets:
            # Now match to sheet names
            if each_ws['name'] not in sheetnames:
                if each_ws['external'] is not None:
                    LOGGER.info(f"Worksheet summary {each_ws['name']} recognized as placeholder for "
                                f"external file {each_ws['external']}")
                else:
                    LOGGER.error(f"Data worksheet {each_ws['name']} not found")

        # bad external files
        external_in_sheet = {ws['external'] for ws in data_worksheets
                                if ws['external'] is not None}
        if self.external_files is not None:
            external_names = {ex['file'] for ex in self.external_files}
        else:
            external_names = set()

        bad_externals = external_in_sheet - external_names
        if bad_externals:
            LOGGER.error('Worksheet descriptions refer to unreported external files: ',
                            extra={'join': bad_externals})

        # Check for existing sheets without description
        extra_names = set(sheetnames) - {'Summary', 'Taxa', 'GBIFTaxa', 'NCBITaxa', 'Locations'} - cited_sheets
        if extra_names:
            LOGGER.error('Undocumented sheets found in workbook: ',
                         extra={'join': extra_names})

        self.data_worksheets = data_worksheets

    @loggerinfo_push_pop('Loading access metadata')
    def _load_access_details(self):

        # Load the ACCESS DETAILS block
        access = self._read_block(*self.fields['access'])
        access = access[0]

        # Access specific validation - bad types handled by _read_block
        # - status must be in list of three accepted values
        if isinstance(access['access'], str):

            status = access['access'].lower()
            embargo_date = access['embargo_date']

            if status not in ['open', 'embargo', 'restricted']:
                LOGGER.error(f"Access status must be Open, Embargo or Restricted not {access['access']}")

            if status == 'embargo':

                if embargo_date is None:
                    LOGGER.error('Dataset embargoed but no embargo date provided')
                elif isinstance(embargo_date, datetime.datetime):
                    now = datetime.datetime.now()

                    if embargo_date < now:
                        LOGGER.error('Embargo date is in the past.')
                    elif embargo_date > now + datetime.timedelta(days=2 * 365):
                        LOGGER.error('Embargo date more than two years in the future.')
                    else:
                        LOGGER.info(f'Dataset access: embargoed until {embargo_date }')

                if access['access_conditions'] is not None:
                    LOGGER.error('Access conditions cannot be set on embargoed data.')

            elif status == 'restricted':
                access_conditions = access['access_conditions']

                if embargo_date is not None:
                    LOGGER.error('Do not set an embargo date with restricted datasets')

                if access_conditions is None:
                    LOGGER.error('Dataset restricted but no access conditions specified')
                else:
                    LOGGER.info(f'Dataset access: restricted with conditions {access_conditions}')
            else:
                LOGGER.info(f'Dataset access: {status}')

        self.access = access

    @loggerinfo_push_pop('Loading core metadata')
    def _load_core(self):

        # Now check core rows
        core = self._read_block(*self.fields['core'])
        core = core[0]

        self.title = core['title']
        self.description = core['description']

        # Project ID specific validation
        pid = core['pid']

        # Check the value is in the provided list
        if pid is not None and self.valid_pid is not None and pid not in self.valid_pid:
            LOGGER.error(f'SAFE Project ID in file ({pid}) does not match any '
                         f'provided project ids: ', extra={'join': self.valid_pid})
        else:
            self.project_id = pid

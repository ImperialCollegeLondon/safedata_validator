from itertools import groupby
from logging import ERROR, WARNING, INFO
from typing import Iterable

from openpyxl import worksheet
from openpyxl.utils import get_column_letter

from safedata_validator.validators import (IsNotBlank, IsNotExcelError, IsNotPadded, HasDuplicates, IsNotNA,
                                           IsNumber, IsNotNumericString, IsString, blank_value, valid_r_name)
from safedata_validator.logger import LOGGER, FORMATTER, CH, loggerinfo_push_pop
from safedata_validator.dataset import Dataset
from safedata_validator.locations import Locations
from safedata_validator.taxa import Taxa


MANDATORY_DESCRIPTORS = {'field_type', 'description', 'field_name'}
OPTIONAL_DESCRIPTORS = {'levels', 'method', 'units', 'taxon_name', 'taxon_field',
                        'interaction_name', 'interaction_field', 'file_container'}

class DataWorksheet:

    """
    This class is used to load and check the formatting and content of a data
    worksheet. It requires the dictionary of metadata for the data worksheet
    that was loaded into a Summary object and then access to various Dataset
    wide resources, such as Locations, Taxa and Extents.
    """

    def __init__(self, sheet_meta: dict, worksheet: worksheet, 
                 dataset: Dataset = None) -> None:
        
        """
        A method to load and check the formatting and content of a
        data worksheet, turning the dictionary for this worksheet
        in the dataset object into a DataWorksheet object.
        Args:
            sheet_meta: The metadata dictionary for this worksheet from
                the summary worksheet, containing the worksheet name, 
                title and description and the names of any external files 
                to be associated with this worksheet.
            worksheet: The openpyxl Worksheet instance containing the data
                to be validated.
            dataset: A safedata_validator Dataset object containing 
            taxa: A safedata_validator Taxa instance for use in validating
                taxon fields and traits
            locations: A safedata_validator Locations instance for use in
                validating locations fields.
            summary: A safedata_validator Summary instance for access to
                dataset level metadata like extents.
        """

        # set defaults
        self.name = sheet_meta['name']
        self.description = sheet_meta['description']
        self.title = sheet_meta['title']
        self.max_row = 0
        self.max_col = 0
        self.descriptors = None
        self.taxa_fields = None
        self.field_name_row = None
        self.fields = []
        self.external = None
        self.worksheet = worksheet

        # For sheet meta with an external file, add in the external file name
        if sheet_meta['external'] is not None:
            self.external = sheet_meta['external']

    def load(self, taxa, locations, summary):

        start_errors = CH.counters['ERROR']
        LOGGER.info(f'Checking data worksheet {self.name}')

        self._load_meta_and_size()

        # TODO - openpyxl read-only mode provides row by row data access and used lazy
        # loading to reduce memory usage. Scanning all of a column will need to
        # load all data (which the old implementation using xlrd did, but which
        # we might be also improve to reduce memory overhead by chunking large sheets).
        # Either way, should use row by row ingestion.
        
    def _load_meta_and_size(self):

        # get the  data dimensions
        self.max_col = self.worksheet.max_column
        self.max_row = self.worksheet.max_row

        # trap completely empty worksheets
        if self.max_row == 0:
            LOGGER.error('Worksheet is empty')
            return

        # Read the field metadata first, using a set of rows that should sample well
        # into any data present
        row_sample = len(MANDATORY_DESCRIPTORS) + len(OPTIONAL_DESCRIPTORS) + 10
        ws_rows = self.worksheet.iter_rows(values_only=True, max_row=row_sample)
        field_meta = [r for r in ws_rows]
        
        # Get the sequence of column types
        first_column_types = [type(v[0]) for v in field_meta]
        first_column_rle = [(v, len(list(g))) for v, g in groupby(first_column_types, None)]
        first_column_type_seq = [v[0] for v in first_column_rle]

        # Forbid anything except a narrow set of type sequences: 
        # * descriptor only, 
        # * descriptors + row numbers, 
        # * descriptors + blanks, and
        # * descriptors + row numbers + blanks. 
        #
        # Basically strings, then maybe integers. However, terminal blanks
        # probably indicate unnumbered data rows, so those are allowed through
        # for checking later.
        if first_column_type_seq == [str]:
            LOGGER.error('Cannot parse data: Column A appears to contain '
                         f'more than {row_sample} decriptor names.')
            return
        elif first_column_type_seq not in ([str],
                                         [str, int],
                                         [str, type(None)],
                                         [str, int, type(None)]):
            LOGGER.error('Cannot parse data: Column A must contain a set of '
                         'field descriptors and then row numbers')
            return

        # Reduce to meta, check descriptors and convert to field dictionaries
        n_field_meta = first_column_rle[0][1]
        field_meta = [field_meta[i] for i in range(n_field_meta)]
        field_descriptors = [r[0] for r in field_meta]
        
        # Whitespace padding?
        clean_descriptors = IsNotPadded(field_descriptors)
        if not clean_descriptors:
            # Report whitespace padding and clean up tuples
            LOGGER.error('Whitespace padding in descriptor names: ', 
                         extra={'join': clean_descriptors.failed})
            
            for idx, cln in enumerate(clean_descriptors):
                field_meta[idx] = tuple([cln] + list(field_meta[idx][1:]))
            
            field_descriptors = clean_descriptors.values

        # * Check for mandatory descriptors
        mandatory_missing = set(MANDATORY_DESCRIPTORS).difference(field_descriptors)
        if mandatory_missing:
            LOGGER.error('Mandatory field descriptors missing: ', 
                         extra={'join': mandatory_missing})
            return

        # * Expected descriptors - do _not_ preclude user defined descriptors
        #   but warn about them to alert to typos.
        unknown_descriptors = set(field_descriptors).difference(
                                MANDATORY_DESCRIPTORS.union(OPTIONAL_DESCRIPTORS))
        if unknown_descriptors:
            LOGGER.warning('Unknown field descriptors:', 
                           extra={'join': unknown_descriptors})

        # * Check for field name as _last_ descriptor. This is a bit picky but
        #   it is a historical restriction and makes the dataframe more
        #   readable.
        if field_descriptors[-1] != 'field_name':
                LOGGER.error('Cannot parse data: field_name row is not the last descriptor')
                return
        
        # Repackage field metadata into a list of per field descriptor
        # dictionaries, dropping the first column containing descriptors names
        self.field_name_row = n_field_meta
        field_meta = [dict(zip(field_descriptors, vals)) for vals in zip(*field_meta)]
        del field_meta[0]

        # Insert column index
        for idx, fld in enumerate(field_meta):
            fld['col_idx'] = idx + 1
        
        # Check field names are unique. This doesn't cause as many problems as
        # duplications in Taxa and Locations, which expect certain fields, so
        # warn and continue. Ignore missing mandatory values - these are handled
        # in field checking.
        field_names = [fld['field_name'] for fld in field_meta 
                       if fld['field_name'] is not None]
        dupes = HasDuplicates(field_names)
        if dupes:
            LOGGER.error('Field names duplicated: ', 
                         extra={'join': dupes.duplicated})

        # Lowercase the field types
        for fld in field_meta:
            fld['field_type'] = None if fld['field_type'] is None else fld['field_type'].lower()

        # get taxa field names for cross checking observation and trait data
        self.taxa_fields = [fld['field_name'] for fld in field_meta if fld['field_type'] == 'taxa']

        self.fields = field_meta
    
    def check_data(taxa, locations, summary):


        
        # Now checking the data itself.
        if True:
            # Check if there is anything below the descriptors
            if len(first_column_rle) == 1:
                if is_blank(sheet_meta['external']):
                    LOGGER.error('No data found below field descriptors.')
                else:
                    LOGGER.info('Data table description associated with '
                                'external file {external}'.format(**sheet_meta))

            elif len(first_column_rle) == 2:

                if first_column_rle[1][0] == xlrd.XL_CELL_NUMBER:
                    # If there are numbers, check they are continuous
                    row_numbers = first_column[first_column_rle[0][1]:]
                    ideal = [float(i) for i in range(1, len(row_numbers) + 1)]
                    if row_numbers != ideal:
                        LOGGER.error("Row numbers not consecutive or do not start with 1")
                else:
                    LOGGER.error("Row numbers missing below field descriptors")

            elif len(first_column_rle) == 3:
                n_terminal_blanks = dwsh.max_row - first_column_rle[0][1] - first_column_rle[1][1]
                row_range = list(range(dwsh.max_row - 1, (dwsh.max_row - 1) - n_terminal_blanks, - 1))
                for check_row in row_range:
                    if all(is_blank(val) for val in worksheet.row_values(check_row)):
                        dwsh.max_row -= 1
                    else:
                        LOGGER.error('Un-numbered rows at end of worksheet contain data')
                        break

        # report on detected size
        dwsh.n_data_row = dwsh.max_row - dwsh.field_name_row
        LOGGER.info('Worksheet contains {n_d} descriptors, {o.n_data_row} data rows and '
                    '{n_f} fields'.format(o=dwsh, n_d=len(dwsh.descriptors), n_f=dwsh.max_col - 1),
                    extra={'indent_after': 2})

 
        # check the data in each field against the metadata
        for meta in metadata:
            # read the values and check them against the metadata
            data = worksheet.col(meta['col_idx'], dwsh.field_name_row, dwsh.max_row)
            self.check_field(dwsh, meta, data)

        # add the new DataWorksheet into the Dataset
        self.dataworksheets.append(dwsh)

        # reporting
        n_errors = CH.counters['ERROR'] - start_errors
        if n_errors > 0:
            LOGGER.info('Dataframe contains {} errors'.format(n_errors),
                        extra={'indent_before': 1})
        else:
            LOGGER.info('Dataframe formatted correctly',
                        extra={'indent_before': 1})


class BaseField:

    # Class instances describe the field_types handled by the class
    # and sets the required descriptors for those field_types. The 
    # no_validation variable is used to suppress data validation for
    # use on e.g. comments fields.

    field_types = None
    required_descriptors = MANDATORY_DESCRIPTORS
    no_validation = False
    check_taxon_meta = False
    check_interaction_meta = False

    def __init__(self, meta: dict, 
                 dwsh: DataWorksheet = None,
                 dataset: Dataset = None,
                 taxa: Taxa = None,
                 locations: Locations = None) -> None:
        
        """
        Base class to check field metadata and then test whether the data
        contents of a field are compatible with the metadata for the field.
        Largely used for checking data worksheets, but also by other methods.
        This base class defines the core checking functions as a base class
        and then subclasses implement data type specific checking.

        In order to work with the row by row data loading from openpyxl,
        the BaseField class does _not_ emit logging messages as they occur 
        but holds a stack of messages created during initiation and as 
        chunks of data are ingested. This logging can then be emitted using 
        the report() method, once all the data has been ingested.

        The base class also has a much more complex signature than the base
        functionality requires. Various subclasses need access to:

            * Dataset level information - extents, taxa and locations.
            * Dataworksheet level information - taxon fields

        Rather than muck around with changing subclass signatures and 
        kwargs, the Base class makes all information available to all
        subclasses.

        Args: 
            meta: A dictionary of field metadata
            dwsh: A DataWorksheet, used to pass worksheet level information,
                which at the moment is just taxa fields.
        """

        self.meta = meta
        self.dwsh = dwsh
        self.dataset = dataset
        self.taxa = taxa
        self.locations = locations

        self.log_stack = []
        self.n_rows = 0
        self.n_na = 0
        self.n_blank = 0
        self.n_excel_error = 0
        self.bad_values = []
        self.bad_rows = [] # TODO check on implementation of this

        # field_name in meta is guaranteed by Dataset loading so 
        # use this to get a reporting name or sub in a column letter
        # _check_meta then reports on bad values.
        if blank_value(self.meta['field_name']):
            field_name = get_column_letter(self.meta['col_idx'])
        else:
            field_name = str(self.meta['field_name'])

        self.log_stack.append((INFO, f'Checking Column {field_name}'))

        # Now check required descriptors - values are present, non-blank and 
        # are unpadded strings (currently all descriptors are strings).
        for dsc in self.required_descriptors:
            self._check_meta(dsc)
        
        # Specific check that field name is valid - the column letter codes are always
        # valid, so missing names won't trigger this.
        if not valid_r_name(field_name):
            self.log_stack.append((ERROR, f"Field name is not valid: {repr(field_name)}. "
                                          "Common errors are spaces and non-alphanumeric "
                                          "characters other than underscore and full stop"))

        # TODO - reimplement these two class attributes and steps as a single
        #        subclass specific extra meta method stub (_check_meta_extra?)
        #        that can be overloaded as required by subclasses?
        if self.check_taxon_meta:
            self._check_taxon_meta()
        
        if self.check_interaction_meta:
            self._check_interaction_meta()

    def _check_meta(self, descriptor):
        """
        A standardised checker to see if a required descriptor is present for
        a field and that it isn't simply empty or whitespace. Messages are 
        added to the log_stack if required.

        Args:
            meta: A dictionary of field metadata descriptors
            descriptor: The name of the descriptor to check.
        Returns:
            A boolean, with True showing no problems and False showing
            that warnings occurred.
        """

        if descriptor not in self.meta:
            self.log_stack.append((ERROR, f'{descriptor} descriptor missing'))
            return False

        val = self.meta[descriptor]
        
        if blank_value(val):
            self.log_stack.append((ERROR, f'{descriptor} descriptor is blank'))
            self.meta[descriptor] = None  # standardise whitestring to None
            return False
        elif not isinstance(val, str):
            self.log_stack.append((ERROR, f'{descriptor} descriptor is not a string: {repr(val)}'))
            return False
        elif val != val.strip():
            self.log_stack.append((ERROR, f'{descriptor} descriptor has whitespace padding: {repr(val)}'))
            self.meta[descriptor] = val.strip()
            return False
        else:
            return True

    def _check_taxon_meta(self, taxa_fields):

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

        tx_nm_prov = ('taxon_name' in self.meta) and (not is_blank(self.field_meta['taxon_name']))
        tx_fd_prov = ('taxon_field' in self.field_meta) and (not is_blank(self.field_meta['taxon_field']))

        if tx_nm_prov and tx_fd_prov:
            LOGGER.error('Taxon name and taxon field both provided, use one only')
            return False
        elif tx_nm_prov and meta['taxon_name'] not in self.taxon_names:
            LOGGER.error('Taxon name not found in the Taxa worksheet')
            return False
        elif tx_fd_prov and meta['taxon_field'] not in taxa_fields:
            LOGGER.error("Taxon field not found in this worksheet")
            return False
        elif not tx_nm_prov and not tx_fd_prov:
            LOGGER.error("One of taxon name or taxon field must be provided")
            return False
        else:
            if tx_nm_prov:
                self.taxon_names_used.update([meta['taxon_name']])
            return True

    def _check_interaction_meta(self, taxa_fields):

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
            LOGGER.error("At least one of interaction name or interaction field must be provided")
            return False
        else:
            if iact_nm_prov:
                # get the taxon names and descriptions from interaction name providers
                int_nm_lab, int_nm_desc = self._parse_levels(meta['interaction_name'])
                # add names to used taxa
                self.taxon_names_used.update(int_nm_lab)
                # check they are found
                if not all([lab in self.taxon_names for lab in int_nm_lab]):
                    LOGGER.error('Unknown taxa in interaction_name descriptor')
                    nm_check = False
                else:

                    nm_check = True
            else:
                int_nm_lab, int_nm_desc = [(), ()]
                nm_check = True

            if iact_fd_prov:
                # check any field labels match to known taxon fields
                int_fd_lab, int_fd_desc = self._parse_levels(meta['interaction_field'])
                if not all([lab in taxa_fields for lab in int_fd_lab]):
                    LOGGER.error('Unknown taxon fields in interaction_field descriptor')
                    fd_check = False
                else:
                    fd_check = True
            else:
                int_fd_lab, int_fd_desc = [(), ()]
                fd_check = True

            if len(int_nm_lab + int_fd_lab) < 2:
                LOGGER.error('At least two interacting taxon labels or fields must be identified')
                num_check = False
            else:
                num_check = True

            if nm_check and fd_check and num_check:
                return True
            else:
                return False

    def _parse_levels(self, txt):
        """
        Splits up category information formatted as label:desc;label:desc, which
        is used in both levels for categorical data and interaction descriptors.
        Args:
            txt: The text string to parse
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
            self.log_stack.append((ERROR, 'Extra colons in level description.'))
            

        # standardise descriptions
        if all([pt == 1 for pt in n_parts]):
            parts = [[pt[0], None] for pt in parts]
        elif all([pt > 1 for pt in n_parts]):
            # truncate extra colons
            parts = [pt[0:2] for pt in parts]
        else:
            self.log_stack.append((ERROR, 'Provide descriptions for either all or none of the categories'))
            parts = [pt[0:2] if len(pt) >= 2 else [pt[0], None] for pt in parts]

        level_labels, level_desc = zip(*parts)

        # - repeated labels?
        if len(set(level_labels)) < len(level_labels):
            self.log_stack.append((ERROR, 'Repeated level labels'))

        # - check for numeric level names: integers would be more common
        #   but don't let floats sneak through either!
        level_labels = IsNotNumericString(level_labels)
        if not level_labels:
            self.log_stack.append((ERROR, 'Numeric level names not permitted'))

        # Remove white space around the labels: simple spacing in the text
        # makes it easier to read and insisting on no space is unnecessary
        level_labels = [vl.strip() for vl in level_labels]

        return level_labels, level_desc

    @staticmethod
    def validate_geo_extent(extent, which):

        if which == 'latitude':
            bnds = [-90, -4, 8, 90]
        elif which == 'longitude':
            bnds = [-180, 108, 120, 180]

        if extent[0] > extent[1]:
            LOGGER.error('{0}: lower bound greater than upper bound'.format(which.capitalize()))
            return False

        out_of_bounds = extent[0] < bnds[0] or extent[1] > bnds[3]
        out_of_borneo = bnds[0] <= extent[0] < bnds[1] or bnds[2] < extent[1] <= bnds[3]

        if out_of_bounds:
            LOGGER.error('{0} values not in valid range[{1[0]}, {1[3]}]: '
                         '{2}'.format(which.capitalize(), bnds, extent))
            return False
        elif out_of_borneo:
            LOGGER.warn('{0} values not in Borneo [{1[1]}, {1[2]}]: '
                        '{2}'.format(which.capitalize(), bnds, extent))
            return True
        else:
            return True






        # run consistency checks where needed and trap unknown field types
        if field_type == 'date':
            self.check_field_datetime(meta, data, which='date')
        elif field_type == 'datetime':
            self.check_field_datetime(meta, data, which='datetime')
        elif field_type == 'time':
            self.check_field_datetime(meta, data, which='time')
        elif field_type == 'taxa':
            self.check_field_taxa(data)
        elif field_type == 'location':
            self.check_field_locations(data)
        elif field_type in ['categorical', 'ordered categorical']:
            self.check_field_categorical(meta, data)
        elif field_type == 'numeric':
            self.check_field_numeric(meta, data)
        elif field_type == 'abundance':
            self.check_field_abundance(meta, data, dwsh.taxa_fields)
        elif field_type in ['categorical trait', 'ordered categorical trait']:
            self.check_field_trait(meta, data, dwsh.taxa_fields, which='categorical')
        elif field_type == 'numeric trait':
            self.check_field_trait(meta, data, dwsh.taxa_fields, which='numeric')
        elif field_type in ['categorical interaction', 'ordered categorical interaction']:
            self.check_field_interaction(meta, data, dwsh.taxa_fields, which='categorical')
        elif field_type == 'numeric interaction':
            self.check_field_interaction(meta, data, dwsh.taxa_fields, which='numeric')
        elif field_type == 'latitude':
            # check field geo expects values in data not xlrd.Cell
            data = [dt.value for dt in data]
            self.check_field_geo(meta, data, which='latitude')
        elif field_type == 'longitude':
            # check field geo expects values in data not xlrd.Cell
            data = [dt.value for dt in data]
            self.check_field_geo(meta, data, which='longitude')
        elif field_type == 'file':
            data = [dt.value for dt in data]
            self.check_field_file(meta, data)
        elif field_type in ['replicate', 'id']:
            # We've looked for missing data, no other constraints.
            pass
        elif field_type == 'comments':
            pass
        elif field_type is None:
            LOGGER.error('Field type is empty')
        else:
            LOGGER.error('Unknown field type {field_type}'.format(**meta))

    def validate_data(self, data: list, **kwargs):
        """Validates a list of data provided for a field. 

        This method runs the common shared validation and returns the cleaned
        values, excluding blanks, NA and any Excel errors.

        The method is overloaded by subclasses to provide field specific
        testing. In this case using the following ensures that only the data
        that passes the common checks is subjected to extra testing:

            data = super().validate_data(data, **kwargs)
        """

        # TODO - handle blank columns/metadata at Dataworksheet level? -etc.
        # # Skip any field with no user provided metadata or data
        # blank_data = [is_blank(cl.value) for cl in data]
        # blank_meta = [is_blank(vl) for ky, vl in list(meta.items())
        #               if ky not in ['col_idx', 'column']]
        # name_is_string = isinstance(meta['field_name'], str)

        # Unless any input is ok (comments fields, basically) filter out 
        # missing and blank data. Also check for Excel formula errors.
        self.n_rows += len(data)

        # Look for NAs
        data = IsNotNA(data, keep_failed=False)
        if not data:
            self.n_na += len(data.failed)

        # Look for blank data
        data = IsNotBlank(data, keep_failed=False)
        if not data:
            self.n_blank += len(data.failed)

        # Look for formula errors: 

        # With openpyxl, using data_only=False loads function cells with
        # data_type = 'f', but using data_only=True returns cells with a 
        # data_type for the resulting value _or_ an 'e' code for a function
        # error. Using value_only=True in iter_row for speed loses that cue,
        # so the values simply come back as a strings that have to be matched
        # to the of Excel error codes. See also:
        #
        # https://groups.google.com/g/openpyxl-users/c/iNi1MKSP-Bc/m/2q_7cHavBQAJ

        data = IsNotExcelError(data, keep_failed=False)
        if not data:
            self.n_excel_error += len(data.failed)

        # Return cleaned data for use in further checking.
        return data.values

    def report(self):

        # Emit the accumulated messages
        for log_level, msg in self.log_stack:
            LOGGER.log(log_level, msg)
        
        # Run any final logging 
        # - warn about NAs
        if self.n_na:
            LOGGER.warning(f'{self.n_na} / {self.n_rows} values missing')

        # We don't tolerate blank cells
        if self.n_blank:
            LOGGER.error(f'{self.n_blank} cells are blank or contain only whitespace text')

        # Report on excel error codes
        if self.n_excel_error:
            LOGGER.error(f'{self.n_excel_error} cells contain Excel formula errors')


class NumericField(BaseField):
    """
    Subclass of BaseField to check for numeric data
    """
    field_types = ('numeric')
    required_descriptors = MANDATORY_DESCRIPTORS.union(['method', 'units'])

    def validate_data(self, data: list, **kwargs):
        
        data = super().validate_data(data, **kwargs)
        
        numeric = IsNumber(data)

        if not numeric:
            self.log_stack.append((ERROR, 'Cells contain non-numeric values'))


class CategoricalField(BaseField):
    """
    Subclass of BaseField to check for categorical data
    """

    field_types = ('categorical', 'ordered_categorical')
    required_descriptors = MANDATORY_DESCRIPTORS.union(['levels'])

    def __init__(self, meta: dict, dwsh: DataWorksheet = None, dataset: Dataset = None, taxa: Taxa = None, locations: Locations = None) -> None:

        super().__init__(meta, dwsh=dwsh, dataset=dataset, taxa=taxa, locations=locations)

        # Additional code to validate the levels metadata and store a set of values
        self.level_labels = set()
        self.reported_levels = set()
        levels = meta.get('levels')
        
        if levels is not None and isinstance(levels, str) and not levels.isspace():
            level_labels, level_desc = self._parse_levels(levels)
            self.level_labels = set(level_labels)

    def validate_data(self, data: list, **kwargs):

        data = super().validate_data(data, **kwargs)

        # Now look for consistency: get the unique values reported in the
        # data, convert to unicode to handle checking of numeric labels and
        # then check the reported levels are a subset of the descriptors.
        # XLRD reads all numbers as floats, so coerce floats back to int
        data = IsString(data, keep_failed=False)
        if not data:
            self.log_stack.append((ERROR, 'Cells contain non-text values'))

        self.reported_levels.update(data)

    def report(self):

        super().report()

        extra = self.reported_levels.difference(self.level_labels)
        unused = self.level_labels.difference(self.reported_levels)

        if extra:
            LOGGER.error('Categories found in data missing from levels descriptor: ',
                            extra={'join': extra})
        if unused:
            LOGGER.error('Categories found in levels descriptor not used in data: ',
                            extra={'join': unused})


class DatetimeField(BaseField):

    field_types = ('date', 'time', 'datetime')

    def ingest(self, data, which='datetime'):

        """
        Checks for data consistency in date and datetime fields. Excel
        doesn't distinguish, both are loaded as datetime.datetime objects
        so we check that the values are compatible with the user provided
        field type. Also handle datetimes provided as strings in common
        ISO formats
        Args:
            meta: The field metadata, to be updated with the range
            data: A list of data values, allegedly containing datetimes or dates
            which: The datetime type to be checked for.
        """

        # Check types and allow all xlrd.XL_CELL_DATE or all xlrd.XL_CELL_STRING
        # but not a mix
        cell_types = {dt.ctype for dt in data}
        extent = None

        if not data:
            # Data is empty - external file description
            pass
        elif cell_types == {xlrd.XL_CELL_DATE}:

            data = [dt.value for dt in data]

            # For date and time options, check decimal components
            if which == 'date':
                contains_time = [dt % 1 > 0 for dt in data]
                if any(contains_time):
                    LOGGER.error('Some values also contain time components')

            elif which == 'time':
                contains_date = [1 <= dt for dt in data]
                if any(contains_date):
                    LOGGER.error('Some values also contain date components')

            # For date containing options, check for actual dates and get extents
            if which in ['date', 'datetime']:
                if any(0.0 <= dt < 1.0 for dt in data):
                    LOGGER.error('Some values _only_  contain time components')

            # get extents
            if len(data):
                extent = (xlrd.xldate_as_datetime(min(data), self.workbook.datemode),
                            xlrd.xldate_as_datetime(max(data), self.workbook.datemode))

                if which == 'time':
                    extent = tuple(vl.time() for vl in extent)

        elif cell_types == {xlrd.XL_CELL_TEXT}:

            data = [dt.value for dt in data]

            # internal function to check dates
            def parse_datetime(dt, which='date'):
                try:
                    if which == 'time':
                        dt = parser.isoparser().parse_isotime(dt)
                    else:
                        dt = parser.isoparse(dt)

                    return dt, True
                except ValueError:
                    return dt, False

            parsed_data = [parse_datetime(dt, which) for dt in data]
            good_data = set((dt[0] for dt in parsed_data if dt[1]))
            bad_data = set((dt[0] for dt in parsed_data if not dt[1]))

            if len(good_data):
                extent = (min(good_data), max(good_data))
                # intentionally scrub out timezones - Excel dates don't have them and can't
                # compare datetimes with mixed tzinfo status
                extent = tuple(dt.replace(tzinfo=None) for dt in extent)

            if len(bad_data):
                LOGGER.error('Problem in parsing date/time strings: ', extra={'join': bad_data})

        elif cell_types == {xlrd.XL_CELL_DATE, xlrd.XL_CELL_TEXT}:
            first_cell_type = data[0].ctype
            first_text = next(idx for idx, val in enumerate(data) if val.ctype != first_cell_type)
            LOGGER.error('Field contains data with mixed test and date formatting. Use either Excel date formats or '
                            'ISO formatted text. Note that text can look _exactly_ like an Excel date or '
                            'time cell, you may need to copy the column and format as numbers to spot '
                            'errors. The number of the first row with formatting not matching the first '
                            'value is: ' + str(first_text))

        else:
            LOGGER.error('Field contains cells formatted as neither text nor date')

        # range and extent setting
        if extent is not None:
            meta['range'] = extent
            if which in ['date', 'datetime']:
                self.update_extent(extent, datetime.datetime, 'temporal_extent')


class TaxaField(BaseField):
    """
    Checks if all the values provided in a Taxon field are found
    in the Taxa worksheet.
    """

    field_types = 'taxa'
    required_descriptors = MANDATORY_DESCRIPTORS

    def __init__(self, meta: dict, dwsh: DataWorksheet = None, dataset: Dataset = None, taxa: Taxa = None, locations: Locations = None) -> None:
        super().__init__(meta, dwsh=dwsh, dataset=dataset, taxa=taxa, locations=locations)

        if self.taxa is None:
            self.log_stack.append((ERROR, 'No taxon details provided for dataset'))
        
        self.taxa_found = set()
    
    def validate_data(self, data: list, **kwargs):

        data = super().validate_data(data, **kwargs)

        data = IsString(data, keep_failed=False)

        if not data:
            self.log_stack.append((ERROR, 'Cells contain non-string values'))

        self.taxa_found.update(data)
    
    def report(self):

        super().report()

        # TODO - not sure about this - no other fields test for emptiness?
        if self.taxa_found == set():
            LOGGER.error('No taxa loaded')
            return
        
        if self.taxa is not None:
            extra_taxa = self.taxa_found.difference(self.taxa.taxon_names)

            if extra_taxa:
                LOGGER.error('Includes unreported taxa: ',
                                extra={'join': extra_taxa})

            # add the found taxa to the list of taxa used
            self.taxa.taxon_names_used.update(self.taxa_found)

def check_field_locations(self, data):

    """
    Checks if all the values provided in a Locations field are
    found in the Locations worksheet, reporting to the Messages instance.
    Args:
        data: A list of data values, allegedly taxon names
    """

    # location names should be strings but we allow integer point numbers
    data = [str(int(dt.value)) if dt.ctype == xlrd.XL_CELL_NUMBER
            else str(dt.value) for dt in data]

    # check if locations are all provided
    found = set(data)
    if self.locations == set():
        LOGGER.error('No locations loaded')
    elif not found.issubset(self.locations):
        LOGGER.error('Includes locations missing from Locations worksheet:',
                        extra={'join': found - self.locations, 'quote': True})

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
    nums = [dt.value for dt in data if dt.ctype == xlrd.XL_CELL_NUMBER]

    if len(nums) < len(data):
        LOGGER.error('Field contains non-numeric data')

    # update the field metadata if there is any data.
    if len(nums):
        meta['range'] = (min(nums), max(nums))


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
    # data is all numeric, as it claims to be. Keep dates separate because
    # they are a pain to find if included as their numeric storage value
    good, bad, ugly = [], [], []
    for dt in data:
        if dt.ctype == xlrd.XL_CELL_NUMBER:
            good.append(dt.value)
        elif dt.ctype == xlrd.XL_CELL_DATE:
            ugly.append(dt.value)
        else:
            bad.append(dt.value)

    if bad:
        LOGGER.error('Field contains non-numeric data: ',
                        extra={'join': set(bad), 'quote': True})

    if ugly:
        ugly = [xlrd.xldate_as_datetime(ug, self.workbook.datemode).isoformat() for ug in ugly]
        LOGGER.error('Field contains date/time formatted data - format will differ from '
                        'from the values shown: ', extra={'join': ugly})

    # update the field metadata if there is any data
    if len(good):
        meta['range'] = (min(good), max(good))

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
    nums = [dt for dt in data if isinstance(dt, float)]

    if len(nums) < len(data):
        LOGGER.error('Field contains non-numeric data')
        if any([RE_DMS.search(str(dt)) for dt in data]):
            LOGGER.warn('Possible degrees minutes and seconds formatting? Use decimal degrees')

    # Check the locations
    if len(nums):
        min_geo = float(min(nums))
        max_geo = float(max(nums))
        extent = (min_geo, max_geo)
        valid = self.validate_geo_extent(extent, which)

        if valid:
            # update the field metadata and the dataset extent
            meta['range'] = extent
            # Look up the extent name to update and then update it
            which_extent = {'latitude': 'latitudinal_extent',
                            'longitude': 'longitudinal_extent'}
            self.update_extent(extent, float, which_extent[which])

def check_field_file(self, meta, data):
    """
    Checks file fields. The data values need to match to an external file
    or can be contained within an archive file provided in the 'file_container'
    descriptor.
    Args:
        meta: A dictionary of metadata descriptors for the field
        data: A list of data values
    """
    
    if self.external_files is None:
        LOGGER.error('No external files listed in Summary')
        return
    
    external_names = {ex['file'] for ex in self.external_files}

    if 'file_container' in meta and meta['file_container'] is not None:
        if meta['file_container'] not in external_names:
            LOGGER.error(
                "Field file_container value not found in external files: {}".format(meta['file_container']))
    else:
        missing_files = set(data) - external_names
        if missing_files:
            LOGGER.error("Field contains files not listed in external files: ",
                            extra={'join': missing_files})
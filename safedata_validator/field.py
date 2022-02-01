import datetime
from itertools import groupby, islice
from collections import OrderedDict
from logging import CRITICAL, ERROR, WARNING, INFO
from openpyxl import worksheet
from openpyxl.utils import get_column_letter
from dateutil import parser
from typing import List
from safedata_validator import dataset

from safedata_validator.validators import (IsInSet, IsNotBlank, IsNotExcelError, IsNotPadded,
                                           HasDuplicates, IsNotNA, IsNumber, 
                                           IsNotNumericString, IsString, IsLocName,
                                           blank_value, valid_r_name, RE_DMS)

from safedata_validator.logger import LOGGER, FORMATTER, CH, loggerinfo_push_pop
from safedata_validator.dataset import Dataset
from safedata_validator.locations import Locations
from safedata_validator.taxa import Taxa

# These are lists, not sets because lists preserve order for preserving logging
# message order in unit testing.
MANDATORY_DESCRIPTORS = ['field_type', 'description', 'field_name']
OPTIONAL_DESCRIPTORS = ['levels', 'method', 'units', 'taxon_name', 'taxon_field',
                        'interaction_name', 'interaction_field', 'file_container']


class DataWorksheet:

    """
    This class is used to load and check the formatting and content of a data
    worksheet. It requires the dictionary of metadata for the data worksheet
    that was loaded into a Summary object and then access to various Dataset
    wide resources, such as Locations, Taxa and Extents.
    """

    @loggerinfo_push_pop('Checking data worksheet')
    def __init__(self, 
                 sheet_meta: dict,
                 dataset: Dataset = None,
                 ) -> None:
        
        """
        This class is used to load and validate the formatting and content of a
        data worksheet. Checking a DataWorksheet has the following workflow:

        1) Create a DataWorksheet instance using the field metadata.
        2) Use calls to the validate_data_rows method to pass in field
           contents - this can be done repeatedly to handle chunked inputs.
        3) Call the report method to obtain all of the logging associated 
           with data validation and any final checks.

        Args:
            sheet_meta: The metadata dictionary for this worksheet from
                the summary worksheet, containing the worksheet name, 
                title and description and the names of any external files 
                to be associated with this worksheet.
            field_meta: The field metadata as an OrderedDict - ordered from
                top to bottom as they appear in the source file - keyed by
                field descriptor with values as a list of field values for
                that descriptor.
            dataset: A safedata_validator Dataset object, providing Taxa, 
                Locations and Summary information.
        """

        # Set initial values
        
        # TODO - checks on sheetmeta 
        self.name = sheet_meta['name']
        self.description = sheet_meta['description']
        self.title = sheet_meta['title']
        # For sheet meta with an external file, add in the external file name
        self.external = sheet_meta.get('external')

        # Links to dataset for dataset level information, and taxa and
        # location validation lists.
        self.dataset = dataset

        # Create the field meta attributes, populated using validate_field_meta
        self.fields_loaded = False
        self.field_meta = None
        self.field_types = None
        self.n_fields = None
        self.n_descriptors = None
        self.fields = None
        self.taxa_fields = None

        # Keep track of row numbering
        self.n_row = 0
        self.current_row = 1
        self.row_numbers_sequential = True
        self.row_numbers_missing = False
        self.row_numbers_noninteger = False
        self.blank_rows = False
        self.trailing_blank_rows = False
        self.start_errors = CH.counters['ERROR']
        self.n_errors = 0

    @loggerinfo_push_pop('Validating field metadata')
    def validate_field_meta(self, field_meta: OrderedDict) -> List[dict]:
        """This method is used to add and validate field_meta. The requirement
        of an OrderedDict for field_meta is to check that the field name is the
        last row in the metadata. A dict should by enough in Python 3.7+ where
        insertion order is guaranteed, but OrderedDict is currently used to
        cover 3.6 and to make the requirement explicit.

        The field meta can include trailing empty fields - these will be checked
        when data rows are added.
        """

        # Checking field_meta - TODO more checking of structure

        # - Descriptors
        # Clean off whitespace padding?
        clean_descriptors = IsNotPadded(field_meta.keys())
        if not clean_descriptors:
            # Report whitespace padding and clean up tuples
            LOGGER.error('Whitespace padding in descriptor names: ', 
                         extra={'join': clean_descriptors.failed})
            
            # Order preserved in dict and validator
            cleaned_entries = [(ky, val) for ky, val 
                               in zip(clean_descriptors, field_meta.values())]
            field_meta = OrderedDict(cleaned_entries)

        # * Expected descriptors - do _not_ preclude user defined descriptors
        #   but warn about them to alert to typos. Missing descriptors vary
        #   with field type so are handled in BaseField.__init__()
        unknown_descriptors = set(field_meta.keys()).difference(
                                set(MANDATORY_DESCRIPTORS).union(OPTIONAL_DESCRIPTORS))
        if unknown_descriptors:
            LOGGER.warning('Unknown field descriptors:', 
                           extra={'join': unknown_descriptors})

        # * Check for field name as _last_ descriptor. This is primarily to
        #   make it easy to read the dataframe - just skip the other descriptors
        if list(field_meta.keys())[-1] != 'field_name':
                LOGGER.error('Field_name row is not the last descriptor')
        
        # Repackage field metadata into a list of per field descriptor
        # dictionaries, _importantly_ preserving the column order.
        field_meta = [dict(zip(field_meta.keys(), val)) 
                        for val in zip(*field_meta.values())]

        # Check provided field names are unique. This doesn't cause as many
        # problems as duplications in Taxa and Locations, which expect certain
        # fields, so warn and continue. Ignore empty mandatory values - these
        # are handled in field checking.
        field_names = [fld['field_name'] for fld in field_meta 
                       if fld['field_name'] is not None]
        dupes = HasDuplicates(field_names)
        if dupes:
            LOGGER.error('Field names duplicated: ', 
                         extra={'join': dupes.duplicated})

        # Lowercase the field types
        for fld in field_meta:
            fld['field_type'] = None if fld['field_type'] is None else fld['field_type'].lower()

        # Populate the instance variables
        self.fields_loaded = True
        self.field_meta = field_meta
        self.n_fields = len(field_meta)
        self.n_descriptors = len(field_meta[0].keys())

        # get taxa field names for cross checking observation and trait data
        self.taxa_fields = [fld['field_name'] for fld in self.field_meta 
                            if fld['field_type'] == 'taxa']

        # Now initialise the Field objects using the mapping of field types
        # to the BaseField subclasses. This needs to handle trailing fields
        # containing _no_ metadata (e.g. from inaccurate Excel sheet bounds or
        # CSV comments) but also fields with missing field_type (but other 
        # metadata).
        #
        # - trailing empty fields are assumed to be blank columns
        #   and validate_data_rows should test for content.
        # - otherwise, these are assumed to be fields and should
        #   trigger warnings about descriptors

        # Get logical flags for empty and trailing empty metadata 
        field_meta_empty = [set(vl.values()) == set([None]) for vl in self.field_meta]
        trailing_empty  = [all(field_meta_empty[- (n + 1):]) 
                           for n in range(0, len(field_meta_empty))]
        trailing_empty.reverse()

        # Get the type map and field list
        field_subclass_map = BaseField.field_type_map()
        unknown_field_types = set()
        self.fields = []

        for col_idx, (tr_empty, fd_empty, fmeta) in \
            enumerate(zip(trailing_empty, field_meta_empty, self.field_meta)):
            
            fmeta['col_idx'] = col_idx + 1

            # Consider cases
            if tr_empty:
                # Hopefully empty trailing field
                self.fields.append(EmptyField(fmeta))
            elif fd_empty:
                # Empty field within other fields
                self.fields.append(BaseField(fmeta, dwsh=self, dataset=self.dataset))
            elif fmeta['field_type'] is None:
                # Non-empty field of None type - use BaseField for basic checks
                self.fields.append(BaseField(fmeta, dwsh=self, dataset=self.dataset))
            elif fmeta['field_type'] not in field_subclass_map:
                #  Non-empty field of unknown type - use BaseField for basic checks
                unknown_field_types.add(fmeta['field_type'])
                self.fields.append(BaseField(fmeta, dwsh=self, dataset=self.dataset))
            else:
                # Known field type.
                fld_class = field_subclass_map[fmeta['field_type']]
                self.fields.append(fld_class(fmeta, dwsh=self, dataset=self.dataset))

        if unknown_field_types:
            LOGGER.error('Unknown field types: ', extra={'join': unknown_field_types})
        
    def validate_data_rows(self, data_rows):
        """Method to pass rows of data into the field checkers for a
        dataworksheet. The data is expected to be as an list of rows
        from a data file, with the first entry in each row being a row
        number.
        """

        if not self.fields_loaded:
            LOGGER.critical('No fields defined - use validate_field_meta')
            return

        # Check the lengths of the rows - this should log critical
        # because it is likely to be a programming error, not a user error.
        if len(data_rows) == 0:
            LOGGER.critical('Empty data_rows passed to validate_data_rows')
            return

        row_lengths = set([len(rw) for rw in data_rows])

        if len(row_lengths) != 1:
            LOGGER.critical('Data rows of unequal length - cannot validate')
            return
        elif row_lengths.pop() != (self.n_fields + 1):
            LOGGER.critical('Data rows not of same length as field metadata - cannot validate')
            return
        
        # Handle empty rows.
        blank_set = set([None])
        blank_row = [set(vals) == blank_set for vals in data_rows]

        trailing_blank_row  = [all(blank_row[- (n + 1):]) 
                                for n in range(0, len(blank_row))]
        trailing_blank_row.reverse()

        internal_blank_row = [blnk & (trail == False) 
                              for blnk, trail in zip(blank_row, trailing_blank_row)]

        # Internals within a set of row are automatically bad, but trailing blank
        # rows are only bad if more data is loaded below them.
        if any(internal_blank_row):
            self.blank_rows = True
        
        if not self.trailing_blank_rows and any(trailing_blank_row):
            self.trailing_blank_rows = True
        elif self.trailing_blank_rows and not all(trailing_blank_row):
            self.blank_rows = True

        # Drop blank rows for further processing
        data_rows = [drw for drw, is_blank in zip(data_rows, blank_row) if not is_blank]

        # Count the rows, then convert the values into columns and extract the row numbers
        self.n_row += len(data_rows)
        data_cols = list(zip(*data_rows))
        row_numbers = data_cols.pop(0)

        # Check for bad values (blanks, non integers) in row numbers
        row_numbers = IsNotBlank(row_numbers, keep_failed=False)

        if not row_numbers:
            self.row_numbers_missing = True
        
        row_numbers = row_numbers.values
        
        if any([not isinstance(vl, int) for vl in row_numbers]):
            self.row_numbers_noninteger = True
        
        # Check the row numbering - skip if this has already failed or
        # the row numbers have already included None or non-integers
        if self.row_numbers_sequential and not (self.row_numbers_noninteger or self.row_numbers_missing):
            expected_numbers = list(range(self.current_row, 
                                           self.current_row + len(row_numbers)))
            self.current_row += len(row_numbers)

            if row_numbers != expected_numbers:
                self.row_numbers_sequential = False
        else:
            self.current_row += len(row_numbers)

        # Now feed the sets of values into the Field validation
        for data, field_inst in zip(data_cols, self.fields):
            field_inst.validate_data(data)

    @loggerinfo_push_pop('Validating field data')
    def report(self):
        """ The validate_data_rows method accumulates logging messages
        individual fields. This method causes these field logs to be 
        emitted and then carries out final checks and reporting across
        all of the field data.
        """
        
        if not self.n_row:
            if self.external is None:
                LOGGER.error('No data passed for validation.')
            else:
                LOGGER.info('Data table description associated with '
                            f'external file {self.external}')
        else:
            for fld in self.fields:
                fld.report()

        # Report on row numbering
        if self.row_numbers_missing:
            LOGGER.error("Missing row numbers in data")
        
        if self.row_numbers_noninteger:
            LOGGER.error("Row numbers contain non-integer values")

        if not (self.row_numbers_noninteger or self.row_numbers_missing) and not self.row_numbers_sequential:
            LOGGER.error("Row numbers not consecutive or do not start with 1")
        
        # Internal blank rows?
        if self.blank_rows:
            LOGGER.error("Data contains empty rows")

         # report on detected size
        LOGGER.info(f"Worksheet '{self.name}' contains {self.n_descriptors} descriptors, "
                    f"{self.n_row} data rows and {self.n_fields} fields")

        # reporting
        self.n_errors  = CH.counters['ERROR'] - self.start_errors
        if self.n_errors > 0:
            LOGGER.info(f'Dataframe contains {self.n_errors} errors')
        else:
            LOGGER.info('Dataframe formatted correctly')

    @loggerinfo_push_pop('Loading from worksheet')
    def load_from_worksheet(self, worksheet: worksheet, row_chunk_size=1000):
        """Populate a Dataworksheet instance from an openpyxl worksheet
        
        This method takes a fresh Dataworksheet instance and populates the
        details using the contents of a SAFE formatted Excel spreadsheet.

        To keep memory requirements low, any data in the worksheet is validated
        by loading sets of rows containing at most `row_chunk_size` rows. Note
        that this will only actually reduce memory use for openpyxl
        ReadOnlyWorksheets, as these load data on demand, but will still work
        on standard Worksheets.

        Args:
            worksheet: 
            row_chunk_size:
        """

        # TODO - openpyxl read-only mode provides row by row data access and used lazy
        # loading to reduce memory usage. Scanning all of a column will need to
        # load all data (which the old implementation using xlrd did, but which
        # we might be also improve to reduce memory overhead by chunking large sheets).
        # Either way, should use row by row ingestion.

        if self.fields_loaded:
            LOGGER.critical('Field metadata already loaded - use fresh instance.')
            return

        # get the data dimensions
        max_row = worksheet.max_row

        # trap completely empty worksheets
        if max_row == 0:
            LOGGER.error('Worksheet is empty')
            return

        # Read the field metadata first: 
        # - There should be rows starting with a string and then rows starting with
        #   a number, or a StopIteration at the end of the row iterator.
        ws_rows = worksheet.iter_rows(values_only=True)
        field_meta = []

        # While there are rows and while they start with a string cell, collect
        # field metadata
        for row in ws_rows:
            if isinstance(row[0], str):
                field_meta.append(row)
            else:
                break


        # Convert field meta to OrderedDict and validate
        field_meta = OrderedDict(((rw[0], rw[1:]) for rw in field_meta))
        self.validate_field_meta(field_meta)

        # Get an iterator on the data rows
        data_rows = worksheet.iter_rows(min_row=len(field_meta) + 1,
                                        values_only=True)

        # Load and validate chunks of data
        n_chunks = ((max_row -  self.n_descriptors) // row_chunk_size) + 1
        for chunk in range(n_chunks):
            
            # itertools.islice handles generators and StopIteration, and also
            # trap empty slices
            data = list(islice(data_rows, row_chunk_size))
            if data:
                self.validate_data_rows(data)
        
        # Finish up
        self.report()

    def _old_load_waste_code():
        """Code from old version that may need to work back into new approach"""
        row_sample = len(MANDATORY_DESCRIPTORS) + len(OPTIONAL_DESCRIPTORS) + 10
        
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
        self.field_name_row = n_field_meta

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
                 dataset: Dataset = None) -> None:
        
        """
        Base class to check field metadata and then test whether the data
        contents of a field are compatible with the metadata for the field.
        Largely used for checking data worksheets, but also by other methods.
        This base class defines the core checking functions as a base class and
        then subclasses implement data type specific checking.

        In order to work with the row by row data loading from openpyxl, the
        BaseField class does _not_ emit logging messages as they occur but holds
        a stack of messages created during initiation and as chunks of data are
        ingested. These logging messages can then be emitted using the report()
        method, along with any final checks, once all the data has been
        ingested.

        The base class also has a much more complex signature than the base
        functionality requires. Various subclasses need access to:

            * Dataset level information - extents, taxa and locations.
            * Dataworksheet level information - taxon fields

        Rather than having complex inheritance with changing subclass signatures
        and kwargs, the BaseField class makes all information available to all
        subclasses.

        Args: 
            meta: A dictionary of field metadata
            dwsh: A DataWorksheet instance, used to pass worksheet level
                information, which at the moment is just taxa fields.
            dataset: An Dataset instance, used to provide dataset level
                information, such as taxa, locations and summary, and update
                dataset level attributes such as data extents.
        """

        self.meta = meta
        self.dwsh = dwsh
        
        self.taxa = None
        self.locations = None
        self.summary = None

        # Shortcuts to main components
        self.dataset = dataset
        if self.dataset:
            self.taxa = getattr(self.dataset, 'taxa')
            self.locations = getattr(self.dataset, 'locations')
            self.summary = getattr(self.dataset, 'summary')

        # Attributes
        self.log_stack = []
        self.n_rows = 0
        self.n_na = 0
        self.n_blank = 0
        self.n_excel_error = 0
        self.bad_values = []
        self.bad_rows = [] # TODO check on implementation of this

        # Get a field name - either from the field_meta, or a column letter
        # from col_idx or 'Unknown'
        self.field_name = self.meta.get('field_name')

        if self.field_name is None:
            idx = self.meta.get('col_idx')
            if idx is None:
                self.field_name = 'Unknown'
            else:
                self.field_name = f"Column_{get_column_letter(idx)}"

        # Now check required descriptors - values are present, non-blank and 
        # are unpadded strings (currently all descriptors are strings).
        for dsc in self.required_descriptors:
            self._check_meta(dsc)
        
        # Specific check that field name is valid - the column letter codes are always
        # valid, so missing names won't trigger this.
        if not valid_r_name(self.field_name):
            self._log(f"Field name is not valid: {repr(self.field_name)}. "
                      "Common errors are spaces and non-alphanumeric "
                      "characters other than underscore and full stop")

        # TODO - rethink implementation? Quite specific behaviour unique to 
        #        some fields, so could implement as overloaded extra stub.
        if self.check_taxon_meta:
            self._check_taxon_meta()
        
        if self.check_interaction_meta:
            self._check_interaction_meta()

    def _log(self, msg, level=ERROR):
        """Helper function for adding to log stack"""
        self.log_stack.append((level, msg))

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
            self._log(f'{descriptor} descriptor missing')
            return False

        val = self.meta[descriptor]
        
        if blank_value(val):
            self._log(f'{descriptor} descriptor is blank')
            self.meta[descriptor] = None  # standardise whitestring to None
            return False
        elif not isinstance(val, str):
            self._log(f'{descriptor} descriptor is not a string: {repr(val)}')
            return False
        elif val != val.strip():
            self._log(f'{descriptor} descriptor has whitespace padding: {repr(val)}')
            self.meta[descriptor] = val.strip()
            return False
        else:
            return True

    def _check_taxon_meta(self):

        """
        Checks the taxonomic metadata of abundance and trait fields. This is
        more involved that the simple _check_meta(), because of the option to
        provide taxon_name or taxon_field descriptors.
        
        Returns:
            A boolean, with True showing no problems and False showing that
            warnings occurred.
        """

        # Are taxon_name and taxon_field provided and not blank: note use of
        # 'and' rather than '&' to allow missing descriptors to short cut

        tx_nm = self.meta.get('taxon_name')
        tx_fd = self.meta.get('taxon_field')
        tx_nm_prov = not blank_value(tx_nm)
        tx_fd_prov = not blank_value(tx_fd)

        if tx_nm_prov and tx_fd_prov:
            self._log('Taxon name and taxon field both provided, use one only')
            return False
        elif not tx_nm_prov and not tx_fd_prov:
            self._log("One of taxon name or taxon field must be provided")
            return False
        elif tx_nm_prov and self.taxa is None:
            self._log('Taxon name provided but no Taxa instance available', CRITICAL)
            return False
        elif tx_nm_prov and self.taxa.is_empty:
            self._log('Taxon name provided but no taxa loaded')
            return False
        elif tx_nm_prov and tx_nm not in self.taxa.taxon_names:
            self._log('Taxon name not found in the Taxa worksheet')
            return False
        elif tx_nm_prov:
            self.taxa.taxon_names_used.add(tx_nm)
            return True
        elif tx_fd_prov and self.dwsh is None:
            self._log(f"Taxon field provided but no dataworksheet provided for this field: {tx_fd}",
                      CRITICAL)
            return False
        elif tx_fd_prov and tx_fd not in self.dwsh.taxa_fields:
            self._log(f"Taxon field not found in this worksheet: {tx_fd}")
            return False
        else:
            return True

    def _check_interaction_meta(self):

        """
        Checks the taxonomic metadata of interaction fields. This is more
        involved that the simple _check_meta(), because of the option to provide
        taxon_name or taxon_field descriptors describing at least two taxonomic
        identifiers.

        Args:
            meta: A dictionary of metadata descriptors for the field
            taxa_fields: A list of Taxa fields in this worksheet.
        
        Returns:
            A boolean, with True showing no problems and False showing that
            warnings occurred.
        """

        # TODO - currently no checking for descriptions being present in
        #        the interaction information. Not sure how many old datasets
        #        would be affected, but this would be good to include.

        # Are interaction_name and/or interaction_field provided and not blank
        iact_nm = self.meta.get('interaction_name')
        iact_fd = self.meta.get('interaction_field')
        iact_nm_prov = not blank_value(iact_nm)
        iact_fd_prov = not blank_value(iact_fd)

        if not iact_nm_prov and not iact_fd_prov:
            self._log("At least one of interaction name or interaction field must be provided")
            return False
        elif iact_nm_prov and self.taxa is None:
            self._log('Interaction name provided but no Taxa instance available', CRITICAL)
            return False
        elif iact_nm_prov and self.taxa.is_empty:
            self._log('Interaction name provided but no taxa loaded')
            return False
        elif iact_fd_prov and self.dwsh is None:
            self._log(f"Interaction field provided but no dataworksheet provided for this field: {iact_fd}",
                      CRITICAL)
            return False

        if iact_nm_prov:
            # get the taxon names and descriptions from interaction name providers
            iact_nm_lab, iact_nm_desc = self._parse_levels(iact_nm)

            # add names to used taxa
            self.taxa.taxon_names_used.update(iact_nm_lab)

            # check they are found
            iact_nm_lab = IsInSet(iact_nm_lab, self.taxa.taxon_names)

            if not iact_nm_lab:
                self._log('Unknown taxa in interaction_name descriptor',
                          extra={'join': iact_nm_lab.failed})
                nm_check = False
            else:

                nm_check = True
            
            iact_nm_lab = iact_nm_lab.values
        else:
            iact_nm_lab = []
            iact_nm_desc = ()
            nm_check = True

        if iact_fd_prov:
            # check any field labels match to known taxon fields
            iact_fd_lab, iact_fd_desc = self._parse_levels(iact_fd)

            iact_fd_lab = IsInSet(iact_fd_lab, self.dwsh.taxa_fields)
            if not iact_fd_lab:
                self._log('Unknown taxon fields in interaction_field descriptor',
                          extra={'join': iact_fd_lab.failed})
                fd_check = False
            else:
                fd_check = True
            
            iact_fd_lab = iact_fd_lab.values
        else:
            iact_fd_lab = []
            iact_fd_desc = ()
            fd_check = True

        if len(iact_nm_lab + iact_fd_lab) < 2:
            self._log('At least two interacting taxon labels or fields must be identified')
            num_check = False
        else:
            num_check = True

        all_desc = iact_nm_desc + iact_fd_desc
        if None in all_desc:
            self._log('Label descriptions for interacting taxa incomplete or missing', 
                      WARNING)

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
            self._log('Extra colons in level description.')

        # standardise descriptions
        if all([pt == 1 for pt in n_parts]):
            parts = [[pt[0], None] for pt in parts]
        elif all([pt > 1 for pt in n_parts]):
            # truncate extra colons
            parts = [pt[0:2] for pt in parts]
        else:
            self._log('Provide descriptions for either all or none of the categories')
            parts = [pt[0:2] if len(pt) >= 2 else [pt[0], None] for pt in parts]

        level_labels, level_desc = zip(*parts)

        # - repeated labels?
        if len(set(level_labels)) < len(level_labels):
            self._log('Repeated level labels')

        # - check for numeric level names: integers would be more common
        #   but don't let floats sneak through either!
        level_labels = IsNotNumericString(level_labels)
        if not level_labels:
            self._log('Numeric level names not permitted')

        # Remove white space around the labels: simple spacing in the text
        # makes it easier to read and insisting on no space is unnecessary
        level_labels = [vl.strip() for vl in level_labels]

        return level_labels, level_desc

    @classmethod
    def field_type_map(cls):
        """This helper function returns a dictionary that maps field types for the
        a Field class and all nested subclasses to the class that handles them"""

        field_type_map = {}

        # Create a stack starting with the calling class
        classes = cls.__subclasses__() 

        # Pop subclasses while the stack is not empty
        while classes:
            # Get a field class
            current_field_class = classes.pop(0)

            # Add any subclasses from that class
            classes.extend(current_field_class.__subclasses__())
            
            # Extend the field map
            for ftype in current_field_class.field_types:
                field_type_map.update({ftype: current_field_class})

        return field_type_map


    def validate_data(self, data: list):
        """Validates a list of data provided for a field. 

        This method runs the common shared validation and returns the cleaned
        values, excluding blanks, NA and any Excel errors.

        The method is overloaded by subclasses to provide field specific
        testing. In this case using the following ensures that only the data
        that passes the common checks is subjected to extra testing:

            data = super().validate_data(data)
        """
        
        if self.no_validation:
            return

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
        """Report on field creation and data validation.

        This  method emits any logging messages accumulated during the creation
        of a field instance and data validation. The method also contains any final
        checking of data values across all validated data.
        """

        LOGGER.info(f'Checking field {self.field_name}')

        FORMATTER.push()

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

        FORMATTER.pop()


class CommentField(BaseField):

    field_types = ('comments', )
    no_validation = True


class ReplicateField(BaseField):

    field_types = ('replicate', 'id')


class NumericField(BaseField):
    """
    Subclass of BaseField to check for numeric data
    """
    field_types = ('numeric',)
    required_descriptors = MANDATORY_DESCRIPTORS + ['method', 'units']

    def validate_data(self, data: list):
        
        data = super().validate_data(data)
        
        numeric = IsNumber(data)

        if not numeric:
            self._log('Cells contain non-numeric values')


class CategoricalField(BaseField):
    """
    Subclass of BaseField to check for categorical data
    """

    field_types = ('categorical', 'ordered_categorical')
    required_descriptors = MANDATORY_DESCRIPTORS + ['levels']

    def __init__(self, meta: dict, dwsh: DataWorksheet = None, 
                 dataset: Dataset = None) -> None:

        super().__init__(meta, dwsh=dwsh, dataset=dataset)

        # Additional code to validate the levels metadata and store a set of
        # values
        self.level_labels = set()
        self.reported_levels = set()
        levels = meta.get('levels')
        
        if levels is not None and isinstance(levels, str) and not levels.isspace():
            level_labels, level_desc = self._parse_levels(levels)
            self.level_labels = set(level_labels)

    def validate_data(self, data: list):

        data = super().validate_data(data)

        # Now look for consistency: get the unique values reported in the
        # data, convert to unicode to handle checking of numeric labels and
        # then check the reported levels are a subset of the descriptors.
        # XLRD reads all numbers as floats, so coerce floats back to int
        data = IsString(data, keep_failed=False)
        if not data:
            self._log('Cells contain non-text values')

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


class TaxaField(BaseField):
    """
    Checks if all the values provided in a taxa field are found
    in the Taxa instance.
    """

    field_types = ('taxa',)

    def __init__(self, meta: dict, dwsh: DataWorksheet = None, 
                 dataset: Dataset = None) -> None:

        super().__init__(meta, dwsh=dwsh, dataset=dataset)

        if self.taxa is None:
            self._log('No taxon details provided for dataset')
        
        self.taxa_found = set()
    
    def validate_data(self, data: list):

        data = super().validate_data(data)

        data = IsString(data, keep_failed=False)

        if not data:
            self._log('Cells contain non-string values')

        if self.taxa is not None:
            self.taxa_found.update(data)
            self.taxa.taxon_names_used.update(data)
    
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


class LocationsField(BaseField):
    """
    Checks if all the values provided in a locations field are found
    in the Locations instance.
    """

    field_types = ('locations', 'location')

    def __init__(self, meta: dict, dwsh: DataWorksheet = None, 
                 dataset: Dataset = None) -> None:

        super().__init__(meta, dwsh=dwsh, dataset=dataset)

        if self.locations is None:
            self._log('No location details provided for dataset')
        
        self.locations_found = set()
    
    def validate_data(self, data: list):

        data = super().validate_data(data)

        data = IsLocName(data, keep_failed=False)

        if not data:
            self._log('Cells contain invalid location values')

        # Now should be strings and integer codes - convert to string
        # representations as used in the locations
        data = [str(v) for v in data]

        if self.locations is not None:
            self.locations_found.update(data)
            self.locations.locations_used.update(data)
    
    def report(self):

        super().report()

        # TODO - not sure about this - no other fields test for emptiness?
        if self.locations_found == set():
            LOGGER.error('No locations loaded')
            return
        
        if self.locations is not None:
            extra_locs = self.locations_found.difference(self.locations.locations)

            if extra_locs:
                LOGGER.error('Includes unreported locations: ',
                                extra={'join': extra_locs})

            # add the found taxa to the list of taxa used
            self.locations.locations_used.update(self.locations_found)


class GeoField(BaseField):

    field_types = ('latitude', 'longitude')

    def __init__(self, meta: dict, dwsh: DataWorksheet = None,
                 dataset: Dataset = None) -> None:
        
        super().__init__(meta, dwsh=dwsh, dataset=dataset)

        if self.dataset is None:
            self._log('No dataset object provided - cannot update extents')

        self.min = None
        self.max = None
        
    def validate_data(self, data: list):
        """Testing of the data range is deferred to report() to avoid
        triggering logging from the Extent objects until data loading is
        completed
        """
        data = super().validate_data(data)

        data = IsNumber(data, keep_failed=False)

        if not data:
            self._log('Field contains non-numeric data')

            if any([RE_DMS.search(str(dt)) for dt in data.failed]):
                self._log('Possible degrees minutes and seconds formatting? Use decimal degrees', 
                          WARNING)
        
        if data.values:
            self.min = min(data.values + [self.min]) if self.min else min(data.values)
            self.max = max(data.values + [self.max]) if self.max else max(data.values)

    def report(self):

        super().report()

        if self.min is None or self.max is None:
            return()
        
        if self.dataset is not None:
            if self.meta['field_type'] == 'latitude':
                self.dataset.latitudinal_extent.update([self.min, self.max])
            elif self.meta['field_type'] == 'longitude':
                self.dataset.longitudinal_extent.update([self.min, self.max])


class NumericTaxonField(NumericField):

    """Checks abundance and numeric trait fields, which are just  NumericFields
    with taxon reporting requirements turned on."""

    # BREAKING CHANGE - abundance did not previosly require 'units' as metadata,
    # which is stupid in retrospect. Individuals? Individuals per m2? Individuals
    # per hour? I do not know how many existing datasets will fall foul of this.

    field_types = ('abundance', 'numeric trait')
    check_taxon_meta = True


class CategoricalTaxonField(CategoricalField):

    """Checks categorical trait fields, which are just CategoricalFields
    with taxon reporting requirements turned on."""

    field_types = ('categorical trait', 'ordered categorical trait')
    check_taxon_meta = True


class NumericInteractionField(NumericField):

    """Checks numeric interaction fields, which are just  NumericFields
    with interaction reporting requirements turned on."""

    field_types = ('numeric interaction', )
    check_interaction_meta = True


class CategoricalInteractionField(CategoricalField):

    """Checks categorical interaction fields, which are just CategoricalFields
    with interaction reporting requirements turned on."""

    field_types = ('categorical interaction', )
    check_interaction_meta = True


class TimeField(BaseField):

    """
    Time field validation is primarily about checking the consistency of the
    provided data. The function accepts either as ISO strings in a consistent
    time format or datetime.time objects. The parsing of timestrings follows
    ISO8601, which includes quite a range of valid inputs: '11:12:13' and
    '111213' are both valid.

    One problem here is that - when stored as Excel datetime cells - dates and
    times in Excel are poorly typed: essentially as a number with a meaning that
    depends in part on the cell format. The openpyxl package returns
    datetime.time objects for cells with time formats.
    """

    field_types = ('time', )

    def __init__(self, meta: dict, dwsh: DataWorksheet = None, 
                 dataset: Dataset = None) -> None:
        
        super().__init__(meta, dwsh=dwsh, dataset=dataset)

        # Defaults
        self.first_data_class_set = None
        self.consistent_class = True
        self.expected_class = True

        self.bad_strings = []
        self.min = None 
        self.max = None

    def validate_data(self, data: list):
        
        data = super().validate_data(data)

        # TODO: report row number of issues - problem of mismatch between
        # self.n_rows (_all_ rows) and index in data (non-NA rows)

        # Check for consistent class formatting, using first row. Do not try and
        # validate further when data is not consistently formatted.
        if self.consistent_class and self.expected_class:

            cell_types = [type(dt) for dt in data]
            cell_type_set = set(cell_types)

            # Are all the values of expected types?
            if not {datetime.time, str}.issuperset(cell_type_set):
                self.expected_class = False
                return
            
            # Are the values internally consistent and consistent with
            # previously loaded data
            if len(cell_type_set) > 1:
                self.consistent_class = False
                return
            elif self.first_data_class_set is None:
                self.first_data_class_set = cell_type_set
            elif self.first_data_class_set != cell_type_set:
                self.consistent_class = False
                return
        else:
            return

        # Value checking only happens while data are consistently formatted
        # There is no need to check time objects passed in, just time formatted
        # strings
        if self.first_data_class_set == {str}:
            
            for val in data:
                try:
                    _ = parser.isoparser().parse_isotime(val)
                except ValueError:
                    self.bad_strings.append(val)
            
    
    def report(self):
        
        super().report()

        if not self.expected_class:
            LOGGER.error('Time data include values neither ISO string nor time formatted')
            return

        if not self.consistent_class:
            LOGGER.error('Time data mixes ISO string and time formatted rows')
            return

        if self.bad_strings:
            LOGGER.error('ISO time strings contain badly formatted values: e.g.',
                         extra={'join': self.bad_strings[:5]})


class DatetimeField(BaseField):

    """
    Date and datetime field validation is primarily about checking the
    consistency of the provided data. The function accepts either as ISO strings
    in a consistent datetime or date format or datetime.datetime objects. The
    parsing of timestrings follows ISO8601, which includes quite a range of
    valid inputs.

    One problem here is that - when stored as Excel datetime cells - dates and
    times in Excel are poorly typed: essentially as a number with a meaning that
    depends in part on the cell format. The openpyxl package returns
    datetime.datetime objects for cells with datetime or date formats.
    """

    field_types = ('datetime', 'date')

    def __init__(self, meta: dict, dwsh: DataWorksheet = None, 
                 dataset: Dataset = None) -> None:
        
        super().__init__(meta, dwsh=dwsh, dataset=dataset)

        if self.dataset is None:
            self._log('No dataset object provided - cannot update extents')

        # Defaults
        self.first_data_class_set = None
        self.consistent_class = True
        self.expected_class = True
        self.all_midnight = True

        self.bad_strings = []
        self.min = None 
        self.max = None

    def validate_data(self, data: list):
        
        data = super().validate_data(data)

        # TODO: report row number of issues - problem of mismatch between
        # self.n_rows (_all_ rows) and index in data (non-NA rows)

        # Check for consistent class formatting, using first row. Do not try and
        # validate further when data is not consistently formatted.
        if self.consistent_class and self.expected_class:

            cell_types = [type(dt) for dt in data]
            cell_type_set = set(cell_types)

            # Are all the values of expected types?
            if not {datetime.datetime, str}.issuperset(cell_type_set):
                self.expected_class = False
                return
            
            # Are the values internally consistent and consistent with
            # previously loaded data
            if len(cell_type_set) > 1:
                self.consistent_class = False
                return
            elif self.first_data_class_set is None:
                self.first_data_class_set = cell_type_set
            elif self.first_data_class_set != cell_type_set:
                self.consistent_class = False
                return
        else:
            return

        # Value checking only happens while data are consistently formatted
        # Look for string formatting, but also presence/absence of non-zero
        # time components.
        midnight = datetime.time(0,0)

        if self.first_data_class_set == {str}:
            parsed_strings = []
            for val in data:
                try:
                    val_parse = parser.isoparse(val)
                    if val_parse.time() != midnight:
                        self.all_midnight = False
                    parsed_strings.append(val_parse)
                except ValueError:
                    self.bad_strings.append(val)
            
            # standardised to datetime objects
            data = parsed_strings

        elif self.first_data_class_set == {datetime.datetime}:

            for val in data:
                if val.time() != midnight:
                    self.all_midnight = False

        # Update range for extents
        if self.min is None:
            self.min = min(data)
        else:
            self.min = min([self.min] + data)
        
        if self.max is None:
            self.max = min(data)
        else:
            self.max = max([self.max] + data)
        
    def report(self):
        
        super().report()

        # INconsistent and bad data classes
        if not self.expected_class:
            LOGGER.error('Date or datetime data include values neither ISO string nor time formatted')
            return

        if not self.consistent_class:
            LOGGER.error('Date or datetime data mixes ISO string and time formatted rows')
            return

        # Bad string formats
        if self.bad_strings:
            LOGGER.error('ISO datetime strings contain badly formatted values: e.g.',
                         extra={'join': self.bad_strings[:5]})

        # Datetime vs Date checking - Datetimes _could_ all be at midnight, so not an error
        # but does look odd.
        if self.meta['field_type'] == 'datetime' and self.all_midnight:
            LOGGER.warning('Field is of type datetime, but only reports dates')
        if self.meta['field_type'] == 'date' and not self.all_midnight:
            LOGGER.error('Field is of type date, but includes time data')

        # Update extent if possible - note that inheritance means that isinstance
        # in extent.Extent is not succesfully testing for datetime.datetime rather
        # than set datatype of datetime.date
        if self.dataset is not None:
            self.dataset.temporal_extent.update([self.min.date(), self.max.date()])


class FileField(BaseField):

    """
    Checks file fields. The data values need to match to an external file or can
    be contained within an archive file provided in the 'file_container'
    descriptor.
    """

    field_types = ('file',)

    def __init__(self, meta: dict, dwsh: DataWorksheet = None, 
                 dataset: Dataset = None) -> None:
        
        super().__init__(meta, dwsh=dwsh, dataset=dataset)

        self.unknown_file_names = set()

        # Check whether filename testing is possible
        if self.summary is None:
            self._log('No Summary instance provided - cannot check file fields', 
                      level=CRITICAL)
            return
        if self.summary.external_files is None:
            self._log('No external files listed in Summary')
            return
        
        self.external_names = {ex['file'] for ex in self.summary.external_files}

        # Look for a file container metadata value, pointing to a single file
        # containing all of the listed values (zip etc)
        self.file_container = self.meta.get('file_container')

        if self.file_container is not None:
            if self.file_container not in self.external_names:
                self._log(f"Field file_container value not found in external files: {self.file_container}")


    def validate_data(self, data: list):

        data = super().validate_data(data)

        # If the files are listed in external and not provided in a file
        # container then check they are all present 
        if not self.file_container:
            unknown_files = set(data) - self.external_names
            if unknown_files:
                self.unknown_file_names.update(unknown_files)

    def report(self):
        
        super().report()

        if self.unknown_file_names:
            LOGGER.error("Field contains external files not provided in Summary: ",
                         extra={'join': self.unknown_file_names})


class EmptyField():
    """This class mocks the interface of BaseField and is only used to keep
    track of trailing empty fields in loaded field data. It does not inherit
    from BaseField, so is never included in the BaseField.field_type_map."""

    def __init__(self, meta: dict) -> None:
        
        self.meta = meta
        self.empty = True

    def validate_data(self, data):

        blank = [blank_value(val) for val in data]

        if not all(blank):
            self.empty = False
    
    def report(self):

        # Get a field name - either a column letter from col_idx if set or 'Unknown'
        idx = self.meta.get('col_idx')
        if idx is None:
            self.field_name = 'Unknown'
        else:
            self.field_name = f"Column_{get_column_letter(idx)}"
        
        if not self.empty:
            LOGGER.info(f'Checking field {self.field_name}')
            FORMATTER.push()
            LOGGER.error('Trailing field with no descriptors contains data.')
            FORMATTER.pop()




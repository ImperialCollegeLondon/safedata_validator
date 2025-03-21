"""Test of logging and exception behaviour of DataWorksheet."""

from collections import OrderedDict
from logging import CRITICAL, ERROR, INFO, WARNING

import openpyxl
import pytest

from safedata_validator.field import Dataset, DataWorksheet


@pytest.mark.parametrize(
    "sheet_meta,expected_log",
    [
        (
            {
                "name": "DF",
                "title": "My data table",
                "description": "This is a test data worksheet",
                "external": None,
            },
            tuple(),  # No logging emitted
        ),
    ],
)
def test_DataWorksheet_init(caplog, fixture_field_meta, sheet_meta, expected_log):
    """Testing init errors - field_meta tested elsewhere.

    TODO - currently no errors throw  in init - but might be added
    """

    DataWorksheet(sheet_meta, fixture_field_meta)

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "field_meta,expected_log",
    [
        (
            OrderedDict(
                field_type=["numeric", "numeric", "numeric"],
                description=["a", "b", "c"],
                field_name=["a", "b", "c"],
            ),
            ((INFO, "Validating field metadata"),),
        ),
        (
            OrderedDict(
                [
                    ("field_type ", ["numeric", "numeric", "numeric"]),
                    (" description", ["a", "b", "c"]),
                    ("  field_name", ["a", "b", "c"]),
                ]
            ),
            (
                (INFO, "Validating field metadata"),
                (ERROR, "Whitespace padding in descriptor names"),
            ),
        ),
        (
            OrderedDict(
                [
                    ("field_type", ["numeric", "numeric", "numeric"]),
                    ("field_name", ["a", "b", "c"]),
                ]
            ),
            ((INFO, "Validating field metadata"),),
        ),
        (
            OrderedDict(
                [
                    ("field_type", ["numeric", "numeric", "numeric"]),
                    ("description", ["a", "b", "c"]),
                    ("whatdoesthisdo", ["a", "b", "c"]),
                    ("field_name", ["a", "b", "c"]),
                ]
            ),
            (
                (INFO, "Validating field metadata"),
                (WARNING, "Unknown field descriptors"),
            ),
        ),
        (
            OrderedDict(
                [
                    ("field_type", ["numeric", "numeric", "numeric"]),
                    ("field_name", ["a", "b", "c"]),
                    ("description", ["a", "b", "c"]),
                ]
            ),
            (
                (INFO, "Validating field metadata"),
                (ERROR, "Field_name row is not the last descriptor"),
            ),
        ),
        (
            OrderedDict(
                [
                    ("field_type", ["numeric", "numeric", "numeric"]),
                    ("description", ["a", "b", "c"]),
                    ("field_name", ["a", "b", "a"]),
                ]
            ),
            (
                (INFO, "Validating field metadata"),
                (ERROR, "Field names duplicated"),
            ),
        ),
        (
            OrderedDict(
                [
                    ("field_type", ["numeric", "bad type", "numeric"]),
                    ("description", ["a", "b", "c"]),
                    ("field_name", ["a", "b", "c"]),
                ]
            ),
            (
                (INFO, "Validating field metadata"),
                (ERROR, "Unknown field types"),
            ),
        ),
        (
            OrderedDict(
                [
                    ("field_type", ["numeric", "numeric", None]),  # trailing empty
                    ("description", ["a", "b", None]),
                    ("field_name", ["a", "b", None]),
                ]
            ),
            ((INFO, "Validating field metadata"),),
        ),
        (
            OrderedDict(
                [
                    ("field_type", ["numeric", None, "numeric"]),  # internal empty
                    ("description", ["a", None, "b"]),
                    ("field_name", ["a", None, "b"]),
                ]
            ),
            ((INFO, "Validating field metadata"),),
        ),
    ],
)
def test_DataWorksheet_validate_field_meta(
    caplog, fixture_dataworksheet, field_meta, expected_log
):
    """Testing errors created when field metadata is validated."""
    fixture_dataworksheet.validate_field_meta(field_meta)

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "load_field_meta,data_rows,expected_log",
    [
        (
            True,
            [
                [1, 1, 2, 3],
                [2, 1, 2, 3],
                [3, 1, 2, 3],
            ],
            ((INFO, "Validating field metadata"),),
        ),
        (
            False,
            [
                [1, 1, 2, 3],
                [2, 1, 2, 3],
                [3, 1, 2, 3],
            ],
            ((CRITICAL, "No fields defined - use validate_field_meta"),),
        ),
        (
            True,
            [
                [1, 1, 2, 3, 4],
                [2, 1, 2, 3],
                [3, 1, 2, 3],
            ],
            (
                (INFO, "Validating field metadata"),
                (CRITICAL, "Data rows of unequal length - cannot validate"),
            ),
        ),
        (
            True,
            [
                [1, 1, 2, 3, 4],
                [2, 1, 2, 3, 4],
                [3, 1, 2, 3, 4],
            ],
            (
                (INFO, "Validating field metadata"),
                (
                    CRITICAL,
                    "Data rows not of same length as field metadata - cannot validate",
                ),
            ),
        ),
        (
            True,
            [],
            (
                (INFO, "Validating field metadata"),
                (CRITICAL, "Empty data_rows passed to validate_data_rows"),
            ),
        ),
    ],
)
def test_DataWorksheet_validate_data_rows(
    caplog,
    fixture_dataworksheet,
    fixture_field_meta,
    load_field_meta,
    data_rows,
    expected_log,
):
    """Testing errors created as data rows are added."""
    caplog.clear()

    if load_field_meta:
        fixture_dataworksheet.validate_field_meta(fixture_field_meta)

    fixture_dataworksheet.validate_data_rows(data_rows)

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "field_meta,data_rows,expected_log",
    [
        (
            OrderedDict(
                field_type=["id", "id", "id"],
                description=["a", "b", "c"],
                field_name=["a", "b", "c"],
            ),
            [
                [1, 1, 2, 3],
                [2, 1, 2, 3],
                [3, 1, 2, 3],
            ],
            (
                (INFO, "Validating field metadata"),
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (INFO, "Worksheet 'DF' contains 3 descriptors, 3 data rows"),
                (INFO, "Dataframe formatted correctly"),
            ),
        ),
        (
            OrderedDict(
                field_type=["id", "id", "id", None],
                description=["a", "b", "c", None],
                field_name=["a", "b", "c", None],
            ),
            [
                [1, 1, 2, 3, None],
                [2, 1, 2, 3, None],
                [3, 1, 2, 3, None],
            ],
            (
                (INFO, "Validating field metadata"),
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (INFO, "Worksheet 'DF' contains 3 descriptors, 3 data rows"),
                (INFO, "Dataframe formatted correctly"),
            ),
        ),
        (
            OrderedDict(
                field_type=["id", "id", "id", None],
                description=["a", "b", "c", None],
                field_name=["a", "b", "c", None],
            ),
            [
                [1, 1, 2, 3, None],
                [2, 1, 2, 3, "argh"],
                [3, 1, 2, 3, None],
            ],
            (
                (INFO, "Validating field metadata"),
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (INFO, "Checking field Column_D"),
                (ERROR, "Trailing field with no descriptors contains data"),
                (INFO, "Worksheet 'DF' contains 3 descriptors, 3 data rows"),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
        (
            OrderedDict(
                field_type=["id", "id", None, "id"],
                description=["a", "b", None, "c"],
                field_name=["a", "b", None, "c"],
            ),
            [
                [1, 1, 2, None, 3],
                [2, 1, 2, "argh", 3],
                [3, 1, 2, None, 3],
            ],
            (
                (INFO, "Validating field metadata"),
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field Column_C"),
                (ERROR, "field_type descriptor is blank"),
                (ERROR, "description descriptor is blank"),
                (ERROR, "field_name descriptor is blank"),
                (ERROR, "2 cells are blank or contain only whitespace text"),
                (INFO, "Checking field c"),
                (INFO, "Worksheet 'DF' contains 3 descriptors, 3 data rows"),
                (INFO, "Dataframe contains 4 errors"),
            ),
        ),
    ],
)
def test_DataWorksheet_empty_meta(
    caplog, fixture_dataworksheet, field_meta, data_rows, expected_log
):
    """Testing errors created as data rows are added."""
    caplog.clear()

    fixture_dataworksheet.validate_field_meta(field_meta)
    fixture_dataworksheet.validate_data_rows(data_rows)
    fixture_dataworksheet.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data_rows,populate_external,expected_log",
    [
        (
            [
                [1, 1, 2, 3],
                [2, 1, 2, 3],
                [3, 1, 2, 3],
            ],
            False,
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 3 data rows and 3 fields",
                ),
                (INFO, "Dataframe formatted correctly"),
            ),
        ),
        (
            [],
            False,
            (
                (ERROR, "No data passed for validation"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 0 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
        (
            [],
            True,
            (
                (INFO, "Data table description associated with external file"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 0 data rows and 3 fields",
                ),
                (INFO, "Dataframe formatted correctly"),
            ),
        ),
        (
            [
                [11, 1, 2, 3],
                [12, 1, 2, 3],
                [13, 1, 2, 3],
            ],
            False,
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (ERROR, "Row numbers not consecutive or do not start with 1"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 3 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
        (
            [
                [1, 1, 2, 3],
                [2, 1, 2, 3],
                [None, 1, 2, 3],
            ],
            False,
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (ERROR, "Missing row numbers in data"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 3 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
        (
            [
                [1, 1, 2, 3],
                [2, 1, 2, 3],
                ["3", 1, 2, 3],
            ],
            False,
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (ERROR, "Row numbers contain non-integer values"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 3 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
        (
            [[1, 1, 2, 3], [2, 1, 2, 3], [3, 1, 2, 3]],
            False,
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 3 data rows and 3 fields",
                ),
                (INFO, "Dataframe formatted correctly"),
            ),
        ),
        (
            [
                [1, 1, 2, 3],
                [2, 1, 2, 3],
                [None, None, None, None],
                [3, 1, 2, 3],
                [None, None, None, None],
            ],
            False,
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (ERROR, "Data contains empty rows"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 3 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
    ],
)
def test_DataWorksheet_report(
    caplog,
    fixture_dataworksheet,
    fixture_field_meta,
    data_rows,
    populate_external,
    expected_log,
):
    """Testing report level errors emission of field errors and dataset level issues."""

    if populate_external:
        fixture_dataworksheet.external = "file_provided.sql"

    fixture_dataworksheet.validate_field_meta(fixture_field_meta)
    caplog.clear()

    if data_rows:
        fixture_dataworksheet.validate_data_rows(data_rows)

    fixture_dataworksheet.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


# TODO - think about the _order_ of different rownumber errors and which
#        get reported when: 1,2,3,4,5,6,9,None,None - both non_seq and missing
#        useful information for user


@pytest.mark.parametrize(
    "data_rows,expected_log",
    [
        (
            [
                [
                    [1, 1, 2, 3],
                    [2, 1, 2, 3],
                    [3, 1, 2, 3],
                ],
                [
                    [4, 1, 2, 3],
                    [5, 1, 2, 3],
                    [6, 1, 2, 3],
                ],
            ],
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 6 data rows and 3 fields",
                ),
                (INFO, "Dataframe formatted correctly"),
            ),
        ),
        (
            [
                [
                    [1, 1, 2, 3],
                    [2, 1, 2, 3],
                    [3, 1, 2, 3],
                ],
                [
                    [7, 1, 2, 3],
                    [8, 1, 2, 3],
                    [9, 1, 2, 3],
                ],
            ],
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (ERROR, "Row numbers not consecutive or do not start with 1"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 6 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
        (
            [
                [
                    [1, 1, 2, 3],
                    [2, 1, 2, 3],
                    [3, 1, 2, 3],
                ],
                [
                    [4, 1, 2, 3],
                    [5, 1, 2, 3],
                    [
                        None,
                        None,
                        None,
                        None,
                    ],
                ],
            ],
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 5 data rows and 3 fields",
                ),
                (INFO, "Dataframe formatted correctly"),
            ),
        ),
        (
            [
                [
                    [1, 1, 2, 3],
                    [2, 1, 2, 3],
                    [3, 1, 2, 3],
                ],
                [
                    [4, 1, 2, 3],
                    [
                        None,
                        None,
                        None,
                        None,
                    ],  # internal blank row within block
                    [5, 1, 2, 3],
                ],
            ],
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (ERROR, "Data contains empty rows"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 5 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
        (
            [
                [
                    [1, 1, 2, 3],
                    [2, 1, 2, 3],
                    [
                        None,
                        None,
                        None,
                        None,
                    ],  # trailing blank row, but with data then loaded
                ],
                [
                    [3, 1, 2, 3],
                    [4, 1, 2, 3],
                    [5, 1, 2, 3],
                ],
            ],
            (
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (ERROR, "Data contains empty rows"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 5 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
    ],
)
def test_DataWorksheet_report_multi_load(
    caplog, fixture_dataworksheet, fixture_field_meta, data_rows, expected_log
):
    """Testing report level errors when more than one set of data loaded."""

    fixture_dataworksheet.validate_field_meta(fixture_field_meta)
    caplog.clear()

    for chunk in data_rows:
        fixture_dataworksheet.validate_data_rows(chunk)

    fixture_dataworksheet.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data_rows,populate_external,expected_log",
    [
        (
            [
                [1, 1, 2, 3],
                [2, 1, 2, 3],
                [3, 1, 2, 3],
            ],
            False,
            (
                (INFO, "Checking worksheet 'DF'"),
                (INFO, "Validating field metadata"),
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 3 data rows and 3 fields",
                ),
                (INFO, "Dataframe formatted correctly"),
            ),
        ),
        (
            [],
            False,
            (
                (INFO, "Checking worksheet 'DF'"),
                (INFO, "Validating field metadata"),
                (ERROR, "No data passed for validation"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 0 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
        (
            [],
            True,
            (
                (INFO, "Checking worksheet 'DF'"),
                (INFO, "Validating field metadata"),
                (INFO, "Data table description associated with external file"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 0 data rows and 3 fields",
                ),
                (INFO, "Dataframe formatted correctly"),
            ),
        ),
        (
            [
                [11, 1, 2, 3],
                [12, 1, 2, 3],
                [13, 1, 2, 3],
            ],
            False,
            (
                (INFO, "Checking worksheet 'DF'"),
                (INFO, "Validating field metadata"),
                (INFO, "Validating field data"),
                (INFO, "Checking field a"),
                (INFO, "Checking field b"),
                (INFO, "Checking field c"),
                (ERROR, "Row numbers not consecutive or do not start with 1"),
                (
                    INFO,
                    "Worksheet 'DF' contains 5 descriptors, 3 data rows and 3 fields",
                ),
                (INFO, "Dataframe contains 1 errors"),
            ),
        ),
    ],
)
def test_DataWorksheet_load_from_worksheet(
    caplog,
    fixture_dataworksheet,
    fixture_field_meta,
    data_rows,
    populate_external,
    expected_log,
):
    """Testing basics of load_from_worksheet on programmatically created worksheet."""

    # Create a worksheet programmatically
    wb = openpyxl.Workbook()
    ws = wb.active

    # Add field meta rows
    field_meta = [[ky, *vl] for ky, vl in fixture_field_meta.items()]
    for rw in field_meta:
        ws.append(rw)

    # Add data
    for rw in data_rows:
        ws.append(rw)

    # Test empty worksheet handling
    if populate_external:
        fixture_dataworksheet.external = "file_provided.sql"

    # Try loading worksheet
    fixture_dataworksheet.load_from_worksheet(ws)

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "example_excel_files, wsheet, n_errors",
    [
        ("good", "DF", 0),
        ("good", "Incidence", 0),
        ("bad", "DF", 41),
        ("bad", "Incidence", 8),
    ],
    indirect=["example_excel_files"],  # take actual params from fixture
)
def test_DataWorksheet_load_from_file(
    caplog, fixture_resources, example_excel_files, wsheet, n_errors
):
    """Test loading a dataworksheet from file.

    This duplicates a lot of Dataset.load_from_workbook
    """

    # Load the taxa and locations
    ds = Dataset(fixture_resources)
    ds.taxa.gbif_taxa.load(example_excel_files["Taxa"])
    ds.locations.load(example_excel_files["Locations"])

    caplog.clear()

    dws = DataWorksheet(
        {
            "name": "DF",
            "title": "My data table",
            "description": "This is a test data worksheet",
            "external": None,
        },
        dataset=ds,
    )

    dws.load_from_worksheet(example_excel_files[wsheet])

    assert dws.n_errors == n_errors

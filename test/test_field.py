"""General tests on the various field types."""
from collections import OrderedDict
from datetime import date, datetime, time
from logging import CRITICAL, ERROR, INFO, WARNING

import pytest

from safedata_validator.field import (
    BaseField,
    CategoricalField,
    CategoricalInteractionField,
    CategoricalTaxonField,
    CommentField,
    DataWorksheet,
    DatetimeField,
    EmptyField,
    FileField,
    GeoField,
    LocationsField,
    NumericField,
    NumericInteractionField,
    NumericTaxonField,
    ReplicateField,
    TaxaField,
    TimeField,
)

# Checking the helper and private methods


@pytest.mark.parametrize(
    "data, exp_log, exp_levels, exp_desc",
    [
        (
            "A simple string",
            ((INFO, "Checking field testy"),),
            ["A simple string"],
            [None],
        ),
        (
            "A simple string;Another string",
            ((INFO, "Checking field testy"),),
            ["A simple string", "Another string"],
            [None, None],
        ),
        (
            "A simple string:description;Another string:description",
            ((INFO, "Checking field testy"),),
            ["A simple string", "Another string"],
            ["description", "description"],
        ),
        (
            "level1:A description: with a colon in it",
            (
                (INFO, "Checking field testy"),
                (ERROR, "Extra colons in level description."),
            ),
            ["level1"],
            ["A description"],
        ),
        (
            "level1;level2:unbalanced descriptions",
            (
                (INFO, "Checking field testy"),
                (
                    ERROR,
                    "Provide descriptions for either all or none of the categories",
                ),
            ),
            ["level1", "level2"],
            [None, "unbalanced descriptions"],
        ),
        (
            "level1;level1",
            ((INFO, "Checking field testy"), (ERROR, "Repeated level labels")),
            ["level1", "level1"],
            [None, None],
        ),
        (
            "1;2",
            (
                (INFO, "Checking field testy"),
                (ERROR, "Numeric level names not permitted"),
            ),
            ["1", "2"],
            [None, None],
        ),
    ],
)
def test_parse_levels(caplog, data, exp_log, exp_levels, exp_desc):
    """Testing behaviour of the NumericField class in using _validate_data."""

    fld = BaseField(
        {
            "field_name": "testy",
            "description": "a test",
            "field_type": "irrelevant for this test",
        }
    )

    obs_lev, obs_desc = fld._parse_levels(data)
    fld.report()

    assert len(exp_log) == len(caplog.records)

    assert all([exp[0] == rec.levelno for exp, rec in zip(exp_log, caplog.records)])
    assert all([exp[1] in rec.message for exp, rec in zip(exp_log, caplog.records)])

    assert len(obs_lev) == len(exp_levels)
    assert all([lo == le for lo, le in zip(obs_lev, exp_levels)])

    assert len(obs_desc) == len(exp_desc)
    assert all([do == de for do, de in zip(obs_desc, exp_desc)])


# TODO -  test errors on empty Taxa


@pytest.mark.parametrize(
    "tx_meta, has_taxa_object, has_dwsh_object, expected_log",
    [
        (
            dict(taxon_name="foo", taxon_field="bar"),
            False,
            False,
            (
                (INFO, "Checking field tx_field"),
                (ERROR, "Taxon name and taxon field both provided, use one only"),
            ),
        ),
        (
            dict(),
            False,
            False,
            (
                (INFO, "Checking field tx_field"),
                (ERROR, "One of taxon name or taxon field must be provided"),
            ),
        ),
        (
            dict(taxon_name="foo"),
            False,
            False,
            (
                (INFO, "Checking field tx_field"),
                (CRITICAL, "Taxon name provided but no Taxa instance available"),
            ),
        ),
        (
            dict(taxon_name="foo"),
            True,
            False,
            (
                (INFO, "Checking field tx_field"),
                (ERROR, "Taxon name not found in the Taxa worksheet"),
            ),
        ),
        (dict(taxon_name="C_born"), True, False, ((INFO, "Checking field tx_field"),)),
        (
            dict(taxon_field="not_provided"),
            False,
            False,
            (
                (INFO, "Checking field tx_field"),
                (
                    CRITICAL,
                    "Taxon field provided but no dataworksheet provided for this field",
                ),
            ),
        ),
        (
            dict(taxon_field="not_provided"),
            False,
            True,
            (
                (INFO, "Checking field tx_field"),
                (ERROR, "Taxon field not found in this worksheet"),
            ),
        ),
        (
            dict(taxon_field="my_taxon_field"),
            False,
            True,
            ((INFO, "Checking field tx_field"),),
        ),
    ],
)
def test_check_taxon_meta(
    caplog, fixture_dataset, tx_meta, has_taxa_object, has_dwsh_object, expected_log
):
    """Testing the use of the BaseField._check_taxon_meta() method."""

    # Set up what information is available for taxon field validation
    ds_obj = fixture_dataset if has_taxa_object else None

    dwsh = DataWorksheet(
        {
            "name": "DF",
            "title": "My data table",
            "description": "This is a test data worksheet",
        }
    )
    dwsh.taxa_fields = ["my_taxon_field"]
    dwsh_obj = dwsh if has_dwsh_object else None

    # Technically, this violates the field_name last requirement, but that is
    # enforced at the dataworksheet level, not the field level.
    field_meta = OrderedDict(
        field_type="numeric", description="description", field_name="tx_field"
    )
    field_meta.update(tx_meta)

    fld = BaseField(field_meta, dataset=ds_obj, dwsh=dwsh_obj)

    caplog.clear()

    # Test the logging from this private method.
    fld._check_taxon_meta()
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "iact_meta, has_taxa_object, has_dwsh_object, expected_log",
    [
        (
            dict(),
            False,
            False,
            (
                (INFO, "Checking field iact_field"),
                (
                    ERROR,
                    "At least one of interaction name or interaction field must be"
                    " provided",
                ),
            ),
        ),
        (
            dict(interaction_name="foo:predator"),
            False,
            False,
            (
                (INFO, "Checking field iact_field"),
                (CRITICAL, "Interaction name provided but no Taxa instance available"),
            ),
        ),
        (
            dict(interaction_name="foo:predator"),
            True,
            False,
            (
                (INFO, "Checking field iact_field"),
                (ERROR, "Unknown taxa in interaction_name descriptor"),
                (
                    ERROR,
                    "At least two interacting taxon labels or fields must be "
                    "identified",
                ),
            ),
        ),
        (
            dict(interaction_name="C_born:predator"),
            True,
            False,
            (
                (INFO, "Checking field iact_field"),
                (
                    ERROR,
                    "At least two interacting taxon labels or fields must be"
                    " identified",
                ),
            ),
        ),
        (
            dict(
                interaction_name="C_born:predator;V_salv:prey"
            ),  # Biologically not very plausible
            True,
            False,
            ((INFO, "Checking field iact_field"),),
        ),
        (
            dict(interaction_field="not_provided:predator"),
            False,
            False,
            (
                (INFO, "Checking field iact_field"),
                (
                    CRITICAL,
                    "Interaction field provided but no dataworksheet provided for this"
                    " field",
                ),
            ),
        ),
        (
            dict(interaction_field="foo:predator"),
            False,
            True,
            (
                (INFO, "Checking field iact_field"),
                (ERROR, "Unknown taxon fields in interaction_field descriptor"),
                (
                    ERROR,
                    "At least two interacting taxon labels or fields must be "
                    "identified",
                ),
            ),
        ),
        (
            dict(interaction_field="predator:eats things"),
            False,
            True,
            (
                (INFO, "Checking field iact_field"),
                (
                    ERROR,
                    "At least two interacting taxon labels or fields must be "
                    "identified",
                ),
            ),
        ),
        (
            dict(interaction_field="predator:eats things;prey:gets eaten"),
            False,
            True,
            ((INFO, "Checking field iact_field"),),
        ),
        (
            dict(
                interaction_field="predator:eats things", interaction_name="C_born:prey"
            ),
            True,
            True,
            ((INFO, "Checking field iact_field"),),
        ),
        (
            dict(
                interaction_field="predator:eats things;prey:gets eaten",
                interaction_name="C_born:decomposer",
            ),  # Tritrophic example
            True,
            True,
            ((INFO, "Checking field iact_field"),),
        ),
        (
            dict(interaction_name="C_born;V_salv"),
            True,
            False,
            (
                (INFO, "Checking field iact_field"),
                (
                    WARNING,
                    "Label descriptions for interacting taxa incomplete or missing",
                ),
            ),
        ),
        (
            dict(interaction_field="predator;prey"),
            False,
            True,
            (
                (INFO, "Checking field iact_field"),
                (
                    WARNING,
                    "Label descriptions for interacting taxa incomplete or missing",
                ),
            ),
        ),
        (
            dict(interaction_field="predator", interaction_name="C_born"),
            True,
            True,
            (
                (INFO, "Checking field iact_field"),
                (
                    WARNING,
                    "Label descriptions for interacting taxa incomplete or missing",
                ),
            ),
        ),
    ],
)
def test_check_interaction_meta(
    caplog, fixture_dataset, iact_meta, has_taxa_object, has_dwsh_object, expected_log
):
    """Testing the use of the BaseField._check_interaction_meta() method."""

    # Set up what information is available for taxon field validation
    ds_obj = fixture_dataset if has_taxa_object else None

    dwsh = DataWorksheet(
        {
            "name": "DF",
            "title": "My data table",
            "description": "This is a test data worksheet",
        }
    )
    dwsh.taxa_fields = ["predator", "prey"]
    dwsh_obj = dwsh if has_dwsh_object else None

    # Technically, this violates the field_name last requirement, but that is
    # enforced at the dataworksheet level, not the field level.
    field_meta = OrderedDict(
        field_type="numeric", description="description", field_name="iact_field"
    )
    field_meta.update(iact_meta)

    fld = BaseField(field_meta, dwsh=dwsh_obj, dataset=ds_obj)

    caplog.clear()

    # Test the logging from this private method.
    fld._check_interaction_meta()
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


# BaseField behaviour


@pytest.mark.parametrize(
    "field_meta, expected_log",
    [
        (
            {
                "field_type": "location",
                "description": "SAFE 2nd order point number",
                "field_name": "Location",
                "col_idx": 1,
            },
            ((INFO, "Checking field Location"),),
        ),
        (
            {
                "description": "SAFE 2nd order point number",
                "field_name": "Location",
                "col_idx": 1,
            },
            (
                (INFO, "Checking field Location"),
                (ERROR, "field_type descriptor missing"),
            ),
        ),
        (
            {
                "field_type": None,
                "description": "SAFE 2nd order point number",
                "field_name": "Location",
                "col_idx": 1,
            },
            (
                (INFO, "Checking field Location"),
                (ERROR, "field_type descriptor is blank"),
            ),
        ),
        (
            {
                "field_type": "      ",
                "description": "SAFE 2nd order point number",
                "field_name": "Location",
                "col_idx": 1,
            },
            (
                (INFO, "Checking field Location"),
                (ERROR, "field_type descriptor is blank"),
            ),
        ),
        (
            {
                "field_type": 123,
                "description": "SAFE 2nd order point number",
                "field_name": "Location",
                "col_idx": 1,
            },
            (
                (INFO, "Checking field Location"),
                (ERROR, "field_type descriptor is not a string"),
            ),
        ),
        (
            {
                "field_type": "  padded   ",
                "description": "SAFE 2nd order point number",
                "field_name": "Location",
                "col_idx": 1,
            },
            (
                (INFO, "Checking field Location"),
                (ERROR, "field_type descriptor has whitespace padding"),
            ),
        ),
        (
            {
                "field_type": "location",
                "description": "SAFE 2nd order point number",
                "field_name": "123_not a valid name",
                "col_idx": 1,
            },
            (
                (INFO, "Checking field 123_not a valid name"),
                (ERROR, "Field name is not valid"),
            ),
        ),
        (
            {
                "field_type": "location",
                "description": "SAFE 2nd order point number",
                "field_name": 123,
                "col_idx": 1,
            },
            (
                (INFO, "Checking field 123"),
                (ERROR, "field_name descriptor is not a string: 123"),
            ),
        ),
        (
            {
                "field_type": "location",
                "description": "SAFE 2nd order point number",
                "field_name": None,
                "col_idx": 12,
            },
            (
                (INFO, "Checking field Column_L"),
                (ERROR, "field_name descriptor is blank"),
            ),
        ),
    ],
)
def test_BaseField_init(caplog, field_meta, expected_log):
    """Testing behaviour of the BaseField class in handling bad descriptors.

    This is done via __init__ and testing within the _check_meta() method.
    """

    fld = BaseField(field_meta, None)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data, expected_log",
    [
        ([1, 2, 3, 4, 5, 6, 7, 8, 9], ((INFO, "Checking field Location"),)),
        (
            [1, "NA", 3, 4, 5, 6, "NA", 8, 9],
            (
                (INFO, "Checking field Location"),
                (WARNING, "2 / 9 values missing"),
            ),
        ),
        (
            [1, None, 3, 4, 5, 6, "   ", 8, 9],
            (
                (INFO, "Checking field Location"),
                (ERROR, "2 cells are blank or contain only whitespace text"),
            ),
        ),
        (
            [1, 2, 3, 4, "#REF!", 6, "#N/A", 8, 9],
            (
                (INFO, "Checking field Location"),
                (ERROR, "2 cells contain Excel formula errors"),
            ),
        ),
    ],
)
def test_BaseField_run_common_validation(caplog, data, expected_log):
    """Testing behaviour of the BaseField class in using _validate_data."""

    fld = BaseField(
        {
            "field_type": "location",
            "description": "SAFE 2nd order point number",
            "field_name": "Location",
        },
        None,
    )

    fld.run_common_validation(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


# Comments field


@pytest.mark.parametrize(
    "data, expected_log",
    [
        ([1, 2, 3, 4, 5, 6, 7, 8, 9], ((INFO, "Checking field comments"),)),
        (
            ["Any", None, None, "old", 5, "rubbish", None, "is", "fine"],
            ((INFO, "Checking field comments"),),
        ),
    ],
)
def test_CommentField_validate_data(caplog, data, expected_log):
    """Testing behaviour of the CommentsField class in using _validate_data."""

    fld = CommentField(
        {
            "field_type": "comments",
            "description": "my comments",
            "field_name": "comments",
        }
    )

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


# Replicate and ID field


@pytest.mark.parametrize(
    "data, expected_log",
    [
        ([1, 2, 3, 4, 5, 6, 7, 8, 9], ((INFO, "Checking field replicate"),)),
        (
            ["Any", "old", 5, "rubbish", "is", "fine", "including", "NA"],
            ((INFO, "Checking field replicate"), (WARNING, "1 / 8 values missing")),
        ),
        (
            ["Any", "old", 5, "rubbish", "is", "fine", None, "but not missing cells"],
            (
                (INFO, "Checking field replicate"),
                (ERROR, "1 cells are blank or contain only whitespace text"),
            ),
        ),
        (
            [
                "Any",
                "old",
                5,
                "rubbish",
                "is",
                "fine",
                "  ",
                "but not missing or empty string cells",
            ],
            (
                (INFO, "Checking field replicate"),
                (ERROR, "1 cells are blank or contain only whitespace text"),
            ),
        ),
    ],
)
def test_ReplicateField_validate_data(caplog, data, expected_log):
    """Testing behaviour of the CommentsField class in using _validate_data."""

    fld = ReplicateField(
        {
            "field_type": "replicate",
            "description": "my comments",
            "field_name": "replicate",
        }
    )

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


# NumericField and derived classes
# - NumericField itself has BaseField init, just test overloaded validate_data
# - Can reuse the same data to also check taxon and interaction classes
#   inheriting from NumericField for validate_data.
# - The __init__ testing duplicates testing of private methods above but this is
#   so that I can be sure the  class inheritance works as expected
# - Can also reuse the __init__ testing for Numeric and Categorical versions
#   of Taxon/Interaction classes.


@pytest.mark.parametrize(
    "tx_meta, has_taxa_object, has_dwsh_object, expected_log",
    [
        (
            dict(taxon_name="foo", taxon_field="bar"),
            False,
            False,
            (
                (INFO, "Checking field"),
                (ERROR, "Taxon name and taxon field both provided, use one only"),
            ),
        ),
        (
            dict(),
            False,
            False,
            (
                (INFO, "Checking field"),
                (ERROR, "One of taxon name or taxon field must be provided"),
            ),
        ),
        (
            dict(taxon_name="foo"),
            False,
            False,
            (
                (INFO, "Checking field"),
                (CRITICAL, "Taxon name provided but no Taxa instance available"),
            ),
        ),
        (
            dict(taxon_name="foo"),
            True,
            False,
            (
                (INFO, "Checking field"),
                (ERROR, "Taxon name not found in the Taxa worksheet"),
            ),
        ),
        (dict(taxon_name="C_born"), True, False, ((INFO, "Checking field"),)),
        (
            dict(taxon_field="not_provided"),
            False,
            False,
            (
                (INFO, "Checking field"),
                (
                    CRITICAL,
                    "Taxon field provided but no dataworksheet provided for this field",
                ),
            ),
        ),
        (
            dict(taxon_field="not_provided"),
            False,
            True,
            (
                (INFO, "Checking field"),
                (ERROR, "Taxon field not found in this worksheet"),
            ),
        ),
        (dict(taxon_field="my_taxon_field"), False, True, ((INFO, "Checking field"),)),
    ],
)
@pytest.mark.parametrize(
    "test_class, field_meta",
    [
        (
            NumericTaxonField,
            {
                "field_type": "abundance",
                "description": "Number of ants",
                "field_name": "num_data",
                "method": "quadrat",
                "units": "individuals per m2",
            },
        ),
        (
            CategoricalTaxonField,
            {
                "field_type": "categorical trait",
                "description": "ant colours",
                "field_name": "ant_colour",
                "levels": "level1;level2",
            },
        ),
    ],
)
def test_CatNumTaxonField_init(
    caplog,
    fixture_dataset,
    test_class,
    field_meta,
    tx_meta,
    has_taxa_object,
    has_dwsh_object,
    expected_log,
):
    """Testing __init__ behaviour of NumericTaxonField and CategoricalTaxonField."""

    # Set up what information is available for taxon field validation
    ds_obj = fixture_dataset if has_taxa_object else None

    dwsh = DataWorksheet(
        {
            "name": "DF",
            "title": "My data table",
            "description": "This is a test data worksheet",
        }
    )
    dwsh.taxa_fields = ["my_taxon_field"]
    dwsh_obj = dwsh if has_dwsh_object else None

    caplog.clear()

    # Do not directly modify the parameterised value - changes for all tests.
    # (Or that could be poor choice of scope)
    field_meta = field_meta.copy()
    field_meta.update(tx_meta)

    fld = test_class(field_meta, dataset=ds_obj, dwsh=dwsh_obj)

    # For Categorical fields, add some data to avoid unused level errors
    if test_class == CategoricalTaxonField:
        fld.validate_data(["level1", "level2"])

    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "iact_meta, has_taxa_object, has_dwsh_object, expected_log",
    [
        (
            dict(),
            False,
            False,
            (
                (INFO, "Checking field"),
                (
                    ERROR,
                    "At least one of interaction name or interaction field must be"
                    " provided",
                ),
            ),
        ),
        (
            dict(interaction_name="foo:predator"),
            False,
            False,
            (
                (INFO, "Checking field"),
                (CRITICAL, "Interaction name provided but no Taxa instance available"),
            ),
        ),
        (
            dict(interaction_name="foo:predator"),
            True,
            False,
            (
                (INFO, "Checking field"),
                (ERROR, "Unknown taxa in interaction_name descriptor"),
                (
                    ERROR,
                    "At least two interacting taxon labels or fields must be "
                    "identified",
                ),
            ),
        ),
        (
            dict(interaction_name="C_born:predator"),
            True,
            False,
            (
                (INFO, "Checking field"),
                (
                    ERROR,
                    "At least two interacting taxon labels or fields must be "
                    "identified",
                ),
            ),
        ),
        (
            dict(
                interaction_name="C_born:predator;V_salv:prey"
            ),  # Biologically not very plausible
            True,
            False,
            ((INFO, "Checking field"),),
        ),
        (
            dict(interaction_field="not_provided:predator"),
            False,
            False,
            (
                (INFO, "Checking field"),
                (
                    CRITICAL,
                    "Interaction field provided but no dataworksheet provided for this"
                    " field",
                ),
            ),
        ),
        (
            dict(interaction_field="foo:predator"),
            False,
            True,
            (
                (INFO, "Checking field"),
                (ERROR, "Unknown taxon fields in interaction_field descriptor"),
                (
                    ERROR,
                    "At least two interacting taxon labels or fields must be "
                    "identified",
                ),
            ),
        ),
        (
            dict(interaction_field="predator:eats things"),
            False,
            True,
            (
                (INFO, "Checking field"),
                (
                    ERROR,
                    "At least two interacting taxon labels or fields must be "
                    "identified",
                ),
            ),
        ),
        (
            dict(interaction_field="predator:eats things;prey:gets eaten"),
            False,
            True,
            ((INFO, "Checking field"),),
        ),
        (
            dict(
                interaction_field="predator:eats things", interaction_name="C_born:prey"
            ),
            True,
            True,
            ((INFO, "Checking field"),),
        ),
        (
            dict(
                interaction_field="predator:eats things;prey:gets eaten",
                interaction_name="C_born:decomposer",
            ),  # Tritrophic example
            True,
            True,
            ((INFO, "Checking field"),),
        ),
        (
            dict(interaction_name="C_born;V_salv"),
            True,
            False,
            (
                (INFO, "Checking field"),
                (
                    WARNING,
                    "Label descriptions for interacting taxa incomplete or missing",
                ),
            ),
        ),
        (
            dict(interaction_field="predator;prey"),
            False,
            True,
            (
                (INFO, "Checking field"),
                (
                    WARNING,
                    "Label descriptions for interacting taxa incomplete or missing",
                ),
            ),
        ),
        (
            dict(interaction_field="predator", interaction_name="C_born"),
            True,
            True,
            (
                (INFO, "Checking field"),
                (
                    WARNING,
                    "Label descriptions for interacting taxa incomplete or missing",
                ),
            ),
        ),
    ],
)
@pytest.mark.parametrize(
    "test_class, field_meta",
    [
        (
            NumericInteractionField,
            {
                "field_type": "numeric interaction",
                "description": "Number of prey eaten",
                "field_name": "num_data",
                "method": "cage experiment",
                "units": "individuals per hour",
            },
        ),
        (
            CategoricalInteractionField,
            {
                "field_type": "categorical interaction",
                "description": "outcome of competition",
                "field_name": "factor_data",
                "levels": "level1;level2",
            },
        ),
    ],
)
def test_CatNumInteractionField_init(
    caplog,
    fixture_dataset,
    test_class,
    field_meta,
    iact_meta,
    has_taxa_object,
    has_dwsh_object,
    expected_log,
):
    """Testing __init__ behaviour of the Numeric and CategoricalInteractionField."""

    # Set up what information is available for taxon field validation
    ds_obj = fixture_dataset if has_taxa_object else None
    dwsh = DataWorksheet(
        {
            "name": "DF",
            "title": "My data table",
            "description": "This is a test data worksheet",
        }
    )
    dwsh.taxa_fields = ["predator", "prey"]
    dwsh_obj = dwsh if has_dwsh_object else None

    # Do not directly modify the parameterised value - changes for all tests.
    field_meta = field_meta.copy()
    field_meta.update(iact_meta)

    caplog.clear()

    # Test the __init__ method
    fld = test_class(field_meta, dataset=ds_obj, dwsh=dwsh_obj)

    # For Categorical fields, add some data to avoid unused level errors
    if test_class == CategoricalInteractionField:
        fld.validate_data(["level1", "level2"])

    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data, expected_log",
    [
        ([1, 2, 3, 4, 5, 6, 7, 8, 9], ((INFO, "Checking field num_data"),)),
        (
            [1, "NA", 3, 4, 5, 6, "NA", 8, 9],
            ((INFO, "Checking field num_data"), (WARNING, "2 / 9 values missing")),
        ),
        (
            [1, None, 3, 4, 5, 6, "   ", 8, 9],
            (
                (INFO, "Checking field num_data"),
                (ERROR, "2 cells are blank or contain only whitespace text"),
            ),
        ),
        (
            [1, 2, 3, 4, "#REF!", 6, "#N/A", 8, 9],
            (
                (INFO, "Checking field num_data"),
                (ERROR, "2 cells contain Excel formula errors"),
            ),
        ),
        (
            [1, 2, 3, 4, "wrong_type", 6, 7, 8, 9],
            (
                (INFO, "Checking field num_data"),
                (ERROR, "Cells contain non-numeric values"),
            ),
        ),
        (
            [1, 2, "NA", 4, "wrong_type", 6, None, 8, 9],
            (
                (INFO, "Checking field num_data"),
                (ERROR, "Cells contain non-numeric values"),
                (WARNING, "1 / 9 values missing"),
                (ERROR, "1 cells are blank or contain only whitespace text"),
            ),
        ),
    ],
)
@pytest.mark.parametrize(
    "test_class, field_meta",
    [
        (
            NumericField,
            {
                "field_type": "numeric",
                "description": "Tree height",
                "field_name": "num_data",
                "method": "looking",
                "units": "metres",
            },
        ),
        (
            NumericTaxonField,
            {
                "field_type": "abundance",
                "description": "Number of ants",
                "field_name": "num_data",
                "method": "quadrat",
                "units": "individuals per m2",
                "taxon_name": "C_born",
            },
        ),
        (
            NumericInteractionField,
            {
                "field_type": "numeric interaction",
                "description": "Number of prey eaten",
                "field_name": "num_data",
                "method": "cage experiment",
                "units": "individuals per hour",
                "interaction_name": "C_born:prey;V_salv:predator",
            },
        ),
    ],
)
def test_NumericField_and_subclasses_validate_data(
    caplog, fixture_dataset, test_class, field_meta, data, expected_log
):
    """Testing behaviour of the NumericField and subclasses in using _validate_data."""

    # Create and instance of the required class
    fld = test_class(field_meta, dataset=fixture_dataset)

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


# CategoricalField and derived classes
# - Can reuse the same data to also check taxon and interaction classes
#   inheriting from CategoricalField for validate_data.


@pytest.mark.parametrize(
    "field_meta, expected_log",
    [
        (
            {
                "field_type": "categorical",
                "description": "a factor",
                "field_name": "factor1",
                "col_idx": 1,
            },
            (
                (INFO, "Checking field factor1"),
                (ERROR, "levels descriptor missing"),
            ),
        ),
        (
            {
                "field_type": "categorical",
                "description": "a factor",
                "field_name": "factor1",
                "levels": None,
                "col_idx": 1,
            },
            (
                (INFO, "Checking field factor1"),
                (ERROR, "levels descriptor is blank"),
            ),
        ),
        (
            {
                "field_type": "categorical",
                "description": "a factor",
                "field_name": "factor1",
                "levels": "   ",
                "col_idx": 1,
            },
            (
                (INFO, "Checking field factor1"),
                (ERROR, "levels descriptor is blank"),
            ),
        ),
        (
            {
                "field_type": "categorical",
                "description": "a factor",
                "field_name": "factor1",
                "levels": 123,
                "col_idx": 1,
            },
            (
                (INFO, "Checking field factor1"),
                (ERROR, "levels descriptor is not a string"),
            ),
        ),
        (
            {
                "field_type": "categorical",
                "description": "a factor",
                "field_name": "factor1",
                "levels": " level1;level2; ",
                "col_idx": 1,
            },
            (
                (INFO, "Checking field factor1"),
                (ERROR, "levels descriptor has whitespace padding"),
                (ERROR, "Categories found in levels descriptor not used in data:"),
            ),
        ),
    ],
)
def test_CategoricalField_init(caplog, field_meta, expected_log):
    """Testing behaviour of the BaseField class in handling bad descriptors.

    This is done via __init__ and testing within the _check_meta() method.
    """

    fld = CategoricalField(field_meta, None)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data, expected_log",
    [
        (
            ["level1", "level2", "level1", "level2", "level1", "level2"],
            ((INFO, "Checking field factor_data"),),
        ),
        (
            ["level1", "NA", "level1", "level2", "NA", "level2"],
            ((INFO, "Checking field factor_data"), (WARNING, "2 / 6 values missing")),
        ),
        (
            ["level1", "    ", "level1", "level2", None, "level2"],
            (
                (INFO, "Checking field factor_data"),
                (ERROR, "2 cells are blank or contain only whitespace text"),
            ),
        ),
        (
            ["level1", "#REF!", "level1", "#N/A", "level1", "level2"],
            (
                (INFO, "Checking field factor_data"),
                (ERROR, "2 cells contain Excel formula errors"),
            ),
        ),
        (
            ["level1", "level2", "level1", 1234, "level1", "level2"],
            (
                (INFO, "Checking field factor_data"),
                (ERROR, "Cells contain non-text values"),
            ),
        ),
        (
            ["level1", "level2", "NA", None, 1234, "level2"],
            (
                (INFO, "Checking field factor_data"),
                (ERROR, "Cells contain non-text values"),
                (WARNING, "1 / 6 values missing"),
                (ERROR, "1 cells are blank or contain only whitespace text"),
            ),
        ),
        (
            ["level1", "level2", "level3", "level2", "level1", "level2"],
            (
                (INFO, "Checking field factor_data"),
                (ERROR, "Categories found in data missing from levels descriptor"),
            ),
        ),
        (
            ["level1", "level1", "level1", "level1", "level1", "level1"],
            (
                (INFO, "Checking field factor_data"),
                (ERROR, "Categories found in levels descriptor not used in data"),
            ),
        ),
    ],
)
@pytest.mark.parametrize(
    "test_class, field_meta",
    [
        (
            CategoricalField,
            {
                "field_type": "categorical",
                "description": "a factor",
                "field_name": "factor_data",
                "levels": "level1;level2;",
            },
        ),
        (
            CategoricalTaxonField,
            {
                "field_type": "categorical trait",
                "description": "ant colours",
                "field_name": "factor_data",
                "levels": "level1;level2",
                "taxon_name": "C_born",
            },
        ),
        (
            CategoricalInteractionField,
            {
                "field_type": "categorical interaction",
                "description": "outcome of competition",
                "field_name": "factor_data",
                "levels": "level1;level2",
                "interaction_name": "C_born:competitor 1;V_salv:competitor 2",
            },
        ),
    ],
)
def test_CategoricalField_and_subclasses_validate_data(
    caplog, fixture_dataset, test_class, field_meta, data, expected_log
):
    """Testing behaviour of the CategoricalField class in using validate_data."""

    fld = test_class(field_meta, dataset=fixture_dataset)

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


# Taxon field - check init and validate_data (overloaded report just emits messages)


@pytest.mark.parametrize(
    "provide_taxa_instance, expected_log",
    [
        (True, ((INFO, "Checking field taxa_field"), (ERROR, "No taxa loaded"))),
        (
            False,
            (
                (INFO, "Checking field taxa_field"),
                (ERROR, "No taxon details provided for dataset"),
                (ERROR, "No taxa loaded"),
            ),
        ),
    ],
)
def test_TaxaField_init(caplog, fixture_dataset, provide_taxa_instance, expected_log):
    """Testing behaviour of the TaxaField class in handling missing taxa."""

    if provide_taxa_instance:
        ds = fixture_dataset
    else:
        ds = None

    fld = TaxaField(
        {"field_type": "taxa", "description": "My taxa", "field_name": "taxa_field"},
        dataset=ds,
    )
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data, expected_log",
    [
        (
            ["C_born", "V_salv", "C_born", "V_salv", "C_born", "V_salv"],
            ((INFO, "Checking field taxa_field"),),
        ),
        (
            ["C_born", 123, "C_born", "V_salv", "C_born", "V_salv"],
            (
                (INFO, "Checking field taxa_field"),
                (ERROR, "Cells contain non-string values"),
            ),
        ),
        ([], ((INFO, "Checking field taxa_field"), (ERROR, "No taxa loaded"))),
        (
            ["C_born", "V_salv", "C_born", "V_salv", "C_born", "V_salv", "P_leo"],
            ((INFO, "Checking field taxa_field"), (ERROR, "Includes unreported taxa")),
        ),
    ],
)
def test_TaxaField_validate_data(caplog, fixture_dataset, data, expected_log):
    """Testing behaviour of the TaxaField class in using validate_data."""

    fld = TaxaField(
        {"field_type": "taxa", "description": "My taxa", "field_name": "taxa_field"},
        dataset=fixture_dataset,
    )

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


# Locations field - check init and validate_data (overloaded report just emits messages)


@pytest.mark.parametrize(
    "provide_loc_instance, expected_log",
    [
        (True, ((INFO, "Checking field locations"), (ERROR, "No locations loaded"))),
        (
            False,
            (
                (INFO, "Checking field locations"),
                (ERROR, "No location details provided for dataset"),
                (ERROR, "No locations loaded"),
            ),
        ),
    ],
)
def test_LocationsField_init(
    caplog, fixture_dataset, provide_loc_instance, expected_log
):
    """Testing behaviour of the LocationsField class in handling missing locations."""

    if provide_loc_instance:
        ds = fixture_dataset
    else:
        ds = None

    fld = LocationsField(
        {
            "field_name": "locations",
            "field_type": "locations",
            "description": "my locations",
        },
        dataset=ds,
    )
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data, expected_log",
    [  # Good data
        (
            ["A_1", "A_2", "A_1", "A_2", "A_1", "A_2"],
            ((INFO, "Checking field locations"),),
        ),
        (["A_1", "A_2", 1, "A_2", "A_1", 2], ((INFO, "Checking field locations"),)),
        # Bad data
        (
            ["A_1", "A_2", 16.2, "A_2", "A_1", datetime.now()],
            (
                (INFO, "Checking field locations"),
                (ERROR, "Cells contain invalid location values"),
            ),
        ),
        ([], ((INFO, "Checking field locations"), (ERROR, "No locations loaded"))),
        (
            ["A_1", "A_2", "A_3", "A_2", "A_1", "A_3"],
            (
                (INFO, "Checking field locations"),
                (ERROR, "Includes unreported locations"),
            ),
        ),
    ],
)
def test_LocationsField_validate_data(caplog, fixture_dataset, data, expected_log):
    """Testing behaviour of the TaxaField class in using validate_data."""

    fld = LocationsField(
        {
            "field_name": "locations",
            "field_type": "locations",
            "description": "my locations",
        },
        dataset=fixture_dataset,
    )

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


# Geographic coordinates field.


@pytest.mark.parametrize(
    "provide_ds_instance, expected_log",
    [
        (True, ((INFO, "Checking field geocoords"),)),
        (
            False,
            (
                (INFO, "Checking field geocoords"),
                (ERROR, "No dataset object provided - cannot update extents"),
            ),
        ),
    ],
)
def test_GeoField_init(caplog, fixture_dataset, provide_ds_instance, expected_log):
    """Testing behaviour of the GeoField class in handling missing dataset."""

    if provide_ds_instance:
        ds = fixture_dataset
    else:
        ds = None

    fld = GeoField(
        {"field_name": "geocoords", "field_type": "latitude", "description": "my gcs"},
        dataset=ds,
    )
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "which,data, expected_log",
    [  # Good data
        ("latitude", [1, 2, 3, 4, 5, 6], ((INFO, "Checking field geocoords"),)),
        (
            "longitude",
            [111, 112, 113, 114, 115, 116],
            ((INFO, "Checking field geocoords"),),
        ),
        # Bad inputs
        (
            "latitude",
            [1, "2", 3, 4, 5, 6],
            (
                (INFO, "Checking field geocoords"),
                (ERROR, "Field contains non-numeric data"),
            ),
        ),
        (
            "longitude",
            [111, "112", 113, 114, 115, 116],
            (
                (INFO, "Checking field geocoords"),
                (ERROR, "Field contains non-numeric data"),
            ),
        ),
        (
            "latitude",
            [1, "2°", 3, 4, 5, 6],
            (
                (INFO, "Checking field geocoords"),
                (ERROR, "Field contains non-numeric data"),
                (WARNING, "Possible degrees minutes and seconds formatting?"),
            ),
        ),
        (
            "longitude",
            [111, "112°", 113, 114, 115, 116],
            (
                (INFO, "Checking field geocoords"),
                (ERROR, "Field contains non-numeric data"),
                (WARNING, "Possible degrees minutes and seconds formatting?"),
            ),
        ),
        (
            "latitude",
            [1, 2, 3, 4, 5, 600],
            ((INFO, "Checking field geocoords"), (ERROR, "exceeds hard bounds")),
        ),
        (
            "longitude",
            [-200, 2, 3, 4, 5, 600],
            ((INFO, "Checking field geocoords"), (ERROR, "exceeds hard bounds")),
        ),
        (
            "latitude",
            [1, 2, 3, 4, 5, 60],
            ((INFO, "Checking field geocoords"), (WARNING, "exceeds soft bounds")),
        ),
        (
            "longitude",
            [111, 112, 113, 114, 115, 160],
            ((INFO, "Checking field geocoords"), (WARNING, "exceeds soft bounds")),
        ),
    ],
)
def test_GeoField_validate_data(caplog, fixture_dataset, which, data, expected_log):
    """Testing behaviour of the TaxaField class in using validate_data."""

    fld = GeoField(
        {"field_name": "geocoords", "field_type": which, "description": "my gcs"},
        dataset=fixture_dataset,
    )

    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data, expected_log",
    [  # Good data
        (
            [
                [time(11, 12), time(11, 12), time(11, 12), time(11, 12), time(11, 12)],
            ],
            ((INFO, "Checking field time"),),
        ),
        (
            [
                ["11:12:13", "11:12:13", "11:12:13", "11:12:13", "11:12:13"],
            ],
            ((INFO, "Checking field time"),),
        ),
        (
            [
                ["111213", "1112", "111213.2", "11:12:13", "11:12:13"],
            ],  # Different ISO8601 string formats
            ((INFO, "Checking field time"),),
        ),
        # Bad inputs
        (
            [
                [time(11, 12), time(11, 12), "11:12:13", "11:12:13", "11:12:13"],
            ],
            (
                (INFO, "Checking field time"),
                (ERROR, "Time data mixes ISO string and time formatted rows"),
            ),
        ),
        (
            [
                [time(11, 12), time(11, 12), 1, 2, 3],
            ],
            (
                (INFO, "Checking field time"),
                (
                    ERROR,
                    "Time data include values neither ISO string nor time formatted",
                ),
            ),
        ),
        (
            [
                ["11twelve13", "11:12:13", "11/12/13", "11:12:13", "11:12:13"],
            ],
            (
                (INFO, "Checking field time"),
                (ERROR, "ISO time strings contain badly formatted values"),
            ),
        ),
        (
            [
                [time(11, 12), time(11, 12), time(11, 12), time(11, 12), time(11, 12)],
                ["11:12:13", "11:12:13", "11:12:13", "11:12:13", "11:12:13"],
            ],
            (
                (INFO, "Checking field time"),
                (ERROR, "Time data mixes ISO string and time formatted rows"),
            ),
        ),
        (
            [
                [time(11, 12), time(11, 12), time(11, 12), time(11, 12), time(11, 12)],
                [1, 2, "11:12:13", "11:12:13", "11:12:13"],
            ],
            (
                (INFO, "Checking field time"),
                (
                    ERROR,
                    "Time data include values neither ISO string nor time formatted",
                ),
            ),
        ),
    ],
)
def test_TimeField_validate_data(caplog, data, expected_log):
    """Testing behaviour of the TaxaField class in using validate_data."""

    fld = TimeField({"field_name": "time", "field_type": "time", "description": "time"})

    for d in data:
        fld.validate_data(d)

    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "field_type, data, expected_log",
    [  # Good data
        (
            "date",
            [
                [
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                ],
            ],
            ((INFO, "Checking field datetimetest"),),
        ),
        (
            "date",
            [
                [
                    "2022-01-06",
                    "2022-01-06",
                    "2022-01-06",
                    "2022-01-06",
                    "2022-01-06",
                ],
            ],
            ((INFO, "Checking field datetimetest"),),
        ),
        (
            "datetime",
            [
                [
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                ],
            ],
            ((INFO, "Checking field datetimetest"),),
        ),
        (
            "datetime",
            [
                [
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                ],
            ],
            ((INFO, "Checking field datetimetest"),),
        ),
        # Bad data - mixed inputs
        (
            "datetime",
            [
                [
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (WARNING, "Field is of type datetime, but only reports dates"),
            ),
        ),
        (
            "date",
            [
                [
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (ERROR, "Field is of type date, but includes time data"),
            ),
        ),
        (
            "datetime",
            [
                [
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (
                    ERROR,
                    "Date or datetime data mixes ISO string and time formatted rows",
                ),
            ),
        ),
        (
            "date",
            [
                [
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    "2022-01-06",
                    "2022-01-06",
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (
                    ERROR,
                    "Date or datetime data mixes ISO string and time formatted rows",
                ),
            ),
        ),
        (
            "date",
            [
                [
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                ],
                ["2022-01-06", "2022-01-06", "2022-01-06", "2022-01-06", "2022-01-06"],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (
                    ERROR,
                    "Date or datetime data mixes ISO string and time formatted rows",
                ),
            ),
        ),
        (
            "datetime",
            [
                [
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                ],
                [
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (
                    ERROR,
                    "Date or datetime data mixes ISO string and time formatted rows",
                ),
            ),
        ),
        # Bad data - bad classes
        (
            "date",
            [
                [
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    1,
                    2,
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (
                    ERROR,
                    "Date or datetime data include values neither ISO string nor time"
                    " formatted",
                ),
            ),
        ),
        (
            "date",
            [
                ["2022-01-06", "2022-01-06", "2022-01-06", 1, 2],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (
                    ERROR,
                    "Date or datetime data include values neither ISO string nor time"
                    " formatted",
                ),
            ),
        ),
        (
            "datetime",
            [
                [
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    1,
                    2,
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (
                    ERROR,
                    "Date or datetime data include values neither ISO string nor time"
                    " formatted",
                ),
            ),
        ),
        (
            "datetime",
            [
                ["2022-01-06 11:12", "2022-01-06 11:12", "2022-01-06 11:12", 1, 2],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (
                    ERROR,
                    "Date or datetime data include values neither ISO string nor time"
                    " formatted",
                ),
            ),
        ),
        # Bad data - string formats
        (
            "date",
            [
                [
                    "2022-01-06",
                    "2022/01/06",
                    "2022-01-06",
                    "20220106",
                    "2022-01-06",
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (ERROR, "ISO datetime strings contain badly formatted values"),
            ),
        ),
        (
            "datetime",
            [
                [
                    "2022-01-06 11:12",
                    "2022/01/06 11:12",
                    "2022-01-06 11:12",
                    "20220106 11:12",
                    "2022-01-06 11:12",
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (ERROR, "ISO datetime strings contain badly formatted values"),
            ),
        ),
    ],
)
def test_DatetimeField_validate_data(
    caplog, fixture_dataset, field_type, data, expected_log
):
    """Testing behaviour of the DatetimeField class in using validate_data."""

    fld = DatetimeField(
        {"field_name": "datetimetest", "field_type": field_type, "description": "time"},
        dataset=fixture_dataset,
    )

    for d in data:
        fld.validate_data(d)

    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "provide_dataset,field_type, data, expected_log",
    [  # Good data - no dataset
        (
            False,
            "date",
            [
                [
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (ERROR, "No dataset object provided - cannot update extents"),
            ),
        ),
        (
            False,
            "date",
            [
                [
                    "2022-01-06",
                    "2022-01-06",
                    "2022-01-06",
                    "2022-01-06",
                    "2022-01-06",
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (ERROR, "No dataset object provided - cannot update extents"),
            ),
        ),
        (
            False,
            "datetime",
            [
                [
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (ERROR, "No dataset object provided - cannot update extents"),
            ),
        ),
        (
            False,
            "datetime",
            [
                [
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                ],
            ],
            (
                (INFO, "Checking field datetimetest"),
                (ERROR, "No dataset object provided - cannot update extents"),
            ),
        ),
        # Good data - with dataset
        (
            True,
            "date",
            [
                [
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                    datetime(2022, 1, 6),
                ],
            ],
            ((INFO, "Checking field datetimetest"),),
        ),
        (
            True,
            "date",
            [
                [
                    "2022-01-06",
                    "2022-01-06",
                    "2022-01-06",
                    "2022-01-06",
                    "2022-01-06",
                ],
            ],
            ((INFO, "Checking field datetimetest"),),
        ),
        (
            True,
            "datetime",
            [
                [
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                    datetime(2022, 1, 6, 11, 12),
                ],
            ],
            ((INFO, "Checking field datetimetest"),),
        ),
        (
            True,
            "datetime",
            [
                [
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                    "2022-01-06 11:12",
                ],
            ],
            ((INFO, "Checking field datetimetest"),),
        ),
    ],
)
def test_DatetimeField_extent(
    caplog, fixture_dataset, provide_dataset, field_type, data, expected_log
):
    """Testing behaviour of the DatetimeField class in setting extents."""

    if provide_dataset:
        ds = fixture_dataset
    else:
        ds = None

    fld = DatetimeField(
        {"field_name": "datetimetest", "field_type": field_type, "description": "time"},
        dataset=ds,
    )

    for d in data:
        fld.validate_data(d)

    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )

    if provide_dataset and field_type == "date":
        assert fixture_dataset.temporal_extent.extent[0] == date(2022, 1, 6)
        assert fixture_dataset.temporal_extent.extent[1] == date(2022, 1, 6)
    if provide_dataset and field_type == "datetime":
        assert fixture_dataset.temporal_extent.extent[0] == date(2022, 1, 6)
        assert fixture_dataset.temporal_extent.extent[1] == date(2022, 1, 6)


@pytest.mark.parametrize(
    "provide_dataset, use_file_container, external_files, expected_log",
    [
        (
            False,
            False,
            [],
            (
                (INFO, "Checking field file"),
                (CRITICAL, "No Summary instance provided - cannot check file fields"),
            ),
        ),
        (
            True,
            False,
            [],
            (
                (INFO, "Checking field file"),
                (ERROR, "No external files listed in Summary"),
            ),
        ),
        (
            False,
            True,  # Never gets to filename checking
            [{"file": "extfile1.sql"}, {"file": "extfile2.zip"}],
            (
                (INFO, "Checking field file"),
                (CRITICAL, "No Summary instance provided - cannot check file fields"),
            ),
        ),
        (
            True,
            True,
            [{"file": "extfile4.sql"}, {"file": "extfile6.zip"}],
            (
                (INFO, "Checking field file"),
                (ERROR, "Field file_container value not found in external files"),
            ),
        ),
        (
            True,
            True,
            [{"file": "extfile1.sql"}, {"file": "extfile2.zip"}],
            ((INFO, "Checking field file"),),
        ),
    ],
)
def test_FileField_init(
    caplog,
    fixture_dataset,
    provide_dataset,
    use_file_container,
    external_files,
    expected_log,
):
    """Testing behaviour of the FileField class in using init."""

    if provide_dataset:
        ds = fixture_dataset
    else:
        ds = None

    field_meta = {"field_name": "file", "field_type": "file", "description": "file"}

    if use_file_container:
        field_meta["file_container"] = "extfile1.sql"

    if provide_dataset and external_files:
        ds.summary.external_files = external_files

    fld = FileField(field_meta, dataset=ds)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data, expected_log",
    [
        (
            ["extfile1.sql", "extfile2.zip", "extfile1.sql", "extfile2.zip"],
            ((INFO, "Checking field file"),),
        ),
        (
            ["extfile1.sql", "extfile2.zip", "extfile6.sql", "extfile2.zip"],
            (
                (INFO, "Checking field file"),
                (ERROR, "Field contains external files not provided in Summary"),
            ),
        ),
    ],
)
def test_FileField_validate_data(caplog, fixture_dataset, data, expected_log):
    """Testing behaviour of the FileField class in using validate_data."""

    field_meta = {"field_name": "file", "field_type": "file", "description": "file"}

    fixture_dataset.summary.external_files = [
        {"file": "extfile1.sql"},
        {"file": "extfile2.zip"},
    ]

    fld = FileField(field_meta, dataset=fixture_dataset)
    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "provide_dataset, use_file_container, external_files, expected_log",
    [
        (
            False,
            False,
            [],
            (
                (INFO, "Checking field file"),
                (CRITICAL, "No Summary instance provided - cannot check file fields"),
            ),
        ),
        (
            True,
            False,
            [],
            (
                (INFO, "Checking field file"),
                (ERROR, "No external files listed in Summary"),
            ),
        ),
        (
            False,
            True,  # Never gets to filename checking
            [{"file": "extfile1.sql"}, {"file": "extfile2.zip"}],
            (
                (INFO, "Checking field file"),
                (CRITICAL, "No Summary instance provided - cannot check file fields"),
            ),
        ),
        (
            True,
            True,
            [{"file": "extfile4.sql"}, {"file": "extfile6.zip"}],
            (
                (INFO, "Checking field file"),
                (ERROR, "Field file_container value not found in external files"),
            ),
        ),
        (
            True,
            True,
            [{"file": "extfile1.sql"}, {"file": "extfile2.zip"}],
            ((INFO, "Checking field file"),),
        ),
    ],
)
def test_FileField_init_v2(
    caplog,
    fixture_dataset,
    provide_dataset,
    use_file_container,
    external_files,
    expected_log,
):
    """Testing behaviour of the FileField class in using init."""

    if provide_dataset:
        ds = fixture_dataset
    else:
        ds = None

    field_meta = {"field_name": "file", "field_type": "file", "description": "file"}

    if use_file_container:
        field_meta["file_container"] = "extfile1.sql"

    if provide_dataset and external_files:
        ds.summary.external_files = external_files

    fld = FileField(field_meta, dataset=ds)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )


@pytest.mark.parametrize(
    "data, expected_log",
    [
        ([None, None, None, None, None, None, None, None, None], tuple()),
        ([None, None, None, None, "    ", None, "   \t", None, None], tuple()),
        (
            [None, None, 1, None, "    ", None, "   \t", None, None],
            (
                (INFO, "Checking field Unknown"),
                (ERROR, "Trailing field with no descriptors contains data."),
            ),
        ),
    ],
)
def test_EmptyField_validate_data(caplog, data, expected_log):
    """Testing behaviour of the EmptyField class in using validate_data."""

    field_meta = {"field_name": None, "field_type": None, "description": None}

    fld = EmptyField(field_meta)
    fld.validate_data(data)
    fld.report()

    assert len(expected_log) == len(caplog.records)

    assert all(
        [exp[0] == rec.levelno for exp, rec in zip(expected_log, caplog.records)]
    )
    assert all(
        [exp[1] in rec.message for exp, rec in zip(expected_log, caplog.records)]
    )

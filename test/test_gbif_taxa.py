"""Tests checking that the GBIF specific classes work as intended."""

from logging import ERROR, INFO, WARNING

import pytest
from dotmap import DotMap

from safedata_validator import taxa

from .conftest import log_check

# ------------------------------------------
# Testing GBIFTaxon
# ------------------------------------------


@pytest.mark.parametrize(
    "test_input,expected_exception",
    [
        (dict(), TypeError),  # no parameters
        (dict(name=1, rank=1), TypeError),  # non string name and rank
        (
            dict(
                name="Bombus bombus", rank="species", gbif_id="error"
            ),  # non-numeric gbif_id
            TypeError,
        ),
        (
            dict(
                name="Bombus bombus", rank="species", gbif_id=3.141
            ),  # non-integer gbif_id
            ValueError,
        ),
    ],
)
def test_taxon_init_errors(test_input, expected_exception):
    """This test checks taxon inputs throw errors as expected."""
    with pytest.raises(expected_exception):
        _ = taxa.GBIFTaxon(**test_input)


# ------------------------------------------
# Testing taxon validators
# ------------------------------------------


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            dict(name="Crematogaster borneensis", rank="species"),
            ("found", True, 1324716, None),
        ),
        (
            dict(name="Morus", rank="genus", gbif_id=2480962),
            ("found", True, 2480962, None),
        ),
        (
            dict(name="Alsomitra simplex", rank="species"),
            ("found", False, 5537041, 3623287),
        ),
    ],
)
def test_validator_search(fixture_gbif_validator, test_input, expected):
    """This test checks inputs against expected outputs for both validator classes.

    Tests the behaviour for deleted taxa as well - these should not pass checks, even
    when an explicit ID is provided to the deleted taxon.
    """

    tx = taxa.GBIFTaxon(**test_input)
    srch_out = fixture_gbif_validator.search(tx)

    assert srch_out.lookup_status == expected[0]
    assert srch_out.is_canon == expected[1]
    assert srch_out.gbif_id == expected[2]

    if srch_out.lookup_status == "found" and not srch_out.is_canon:
        assert srch_out.canon_usage.gbif_id == expected[3]


@pytest.mark.parametrize(
    argnames=["test_input", "expected"],
    argvalues=[
        (
            dict(name="Acladium atrum", rank="species"),
            ("Deleted taxon", False, 3367489, None),
        ),
        (
            dict(name="Acladium atrum", rank="species", gbif_id=3367489),
            ("Deleted taxon", False, 3367489, None),
        ),
    ],
)
def test_validator_search_deleted(fixture_gbif_validator, test_input, expected):
    """This test checks behaviour for searching for deleted taxa."""

    tx = taxa.GBIFTaxon(**test_input)
    srch_out = fixture_gbif_validator.search(tx)

    assert srch_out.lookup_status == expected[0]
    assert srch_out.is_canon == expected[1]
    assert srch_out.gbif_id == expected[2]

    if srch_out.lookup_status == "found" and not srch_out.is_canon:
        assert srch_out.canon_usage.gbif_id == expected[3]


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (1324716, ("found", "Crematogaster borneensis")),
        (2480962, ("found", "Morus")),
        (3367489, ("Deleted taxon", "Acladium atrum")),
    ],
)
def test_validator_gbif_lookup_outputs(fixture_gbif_validator, test_input, expected):
    """This test checks inputs against expected outputs for both validator classes."""

    found = fixture_gbif_validator.id_lookup(test_input)

    assert found.lookup_status == expected[0]
    assert found.name == expected[1]


@pytest.mark.parametrize(
    "test_input,expected_exception",
    [
        (None, ValueError),  # no parameters
        ("invalid_string", ValueError),  # a string
        (-1, ValueError),  # bad ID
        (100000000000000, taxa.GBIFError),  # bad ID
    ],
)
def test_validator_gbif_lookup_errors(
    fixture_gbif_validator, test_input, expected_exception
):
    """This test checks validator.id_lookup inputs throw errors as expected."""

    with pytest.raises(expected_exception):
        _ = fixture_gbif_validator.id_lookup(test_input)


# ------------------------------------------
# Testing GBIFTaxa validate_and_add_taxon
# - This  method tests the individual tuples of taxon and parent
#   provided in the GBIFTaxa sheet, separating loading from testing. The test
#   cases test the various combinations of good and bad inputs and checks
#   that the set of emitted logging records is as expected.
# ------------------------------------------


@pytest.mark.parametrize(
    "taxon_tuple,expected_log_entries",
    [  # BAD INPUTS: Testing checking of sanitisation
        # Missing workname
        (
            (
                None,
                ["mspecies", "Morphospecies", None, None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "Worksheet name missing, whitespace only or not text"),
                (INFO, "accepted"),
                (INFO, "has valid parent information"),
            ),
        ),
        # Non string workname
        (
            (
                123,
                ["mspecies", "Morphospecies", None, None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "Worksheet name missing, whitespace only or not text"),
                (INFO, "accepted"),
                (INFO, "has valid parent information"),
            ),
        ),
        # Whitespace workname
        (
            (
                "    \n",
                ["mspecies", "Morphospecies", None, None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "Worksheet name missing, whitespace only or not text"),
                (INFO, "accepted"),
                (INFO, "has valid parent information"),
            ),
        ),
        # Padded workname
        (
            (
                " Morphospecies",
                ["mspecies", "Morphospecies", None, None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "Worksheet name has whitespace padding"),
                (INFO, "accepted"),
                (INFO, "has valid parent information"),
            ),
        ),
        # Missing taxon name
        (
            (
                "mspecies 1",
                [None, "Morphospecies", None, None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "Taxon name missing, whitespace only or not text"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Non string taxon name
        (
            (
                "mspecies 1",
                [123, "Morphospecies", None, None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "Taxon name missing, whitespace only or not text"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Whitespace taxon name
        (
            (
                "mspecies 1",
                ["    \n", "Morphospecies", None, None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "Taxon name missing, whitespace only or not text"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Non-numeric GBIF ID
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", "abc", None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "GBIF ID contains value that is not an integer"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Non-numeric GBIF ID
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", 123.4, None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "GBIF ID contains value that is not an integer"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Unknown GBIF ID
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, None],
                ["Formicidae", "Family", 15, None],
            ),
            (
                (ERROR, "GBIF ID problem: GBIF ID not found"),
                (ERROR, "Taxon of type morphospecies has invalid parent information."),
            ),
        ),
        # Non-numeric Ignore ID
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, "abc"],
                ["Formicidae", "Family", None, None],
            ),
            (
                (ERROR, "Ignore ID contains value that is not an integer"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Parent details bad - missing name
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, None],
                [None, "Family", None, None],
            ),
            (
                (ERROR, "Parent name missing or not text"),
                (ERROR, "Parent taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Parent details bad - non-string name
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, None],
                [123, "Family", None, None],
            ),
            (
                (ERROR, "Parent name missing or not text"),
                (ERROR, "Parent taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Parent details bad - padded name
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, None],
                ["Formicidae ", "Family", None, None],
            ),
            (
                (ERROR, "Parent name has whitespace padding"),
                (INFO, "accepted"),
                (INFO, "Taxon of type morphospecies has valid parent information"),
            ),
        ),
        # Parent details bad - missing rank
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, None],
                ["Formicidae", None, None, None],
            ),
            (
                (ERROR, "Parent rank missing or not text"),
                (ERROR, "Parent taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Parent details bad - non-string rank
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, None],
                ["Formicidae", 123, None, None],
            ),
            (
                (ERROR, "Parent rank missing or not text"),
                (ERROR, "Parent taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Parent details bad - padded name
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, None],
                ["Formicidae", "Family ", None, None],
            ),
            (
                (ERROR, "Parent rank has whitespace padding"),
                (INFO, "accepted"),
                (INFO, "Taxon of type morphospecies has valid parent information"),
            ),
        ),
        # Parent details bad - bad gbif id
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, None],
                ["Formicidae", "Family", "abc", None],
            ),
            (
                (ERROR, "Parent GBIF ID contains value that is not an integer"),
                (ERROR, "Parent taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Parent details bad - bad gbif id
        (
            (
                "mspecies 1",
                ["mspecies 1", "Morphospecies", None, None],
                ["Formicidae", "Family", 123.4, None],
            ),
            (
                (ERROR, "Parent GBIF ID contains value that is not an integer"),
                (ERROR, "Parent taxon details not properly formatted, cannot validate"),
            ),
        ),
    ],
)
def test_validate_taxon_sanitise(
    caplog, fixture_resources, taxon_tuple, expected_log_entries
):
    # These tests check that bad inputs are caught correctly

    taxa_instance = taxa.GBIFTaxa(fixture_resources)
    taxa_instance.validate_and_add_taxon(taxon_tuple)

    log_check(caplog, expected_log_entries)


@pytest.mark.parametrize(
    "taxon_tuple,expected_log_entries",
    [  # IGNORED: Testing ignored backbone matches
        # Ignored non backbone (no parent)
        (
            ("Morphospecies", ["mspecies", "Morphospecies", None, 123456789], None),
            (  # (INFO, 'No parent taxon provided'),
                (ERROR, "Ignore ID can only be used with GBIF backbone taxon ranks"),
                (ERROR, "Taxa with Ignore ID must provide parent information"),
            ),
        ),
        # Ignored unknown (no parent)
        (
            ("Unknown species", ["Unknown species", "species", None, 123456789], None),
            (  # (INFO, 'No parent taxon provided'),
                (ERROR, "Taxon with Ignore ID not found in GBIF backbone"),
                (ERROR, "Taxa with Ignore ID must provide parent information"),
            ),
        ),
        # - Ignoring good GBIF match with no parent
        (
            (
                "Crematogaster borneensis",
                ["Crematogaster borneensis", "Species", None, 1324716],
                None,
            ),
            (  # (INFO, 'No parent taxon provided'),
                (INFO, "Canon GBIF usage ignored"),
                (ERROR, "Taxa with Ignore ID must provide parent information"),
            ),
        ),
        # - Ignoring GBIF match with bad parent
        (
            (
                "Crematogaster borneensis",
                ["Crematogaster borneensis", "Species", None, 1324716],
                ["Madeupnotarealfamilyidae", "Family", None, None],
            ),
            (
                (ERROR, "No match found"),
                (INFO, "Canon GBIF usage ignored"),
                (ERROR, "Taxon with Ignore ID has invalid parent information"),
            ),
        ),
        # - Ignoring GBIF match on canon taxa with good parent and correct ID
        (
            (
                "Panthera leo",
                ["Panthera leo", "Species", None, 5219404],
                ["Felis", "Genus", None, None],
            ),
            (
                (INFO, "accepted"),
                (INFO, "Canon GBIF usage ignored"),
                (INFO, "Taxon with ignored canon usage has valid parent information"),
            ),
        ),
        # - Ignoring GBIF match on canon taxa with good parent and bad ID
        (
            (
                "Panthera leo",
                ["Panthera leo", "Species", None, 123456789],
                ["Felis", "Genus", None, None],
            ),
            (
                (INFO, "accepted"),
                (ERROR, "Ignore ID does not match the canon GBIF usage"),
                (INFO, "Taxon with ignored canon usage has valid parent information"),
            ),
        ),
        # - Ignoring GBIF match on synonym with good parent but non canon ID (synonym)
        (
            (
                "Felis imperialis",
                ["Felis imperialis", "Species", None, 8531298],
                ["Felis", "Genus", None, None],
            ),
            (
                (INFO, "accepted"),
                (
                    ERROR,
                    "Taxon is non-canon and Ignore ID does not match the canon GBIF "
                    "usage",
                ),
                (INFO, "Taxon with ignored canon usage has valid parent information"),
            ),
        ),
        # - Ignoring GBIF match on synonym with good parent using correct canon ID
        (
            (
                "Felis imperialis",
                ["Felis imperialis", "Species", None, 5219404],
                ["Felis", "Genus", None, None],
            ),
            (
                (INFO, "accepted"),
                (INFO, "Canon GBIF usage ignored"),
                (INFO, "Taxon with ignored canon usage has valid parent information"),
            ),
        ),
        # NON-BACKBONE Testing non backbone outcomes
        # No parent on non-backbone
        (
            ("Morphospecies", ["mspecies", "Morphospecies", None, None], None),
            (  # (INFO, 'No parent taxon provided'),
                (ERROR, " must provide parent information"),
            ),
        ),
        # Invalid and valid parent cases tested above already but check main taxon
        # logging
        # Unknown parent
        (
            (
                "Morphospecies",
                ["mspecies", "Morphospecies", None, None],
                ["Madeupnotarealfamilyidae", "Family", None, None],
            ),
            (
                (ERROR, "No match found"),
                (ERROR, "has invalid parent information"),
            ),
        ),
        # Normal valid parent
        (
            (
                "Morphospecies",
                ["mspecies", "Morphospecies", None, None],
                ["Formicidae", "Family", None, None],
            ),
            (
                (INFO, "accepted"),
                (INFO, "has valid parent information"),
            ),
        ),
        # BACKBONE: Testing backbone outcomes
        # - Good match - no parent
        (
            (
                "Crematogaster borneensis",
                ["Crematogaster borneensis", "Species", None, None],
                None,
            ),
            (  # (INFO, 'No parent taxon provided'),
                (INFO, "Taxon found in GBIF backbone (accepted)"),
            ),
        ),
        # - Synonymous match - no parent
        (
            ("Felis imperialis", ["Felis imperialis", "Species", None, None], None),
            (  # (INFO, 'No parent taxon provided'),
                (WARNING, "Taxon considered a synonym"),
            ),
        ),
        # - Good match, valid parent
        (
            (
                "Crematogaster borneensis",
                ["Crematogaster borneensis", "Species", None, None],
                ["Crematogaster", "Genus", None, None],
            ),
            (
                (INFO, "accepted"),
                (
                    INFO,
                    "Taxon in GBIF backbone (accepted) with compatible parent "
                    "information",
                ),
            ),
        ),
        # - Good match, incompatible parent
        (
            (
                "Crematogaster borneensis",
                ["Crematogaster borneensis", "Species", None, None],
                ["Pinus", "Genus", None, None],
            ),
            (
                (INFO, "accepted"),
                (
                    ERROR,
                    "Taxon in GBIF backbone (accepted) with incompatible parent "
                    "information",
                ),
            ),
        ),
        # - Good match, incompatible parent
        (
            (
                "Crematogaster borneensis",
                ["Crematogaster borneensis", "Species", None, None],
                ["Madeupnotarealgenus", "Genus", None, None],
            ),
            (
                (ERROR, "No match found"),
                (
                    ERROR,
                    "Taxon in GBIF backbone (accepted) but with invalid parent "
                    "information",
                ),
            ),
        ),
        #  - No match, no parent
        (
            (
                "Crematogaster ormei",
                ["Crematogaster ormei", "Species", None, None],
                None,
            ),
            (  # (INFO, 'No parent taxon provided'),
                (ERROR, "Taxon name and rank combination not found"),
            ),
        ),
        #  - No match, bad parent
        (
            (
                "Crematogaster ormei",
                ["Crematogaster ormei", "Species", None, None],
                ["Madeupnotarealgenus", "Genus", None, None],
            ),
            (
                (ERROR, "No match found"),
                (ERROR, "Taxon not found in GBIF and has invalid parent information"),
            ),
        ),
        # - No match, valid parent
        (
            (
                "Crematogaster ormei",
                ["Crematogaster ormei", "Species", None, None],
                ["Crematogaster", "Genus", None, None],
            ),
            (
                (INFO, "accepted"),
                (INFO, "Taxon not found in GBIF but has valid parent information"),
            ),
        ),
        # ODDITIES
        # - Subspecies
        (
            (
                "Water monitor",
                ["Varanus salvator macromaculatus", "Subspecies", None, None],
                None,
            ),
            (  # (INFO, 'No parent taxon provided'),
                (INFO, "Taxon found in GBIF backbone (accepted)"),
            ),
        ),
        # - Ambiguous
        (
            ("Gannets", ["Morus", "Genus", None, None], None),
            ((ERROR, "GBIF issue:"),),  # (INFO, 'No parent taxon provided'),
        ),
        # - Ambiguous resolved with GBIF ID
        (
            ("Gannets", ["Morus", "Genus", 2480962, None], None),
            (  # (INFO, 'No parent taxon provided'),
                (INFO, "Taxon found in GBIF backbone (accepted)"),
            ),
        ),
        #  - Ambiguous parent
        (
            (
                "Solenopsis #1",
                ["NA", "Morphospecies", None, None],
                ["Solenopsis", "Genus", None, None],
            ),
            (
                (ERROR, "Multiple equal matches"),
                (ERROR, "invalid parent information"),
            ),
        ),
        #  - Ambiguous parent resolved
        (
            (
                "Solenopsis #1",
                ["NA", "Morphospecies", None, None],
                ["Solenopsis", "Genus", 9107211, None],
            ),
            (
                (INFO, "accepted"),
                (INFO, "valid parent information"),
            ),
        ),
    ],
)
def test_validate_taxon_lookup(
    caplog, fixture_resources, taxon_tuple, expected_log_entries
):
    taxa_instance = taxa.GBIFTaxa(fixture_resources)
    taxa_instance.validate_and_add_taxon(taxon_tuple)

    log_check(caplog, expected_log_entries)


@pytest.mark.parametrize(
    "taxon_tuple,expected_log_entries",
    [
        #  - Deleted taxon - need to pick an example taxon here with no accepted taxon
        #    of the same name and rank.
        (
            (
                "Deleted species",
                ["Acladium atrum", "species", None, None],
                None,
            ),
            ((ERROR, "Deleted taxon"),),
        ),
        (
            (
                "Deleted species",
                ["Acladium atrum", "species", 3367489, None],
                None,
            ),
            ((ERROR, "Deleted taxon"),),
        ),
    ],
)
def test_validate_deleted_taxon_lookup(
    caplog, fixture_resources, taxon_tuple, expected_log_entries
):
    taxa_instance = taxa.GBIFTaxa(fixture_resources)
    taxa_instance.validate_and_add_taxon(taxon_tuple)

    log_check(caplog, expected_log_entries)


# TODO - add test to check that duplicate worksheet names are caught.


@pytest.mark.parametrize(
    argnames=["mock_output", "expected_log_entries"],
    argvalues=[
        pytest.param(
            DotMap({"data_columns": []}),
            (
                (INFO, "Loading GBIFTaxa worksheet"),
                (INFO, "Reading taxa data"),
                (ERROR, "No data or only headers in GBIFTaxa worksheet"),
            ),
            id="No data in sheet",
        ),
        pytest.param(
            DotMap({"data_columns": ["some_columns"], "bad_headers": "duplicated"}),
            (
                (INFO, "Loading GBIFTaxa worksheet"),
                (INFO, "Reading taxa data"),
                (ERROR, "Cannot parse taxa with duplicated headers"),
            ),
            id="Duplicated headers",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple()],
                    "headers": ["name", "taxon name"],
                }
            ),
            (
                (INFO, "Loading GBIFTaxa worksheet"),
                (INFO, "Reading taxa data"),
                (ERROR, "Missing core fields:"),
            ),
            id="Missing core field",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple(), tuple(), tuple()],
                    "headers": [
                        "name",
                        "taxon name",
                        "taxon type",
                        "additional header",
                    ],
                }
            ),
            (
                (INFO, "Loading GBIFTaxa worksheet"),
                (INFO, "Reading taxa data"),
                (INFO, "Additional fields provided:"),
                (INFO, "No taxon rows found"),
            ),
            id="Additional header",
        ),
    ],
)
def test_load_worksheet_headers(
    caplog, mocker, fixture_resources, mock_output, expected_log_entries
):
    """Test that unexpected header names are caught by load."""
    from safedata_validator.logger import FORMATTER

    # Create GBIFTaxa class
    tx = taxa.GBIFTaxa(fixture_resources)

    # Setup mocking
    mock_get = mocker.patch("safedata_validator.taxa.GetDataFrame")
    mock_get.return_value = mock_output

    tx.load("meaningless_string")
    # This is needed to ensure that the logging depth is not altered by this test
    FORMATTER.pop()

    log_check(caplog, expected_log_entries)


@pytest.mark.parametrize(
    "example_excel_files, n_errors, n_taxa",
    [("good", 0, 20), ("bad", 20, 20)],
    indirect=["example_excel_files"],  # take actual params from fixture
)
def test_taxa_load(example_excel_files, fixture_resources, n_errors, n_taxa):
    """This tests the ensemble loading of taxa from a file.

    It uses indirect parameterisation to access the fixtures containing the sample excel
    files.
    """

    tx = taxa.GBIFTaxa(fixture_resources)
    tx.load(example_excel_files["Taxa"])

    assert tx.n_errors == n_errors

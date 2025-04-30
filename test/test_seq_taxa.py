"""Tests checking that the SeqTaxa specific classes work as intended."""

from logging import ERROR, INFO

import pytest
from dotmap import DotMap

from .conftest import log_check


@pytest.mark.parametrize(
    argnames="taxon_name,expected_name",
    argvalues=[
        pytest.param(
            "coli",
            "coli",
            id="simple no change ",
        ),
        pytest.param("Escherichia coli", "Escherichia coli", id="binomial no change"),
        pytest.param(
            "Escherichia coli B21", "Escherichia coli B21", id="trinomial no change"
        ),
        pytest.param(
            "Candidatus Koribacter",
            "Koribacter",
            id="candidatus simple",
        ),
        pytest.param(
            "Candidatus Koribacter versatilis",
            "Koribacter versatilis",
            id="candidatus binomial",
        ),
    ],
)
def test_remove_additional_tags(taxon_name, expected_name):
    """Test that the function to remove additional tags from taxa names works."""
    from safedata_validator.taxa import remove_additional_tags

    new_name = remove_additional_tags(taxon_name=taxon_name)

    assert new_name == expected_name


@pytest.mark.parametrize(
    "test_input,exp_name,exp_log",
    [
        (dict(name="Bacteria", rank="Kingdom"), "Bacteria", ()),
        (dict(name="k__Bacteria", rank="Kingdom"), "Bacteria", ()),
        (
            dict(name="k__Bacteria", rank="Phylum"),
            "Bacteria",
            ((ERROR, "Prefix of taxon k__Bacteria inconsistent with rank Phylum"),),
        ),
        (dict(name="p__Acidobacteria", rank="Phylum"), "Acidobacteria", ()),
        (dict(name="s__", rank="Species"), None, ()),
    ],
)
def test_taxa_strip(caplog, test_input, exp_name, exp_log):
    """Checks that the function to remove k__ type notation is functioning properly.

    This function also checks that the supplied rank matches the rank implied by the
    notation.
    """
    from safedata_validator.taxa import taxa_strip

    s_taxa = taxa_strip(**test_input)

    assert s_taxa == exp_name

    log_check(caplog, exp_log)


@pytest.mark.parametrize(
    argnames=["mock_output", "expected_log_entries"],
    argvalues=[
        pytest.param(
            DotMap({"data_columns": []}),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (ERROR, "No data or only headers in Taxa worksheet"),
            ),
            id="No data in sheet",
        ),
        pytest.param(
            DotMap({"data_columns": ["some_columns"], "bad_headers": "duplicated"}),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (ERROR, "Cannot parse taxa with duplicated headers"),
            ),
            id="Duplicated headers",
        ),
        pytest.param(
            DotMap({"data_columns": [tuple()], "headers": ["genus"]}),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (ERROR, "Sequencing taxa sheet is missing the name field"),
            ),
            id="Missing name field",
        ),
        pytest.param(
            DotMap(
                {"data_columns": [tuple(), tuple()], "headers": ["name", "comments"]}
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (ERROR, "At least one top-level taxonomic rank must be provided!"),
            ),
            id="No rank fields",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple(), tuple(), tuple(), tuple()],
                    "headers": [
                        "name",
                        "kingdom",
                        "phylum",
                        "random extra header",
                    ],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Additional fields provided:"),
                (INFO, "No taxon rows found"),
            ),
            id="Extra header",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple(), tuple()],
                    "headers": ["name", "domain", "kingdom"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "No taxon rows found"),
            ),
            id="domain + kingdom",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple(), tuple()],
                    "headers": ["name", "superkingdom", "kingdom"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "No taxon rows found"),
            ),
            id="superkingdom + kingdom",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple()],
                    "headers": ["name", "domain"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "No taxon rows found"),
            ),
            id="just domain",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple()],
                    "headers": ["name", "superkingdom"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "No taxon rows found"),
            ),
            id="just superkingdom",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple()],
                    "headers": ["name", "kingdom"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "No taxon rows found"),
            ),
            id="just kingdom",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple(), tuple()],
                    "headers": ["name", "domain", "superkingdom"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (
                    ERROR,
                    "Cannot provide both 'domain' and 'superkingdom' as taxonomic "
                    "ranks!",
                ),
            ),
            id="both domain and superkingdom",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple()],
                    "headers": ["name", "phylum"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (ERROR, "At least one top-level taxonomic rank must be provided!"),
            ),
            id="no top rank",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple(), tuple(), tuple()],
                    "headers": ["name", "superkingdom", "phylum", "order"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (
                    ERROR,
                    "Need to provide all taxonomic ranks higher than current lowest "
                    "rank (order) in SeqTaxa, missing ranks are as follows:",
                ),
            ),
            id="taxon hierarchy gap",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="good taxon",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("NA",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded phylum: Basidiomycota"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="NA leaf",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("NA",),
                        ("Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="NA mid rank",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                        ("Tremellales",),
                        ("Trimorphomycetaceae",),
                        ("NA",),
                        ("podzolica",),
                    ],
                    "headers": [
                        "name",
                        "kingdom",
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "genus",
                        "species",
                    ],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded species: <genus unknown> podzolica"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="NA genus",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                        ("Tremellales",),
                        ("Trimorphomycetaceae",),
                        (23,),
                        ("podzolica",),
                    ],
                    "headers": [
                        "name",
                        "kingdom",
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "genus",
                        "species",
                    ],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (ERROR, "Rank genus has non-string or empty string value: 23"),
                (INFO, "Loaded species: <genus unknown> podzolica"),
                (INFO, "SeqTaxa contains 1 errors"),
            ),
            id="whitespace genus",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("k__Fungi",),
                        ("p__Basidiomycota",),
                        ("c__Tremellomycetes",),
                        ("o__Tremellales",),
                        ("f__Trimorphomycetaceae",),
                        ("g__",),
                        ("s__podzolica",),
                    ],
                    "headers": [
                        "name",
                        "kingdom",
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "genus",
                        "species",
                    ],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded species: <genus unknown> podzolica"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="g__ genus",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("k__Fungi",),
                        ("p__Basidiomycota",),
                        ("c__Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="k__ notation good",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("k__Fungi",),
                        ("p__Basidiomycota",),
                        ("s__Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (
                    ERROR,
                    "Prefix of taxon s__Tremellomycetes inconsistent with rank class",
                ),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "SeqTaxa contains 1 errors"),
            ),
            id="k__ notation bad",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("k__Fungi",),
                        ("Basidiomycota",),
                        ("c__Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="k__ notation mixed",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Candidatus Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="candidatus single",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Candidatus Basidiomycota",),
                        ("Candidatus Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="candidatus multiple",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                        ("Tremellales",),
                        ("Trimorphomycetaceae",),
                        ("Saitozyma",),
                        ("Saitozyma podzolica",),
                    ],
                    "headers": [
                        "name",
                        "kingdom",
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "genus",
                        "species",
                    ],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (
                    ERROR,
                    "Provided species name appears to be a binomial (which isn't "
                    "allowed): Saitozyma podzolica",
                ),
                (INFO, "Loaded genus: Saitozyma"),
                (INFO, "SeqTaxa contains 1 errors"),
            ),
            id="species binomial",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                        ("Tremellales",),
                        ("Trimorphomycetaceae",),
                        ("NA",),
                        ("Saitozyma podzolica",),
                    ],
                    "headers": [
                        "name",
                        "kingdom",
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "genus",
                        "species",
                    ],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (
                    ERROR,
                    "Provided species name appears to be a binomial (which isn't "
                    "allowed): Saitozyma podzolica",
                ),
                (INFO, "Loaded family: Trimorphomycetaceae"),
                (INFO, "SeqTaxa contains 1 errors"),
            ),
            id="species binomial with NA",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                        ("Tremellales",),
                        ("Trimorphomycetaceae",),
                        ("Saitozyma",),
                        ("Candidatus podzolica",),
                    ],
                    "headers": [
                        "name",
                        "kingdom",
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "genus",
                        "species",
                    ],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (INFO, "Loaded species: Saitozyma podzolica"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="species with candidatus tag",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                        ("Tremellales",),
                        ("Trimorphomycetaceae",),
                        ("Saitozyma",),
                        ("Candidatus Saitozyma podzolica",),
                    ],
                    "headers": [
                        "name",
                        "kingdom",
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "genus",
                        "species",
                    ],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (
                    ERROR,
                    "Provided species name appears to be a binomial (which isn't "
                    "allowed): Saitozyma podzolica",
                ),
                (INFO, "Loaded genus: Saitozyma"),
                (INFO, "SeqTaxa contains 1 errors"),
            ),
            id="species binomial with candidatus tag",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        (562,),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: 562"),
                (ERROR, "Worksheet name is not a string: 562"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "SeqTaxa contains 1 errors"),
            ),
            id="numeric name",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100 ",),
                        ("Fungi",),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (ERROR, "Worksheet name has whitespace padding: 'ASV_100 '"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "SeqTaxa contains 1 errors"),
            ),
            id="padded name",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        (123,),
                        ("Basidiomycota",),
                        ("Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (ERROR, "Rank kingdom has non-string or empty string value: 123"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "SeqTaxa contains 1 errors"),
            ),
            id="non string taxon",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("ASV_100",),
                        ("Fungi",),
                        ("Basidiomycota ",),
                        ("Tremellomycetes",),
                    ],
                    "headers": ["name", "kingdom", "phylum", "class"],
                }
            ),
            (
                (INFO, "Loading SeqTaxa worksheet"),
                (INFO, "Reading bioinformatics taxon data"),
                (INFO, "Loading row 1: ASV_100"),
                (ERROR, "Rank phylum has whitespace padding: 'Basidiomycota '"),
                (INFO, "Loaded class: Tremellomycetes"),
                (INFO, "SeqTaxa contains 1 errors"),
            ),
            id="padded taxon",
        ),
    ],
)
def test_load_worksheet(
    caplog, mocker, fixture_resources, mock_output, expected_log_entries
):
    """Test that unexpected header names are caught by load."""
    from safedata_validator.logger import FORMATTER
    from safedata_validator.taxa import SeqTaxa

    # Create GBIFTaxa class
    tx = SeqTaxa(fixture_resources)

    # Setup mocking
    mock_get = mocker.patch("safedata_validator.taxa.GetDataFrame")
    mock_get.return_value = mock_output

    tx.load("meaningless_string")
    # This is needed to ensure that the logging depth is not altered by this test
    FORMATTER.pop()

    log_check(caplog, expected_log_entries)


@pytest.mark.parametrize(
    "example_seq_files, n_errors, n_taxa, t_taxa",
    [
        pytest.param("good", 0, 14, 49, id="good"),
        pytest.param("bad", 3, 14, 47, id="bad"),
    ],
    indirect=["example_seq_files"],  # take actual params from fixture
)
def test_taxa_load(fixture_resources, example_seq_files, n_errors, n_taxa, t_taxa):
    """This tests the ensemble loading of taxa from sequencing from a file.

    It uses indirect parametrisation to access the fixtures containing the sample excel
    files.
    """
    from safedata_validator.taxa import SeqTaxa

    tx = SeqTaxa(fixture_resources)
    tx.load(example_seq_files["SeqTaxa"])

    assert tx.n_errors == n_errors
    # Compare both named taxa and total taxa
    assert len(tx.taxon_names) == n_taxa
    assert len(tx.taxon_index) == t_taxa

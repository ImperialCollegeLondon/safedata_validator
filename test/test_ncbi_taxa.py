"""Tests checking that the NCBI specific classes work as intended."""
import copy
from contextlib import nullcontext as does_not_raise
from logging import CRITICAL, ERROR, INFO, WARNING

import pytest
from dotmap import DotMap

from safedata_validator import taxa

from .conftest import log_check

# ------------------------------------------
# Testing NCBITaxon
# ------------------------------------------


@pytest.mark.parametrize(
    "test_input,raises,message",
    [
        pytest.param(
            dict(
                name=1,
                rank="genus",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier=["genus", "Morus", 37577, 30446],
            ),
            pytest.raises(TypeError),
            "Provided taxon name not a string",
            id="non string name",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id="error",
                parent_ncbi_id=30446,
                taxa_hier=[("genus", "Morus", 37577, 30446)],
            ),
            pytest.raises(TypeError),
            "NCBI ID is not an integer",
            id="string ncbi_id",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=3.141,
                parent_ncbi_id=30446,
                taxa_hier=[("genus", "Morus", 37577, 30446)],
            ),
            pytest.raises(TypeError),
            "NCBI ID is not an integer",
            id="non-integer ncbi_id",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier="test",
            ),
            pytest.raises(TypeError),
            "Taxonomic hierarchy not a list",
            id="string instead of dictionary",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier=[],
            ),
            pytest.raises(ValueError),
            "Taxon hierarchy empty",
            id="empty dictionary",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier=[27],
            ),
            pytest.raises(ValueError),
            "Taxon hierarchy values not all tuples",
            id="non-tuple value",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier=[("genus", 37577, "Morus")],
            ),
            pytest.raises(ValueError),
            "Taxon hierarchy tuples malformed",
            id="Bad tuple length",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier=[("genus", "Morus", 37577, "extra")],
            ),
            pytest.raises(ValueError),
            "Taxon hierarchy tuples malformed",
            id="Bad tuple types",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="species",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier=[("genus", "Morus", 37577, 30446)],
            ),
            pytest.raises(ValueError),
            "Provided rank (species) not in hierarchy",
            id="rank missing",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="species",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier=[
                    ("genus", "Morus", 37577, 30446),
                    ("species", "Morus bassanus", 1, 37578),
                ],
            ),
            pytest.raises(ValueError),
            "Provided rank (species) does not match first rank in hierarchy (genus)",
            id="rank not first",
        ),
        pytest.param(
            dict(
                name="Bombus bombus",
                rank="genus",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier=[("genus", "Morus", 37577, 30446)],
            ),
            pytest.raises(ValueError),
            "Provided taxon name (Bombus bombus) does not match "
            "first name in hierarchy (Morus)",
            id="name mismatch",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=27,
                parent_ncbi_id=30446,
                taxa_hier=[("genus", "Morus", 37577, 30446)],
            ),
            pytest.raises(ValueError),
            "Provided NCBI ID (27) does not match first ID in hierarchy (37577)",
            id="taxid mismatch",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=27,
                parent_ncbi_id=666,
                taxa_hier=[("genus", "Morus", 37577, 30446)],
            ),
            pytest.raises(ValueError),
            "Provided parent NCBI ID (666) does not match "
            "first parent ID in hierarchy (30446)",
            id="parent taxid mismatch",
        ),
        pytest.param(
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=37577,
                parent_ncbi_id=30446,
                taxa_hier=[("genus", "Morus", 37577, 30446)],
            ),
            does_not_raise(),
            None,
            id="passes",
        ),
    ],
)
def test_taxon_init_errors(test_input, raises, message):
    """This test checks NCBI taxon inputs expected errors."""

    with raises as excep:
        _ = taxa.NCBITaxon(**test_input)

        if not isinstance(raises, does_not_raise):
            assert excep.value == message


# ------------------------------------------
# Testing taxa_strip
# ------------------------------------------
# Only need to test that output is sensible here
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

    s_taxa = taxa.taxa_strip(**test_input)

    assert s_taxa == exp_name

    log_check(caplog, exp_log)


# ------------------------------------------
# Testing construct_bi_or_tri
# ------------------------------------------
@pytest.mark.parametrize(
    "test_input,raises,expected,expected_log",
    [
        pytest.param(
            dict(higher_nm="Escherichia", lower_nm="coli", tri=False),
            does_not_raise(),
            "Escherichia coli",
            (),
            id="species",
        ),
        pytest.param(
            dict(higher_nm="Escherichia", lower_nm="Escherichia coli", tri=False),
            does_not_raise(),
            "Escherichia coli",
            (),
            id="species already binomial",
        ),
        pytest.param(
            dict(higher_nm="Gorilla", lower_nm="gorilla", tri=False),
            does_not_raise(),
            "Gorilla gorilla",
            (),
            id="species 2",
        ),
        pytest.param(
            dict(
                higher_nm="Candidatus Koribacter",
                lower_nm="Candidatus versatilis",
                tri=False,
            ),
            does_not_raise(),
            "Candidatus Koribacter versatilis",
            (),
            id="species both candidatus",
        ),
        pytest.param(
            dict(higher_nm="Candidatus Koribacter", lower_nm="versatilis", tri=False),
            does_not_raise(),
            "Candidatus Koribacter versatilis",
            (),
            id="species genus candidatus",
        ),
        pytest.param(
            dict(higher_nm="Over long genus name", lower_nm="vulpes", tri=False),
            pytest.raises(ValueError),
            None,
            ((ERROR, "Genus name (Over long genus name) appears to be too long"),),
            id="species bad genus",
        ),
        pytest.param(
            dict(higher_nm="Canis", lower_nm="Vulpes vulpes", tri=False),
            pytest.raises(ValueError),
            None,
            (
                (
                    ERROR,
                    "Species name (Vulpes vulpes) appears to be binomial but does not"
                    " contain genus name (Canis)",
                ),
            ),
            id="species inconsistent genus",
        ),
        pytest.param(
            dict(higher_nm="Vulpes vulpes", lower_nm="japonica", tri=True),
            does_not_raise(),
            "Vulpes vulpes japonica",
            (),
            id="subsp",
        ),
        pytest.param(
            dict(
                higher_nm="Candidatus Koribacter versatilis",
                lower_nm="Ellin345",
                tri=True,
            ),
            does_not_raise(),
            "Candidatus Koribacter versatilis Ellin345",
            (),
            id="subsp species candidatus",
        ),
        pytest.param(
            dict(
                higher_nm="Candidatus Koribacter versatilis",
                lower_nm="Candidatus Ellin345",
                tri=True,
            ),
            does_not_raise(),
            "Candidatus Koribacter versatilis Ellin345",
            (),
            id="subsp both candidatus",
        ),
        pytest.param(
            dict(
                higher_nm="Vulpes vulpes", lower_nm="Vulpes vulpes schrenckii", tri=True
            ),
            does_not_raise(),
            "Vulpes vulpes schrenckii",
            (),
            id="subsp already trinomial",
        ),
        pytest.param(
            dict(
                higher_nm="Canis vulpes", lower_nm="Vulpes vulpes schrenckii", tri=True
            ),
            pytest.raises(ValueError),
            None,
            (
                (
                    ERROR,
                    "Subspecies name (Vulpes vulpes schrenckii) appears to be trinomial"
                    " but does not contain species name (Canis vulpes)",
                ),
            ),
            id="subsp inconsistent species",
        ),
        pytest.param(
            dict(higher_nm="Over long name", lower_nm="schrenckii", tri=True),
            pytest.raises(ValueError),
            None,
            ((ERROR, "Species name (Over long name) appears to be too long"),),
            id="subsp long species",
        ),
        pytest.param(
            dict(higher_nm="Vulpes", lower_nm="Vulpes vulpes schrenckii", tri=True),
            pytest.raises(ValueError),
            None,
            ((ERROR, "Species name (Vulpes) too short"),),
            id="subsp short species",
        ),
    ],
)
def test_construct_bi_or_tri(caplog, test_input, raises, expected, expected_log):
    """Test function that constructs species binomials from species and genus names.

    We test that it can catch when the species name is already a binomial, and that it
    catches Candidatus names and handles them properly.
    """
    with raises:
        s_nm = taxa.construct_bi_or_tri(**test_input)

    if isinstance(raises, does_not_raise):
        assert s_nm == expected

    log_check(caplog, expected_log)


# ------------------------------------------
# Testing taxon validators
# ------------------------------------------
@pytest.mark.parametrize(
    "tax_id,raises, leaf,n_tax",
    [
        pytest.param(
            562,
            does_not_raise(),
            ("Escherichia coli", "species", 562, 561),
            8,
            id="Backbone leaf",
        ),
        pytest.param(
            562,
            does_not_raise(),
            ("Escherichia coli 1-110-08_S1_C1", "strain", 1444049, 562),
            9,
            id="Non backbone leaf",
        ),
        pytest.param(
            131567,
            does_not_raise(),
            ("cellular organisms", "no rank", 131567, 1),
            1,
            id="No backbone ranks",
        ),
        pytest.param(1000, pytest.raises(taxa.NCBIError), None, None, id="Merged taxa"),
        pytest.param(-1000, pytest.raises(taxa.NCBIError), None, None, id="Bad id"),
    ],
)
def test_get_canon_hierarchy(fixture_ncbi_validator, raises, tax_id, leaf, n_tax):
    """Test structure and failure modes of _get_hierarchy."""

    with raises as err:
        hier = fixture_ncbi_validator._get_canon_hierarchy(tax_id)

        if isinstance(err, does_not_raise):
            assert hier[0] == leaf
            assert len(hier) == n_tax


@pytest.mark.parametrize(
    "tax_id,expected_hier,",
    [
        pytest.param(
            562,
            [
                ("species", "Escherichia coli", 562, 561),
                ("genus", "Escherichia", 561, 543),
                ("family", "Enterobacteriaceae", 543, 91347),
                ("order", "Enterobacterales", 91347, 1236),
                ("class", "Gammaproteobacteria", 1236, 1224),
                ("phylum", "Pseudomonadota", 1224, 2),
                ("superkingdom", "Bacteria", 2, None),
            ],
            id="Backbone leaf",
        ),
        pytest.param(
            9627,
            [
                ("species", "Vulpes vulpes", 9627, 9625),
                ("genus", "Vulpes", 9625, 9608),
                ("family", "Canidae", 9608, 33554),
                ("order", "Carnivora", 33554, 40674),
                ("class", "Mammalia", 40674, 7711),
                ("phylum", "Chordata", 7711, 33208),
                ("kingdom", "Metazoa", 33208, 2759),
                ("superkingdom", "Eukaryota", 2759, None),
            ],
            id="Backbone leaf with more complexity",
        ),
        pytest.param(
            1444049,
            [
                ("strain", "Escherichia coli 1-110-08_S1_C1", 1444049, 562),
                ("species", "Escherichia coli", 562, 561),
                ("genus", "Escherichia", 561, 543),
                ("family", "Enterobacteriaceae", 543, 91347),
                ("order", "Enterobacterales", 91347, 1236),
                ("class", "Gammaproteobacteria", 1236, 1224),
                ("phylum", "Pseudomonadota", 1224, 2),
                ("superkingdom", "Bacteria", 2, None),
            ],
            id="Non backbone leaf",
        ),
        pytest.param(
            2,
            [
                ("superkingdom", "Bacteria", 2, None),
            ],
            id="Only one backbone rank",
        ),
        pytest.param(
            131567,
            [
                ("no rank", "cellular organisms", 131567, None),
            ],
            id="No backbone ranks",
        ),
    ],
)
def test_canon_to_backbone_hierarchy(fixture_ncbi_validator, tax_id, expected_hier):
    """Test structure and failure modes of _get_hierarchy."""

    canon_hier = fixture_ncbi_validator._get_canon_hierarchy(tax_id)

    bb_hier = fixture_ncbi_validator._canon_to_backbone_hierarchy(canon_hier)

    assert bb_hier == expected_hier


@pytest.mark.parametrize(
    "provided, report, exp_congruent, exp_log",
    [
        pytest.param(
            (
                ("family", "Enterobacteriaceae"),
                ("genus", "Escherichia"),
                ("species", "Escherichia coli"),
            ),
            False,
            True,
            (),
            id="good no reporting",
        ),
        pytest.param(
            (
                ("family", "Enterobacteriaceae"),
                ("genus", "Escherichia"),
                ("species", "Escherichia coli"),
            ),
            True,
            True,
            (),
            id="good with reporting",
        ),
        pytest.param(
            (),
            False,
            True,
            (),
            id="provided is empty",
        ),
        pytest.param(
            (
                ("family", "Enterobacteriaceae"),
                ("genus", "Escherichia"),
                ("species", "Escherichia coli"),
                ("strain", "Escherichia coli 1-110-08_S1_C1"),
            ),
            False,
            False,
            (),
            id="bad rank no reporting",
        ),
        pytest.param(
            (
                ("family", "Enterobacteriaceae"),
                ("genus", "Escherichia"),
                ("species", "Escherichia coli"),
                ("strain", "Escherichia coli 1-110-08_S1_C1"),
            ),
            True,
            False,
            (
                (
                    ERROR,
                    "Taxonomy mismatch for Escherichia coli 1-110-08_S1_C1 at rank "
                    "strain: rank not found in expected hierarchy",
                ),
            ),
            id="bad rank with reporting",
        ),
        pytest.param(
            (
                ("family", "Enterobacterales"),
                ("genus", "Escherichia"),
                ("species", "Escherichia coli"),
            ),
            False,
            False,
            (),
            id="bad no reporting",
        ),
        pytest.param(
            (
                ("family", "Enterobacterales"),
                ("genus", "Escherichia"),
                ("species", "Escherichia coli"),
            ),
            True,
            False,
            (
                (
                    ERROR,
                    "Taxonomy mismatch for Enterobacterales at rank family: "
                    "expecting Enterobacteriaceae",
                ),
            ),
            id="bad with reporting",
        ),
        pytest.param(
            (
                ("family", "Enterobacteraceae"),
                ("genus", "Escherichia"),
                ("species", "Enterococcus coli"),
            ),
            False,
            True,
            (),
            id="synonyms no reporting",
        ),
        pytest.param(
            (
                ("family", "Enterobacteraceae"),
                ("genus", "Escherichia"),
                ("species", "Enterococcus coli"),
            ),
            True,
            True,
            (
                (
                    WARNING,
                    "Non-canon name Enterobacteraceae at rank family: "
                    "synonym for Enterobacteriaceae",
                ),
                (
                    WARNING,
                    "Non-canon name Enterococcus coli at rank species: "
                    "synonym for Escherichia coli",
                ),
            ),
            id="synonyms with reporting",
        ),
    ],
)
def test_check_congruent_hierarchies(
    fixture_ncbi_validator, caplog, provided, report, exp_congruent, exp_log
):
    """Test  _check_congruent_hierarchies."""

    expected = fixture_ncbi_validator._get_canon_hierarchy(562)

    congruent = fixture_ncbi_validator._check_congruent_hierarchies(
        expected_hier=expected,
        provided_hier=provided,
        report=report,
    )
    assert congruent == exp_congruent

    log_check(caplog, exp_log)


# First test the search function
@pytest.mark.parametrize(
    "test_input,expected",
    [
        pytest.param(
            dict(nnme="E coli", ncbi_id=562),
            ("species", "Escherichia coli", 562, 561),
            id="species",
        ),
        pytest.param(
            dict(nnme="E coli strain", ncbi_id=1444049),
            ("strain", "Escherichia coli 1-110-08_S1_C1", 1444049, 562),
            id="strain",
        ),
        pytest.param(
            dict(nnme="Streptophytina", ncbi_id=131221),
            ("subphylum", "Streptophytina", 131221, 35493),
            id="subphylum",
        ),
        pytest.param(
            dict(nnme="Opisthokonta", ncbi_id=33154),
            ("clade", "Opisthokonta", 33154, 2759),
            id="clade",
        ),
    ],
)
def test_id_lookup(fixture_ncbi_validator, test_input, expected):
    """This test checks the results of looking up a specific NCBI taxonomy ID."""

    fnd_tx = fixture_ncbi_validator.id_lookup(**test_input)

    # Check the properties and that the leaf of the hierarchy match
    assert (fnd_tx.rank, fnd_tx.name, fnd_tx.ncbi_id) == expected[:-1]
    assert fnd_tx.taxa_hier[0] == expected


# Third function that checks that id_lookup throws the appropriate errors
@pytest.mark.parametrize(
    "test_input,expected_exception",
    [
        (None, TypeError),  # no parameters
        (dict(nnme="E coli", ncbi_id="invalid_string"), TypeError),  # a string
        (dict(nnme="E coli", ncbi_id=27.5), TypeError),  # non-integer ID
        (dict(nnme=27, ncbi_id=27), TypeError),  # named using a number
        (dict(nnme="E coli", ncbi_id=-1), ValueError),  # bad ID
        (dict(nnme="E coli", ncbi_id=100000000000000), taxa.NCBIError),  # bad ID
    ],
)
def test_id_lookup_errors(fixture_ncbi_validator, test_input, expected_exception):
    """This test checks that validator.id_lookup inputs throw errors as expected."""

    with pytest.raises(expected_exception):
        _ = fixture_ncbi_validator.id_lookup(**test_input)


# Then do the same for the taxa search function
@pytest.mark.parametrize(
    "test_input,leaf,raises,non_canon,exp_log",
    [
        pytest.param(
            dict(
                nnme="nope",
                taxon_hier=[("species", "Not a species")],
            ),
            None,
            pytest.raises(taxa.NCBIError),
            None,
            ((ERROR, "Taxa nope not found with name Not a species and rank species"),),
            id="Nothing found",
        ),
        pytest.param(
            dict(
                nnme="E coli",
                taxon_hier=[("genus", "Escherichia"), ("species", "Escherichia coli")],
            ),
            ("species", "Escherichia coli", 562, 561),
            does_not_raise(),
            (None, None),
            ((INFO, "Match found for E coli"),),
            id="Simple species",
        ),
        pytest.param(
            dict(
                nnme="Entero",
                taxon_hier=[
                    ("order", "Enterobacterales"),
                    ("family", "Enterobacteriaceae"),
                ],
            ),
            ("family", "Enterobacteriaceae", 543, 91347),
            does_not_raise(),
            (None, None),
            ((INFO, "Match found for Entero"),),
            id="Simple family",
        ),
        pytest.param(
            dict(
                nnme="E coli strain",
                taxon_hier=[
                    ("species", "Escherichia coli"),
                    ("strain", "Escherichia coli 1-110-08_S1_C1"),
                ],
            ),
            ("strain", "Escherichia coli 1-110-08_S1_C1", 1444049, 562),
            does_not_raise(),
            (None, None),
            ((INFO, "Match found for E coli strain"),),
            id="Simple strain",
        ),
        pytest.param(
            dict(
                nnme="Strepto",
                taxon_hier=[
                    ("phylum", "Streptophyta"),
                    ("subphylum", "Streptophytina"),
                ],
            ),
            ("subphylum", "Streptophytina", 131221, 35493),
            does_not_raise(),
            (None, None),
            ((INFO, "Match found for Strepto"),),
            id="Simple phylum",
        ),
        # pytest.param(
        #     dict(
        #         nnme="Opistho",
        #         taxon_hier=[("superkingdom", "Eukaryota"), ("clade", "Opisthokonta")],
        #     ),
        #     ("clade", "Opisthokonta", 33154, 2759),
        #     does_not_raise(),
        #     (None, None),
        #     ((INFO, "Match found for Opistho"),),
        #     id="Clade anchored by superkingdom",
        # ),
        pytest.param(
            dict(
                nnme="fox",
                taxon_hier=[("genus", "Vulpes"), ("species", "Canis vulpes")],
            ),
            ("species", "Vulpes vulpes", 9627, 9625),
            does_not_raise(),
            ("Canis vulpes", "synonym"),
            (
                (WARNING, "Non-canon usage: Canis vulpes is synonym for Vulpes vulpes"),
                (INFO, "Match found for fox"),
            ),
            id="Synonym with congruent taxonomy",
        ),
        pytest.param(
            dict(
                nnme="fox",
                taxon_hier=[("genus", "Canis"), ("species", "Canis vulpes")],
            ),
            None,
            pytest.raises(taxa.NCBIError),
            None,
            (
                (WARNING, "Non-canon usage: Canis vulpes is synonym for Vulpes vulpes"),
                (ERROR, "Taxonomy mismatch for Canis at rank genus: expecting Vulpes"),
                (ERROR, "Match found for fox with incongruent taxonomy"),
            ),
            id="Synonym with incongruent taxonomy",
        ),
        pytest.param(
            dict(
                nnme="T maritimum",
                taxon_hier=[
                    ("genus", "Tenacibaculum"),
                    ("species", "Tenacibaculum maritimum"),
                ],
            ),
            ("species", "Tenacibaculum maritimum", 107401, 104267),
            does_not_raise(),
            (None, None),
            ((INFO, "Match found for T maritimum"),),
            id="T maritimum simple",
        ),
        pytest.param(
            dict(
                nnme="C marina",
                taxon_hier=[
                    ("genus", "Tenacibaculum"),
                    ("species", "Cytophaga marina"),
                ],
            ),
            ("species", "Tenacibaculum maritimum", 107401, 104267),
            does_not_raise(),
            ("Cytophaga marina", "synonym"),
            (
                (
                    WARNING,
                    "Non-canon usage: Cytophaga marina is synonym for "
                    "Tenacibaculum maritimum",
                ),
                (INFO, "Match found for C marina"),
            ),
            id="T maritimum from C marina",
        ),
        pytest.param(
            dict(
                nnme="C marina",
                taxon_hier=[("genus", "Cytophaga"), ("species", "Cytophaga marina")],
            ),
            None,
            pytest.raises(taxa.NCBIError),
            None,
            (
                (
                    WARNING,
                    "Non-canon usage: Cytophaga marina is synonym for "
                    "Tenacibaculum maritimum",
                ),
                (
                    ERROR,
                    "Taxonomy mismatch for Cytophaga at rank genus: expecting "
                    "Tenacibaculum",
                ),
                (ERROR, "Match found for C marina with incongruent taxonomy"),
            ),
            id="T maritimum from C marina with incongruent taxonomy",
        ),
        pytest.param(
            dict(nnme="Morus", taxon_hier=[("family", "Moraceae"), ("genus", "Morus")]),
            ("genus", "Morus", 3497, 3487),
            does_not_raise(),
            (None, None),
            ((INFO, "Match found for Morus"),),
            id="Plant Morus resolved by family",
        ),
        pytest.param(
            dict(nnme="Morus", taxon_hier=[("family", "Sulidae"), ("genus", "Morus")]),
            ("genus", "Morus", 37577, 30446),
            does_not_raise(),
            (None, None),
            ((INFO, "Match found for Morus"),),
            id="Bird morus resolved by family",
        ),
        pytest.param(
            dict(nnme="Morus", taxon_hier=[("phylum", "Chordata"), ("genus", "Morus")]),
            ("genus", "Morus", 37577, 30446),
            does_not_raise(),
            (None, None),
            ((INFO, "Match found for Morus"),),
            id="Bird Morus resolved by phylum",
        ),
        pytest.param(
            dict(
                nnme="Morus",
                taxon_hier=[("superkingdom", "Eukaryota"), ("genus", "Morus")],
            ),
            None,
            pytest.raises(taxa.NCBIError),
            None,
            (
                (
                    ERROR,
                    "Multiple matches for taxon Morus: "
                    "provided taxonomy does not resolve candidates",
                ),
            ),
            id="Morus not resolved by taxonomy",
        ),
        pytest.param(
            dict(
                nnme="Morus",
                taxon_hier=[("family", "Canidae"), ("genus", "Morus")],
            ),
            None,
            pytest.raises(taxa.NCBIError),
            None,
            (
                (
                    ERROR,
                    "Multiple matches for taxon Morus: "
                    "provided taxonomy not congruent with candidates",
                ),
            ),
            id="Morus no candidates match taxonomy",
        ),
        pytest.param(
            dict(
                nnme="Morus",
                taxon_hier=[("family", "Moraceae"), ("genus", "mulberries")],
            ),
            ("genus", "Morus", 3497, 3487),
            does_not_raise(),
            ("mulberries", "genbank common name"),
            (
                (
                    WARNING,
                    "Non-canon usage: mulberries is genbank common name for Morus",
                ),
                (INFO, "Match found for Morus"),
            ),
            id="Plant Morus resolved by family with non-canon name",
        ),
        pytest.param(
            dict(
                nnme="Bacteria",
                taxon_hier=[
                    ("superkingdom", "Bacteria"),
                ],
            ),
            ("superkingdom", "Bacteria", 2, None),
            does_not_raise(),
            (None, None),
            ((INFO, "Match found for Bacteria"),),
            id="Bacteria superkingdom",
        ),
        pytest.param(
            dict(
                nnme="Bacteria",
                taxon_hier=[
                    ("kingdom", "Bacteria"),
                ],
            ),
            None,
            pytest.raises(taxa.NCBIError),
            None,
            (
                (
                    (
                        ERROR,
                        "Taxa Bacteria not found with name Bacteria and rank kingdom",
                    ),
                ),
            ),
            id="Bacteria kingdom",
        ),
    ],
)
def test_taxa_search(
    caplog, fixture_ncbi_validator, test_input, leaf, raises, non_canon, exp_log
):
    """This test checks the results of searching for a specific taxon."""

    with raises:
        fnd_tx = fixture_ncbi_validator.taxa_search(**test_input)

        if isinstance(raises, does_not_raise):
            # Check the properties and first taxon in hierarchy matches
            assert (fnd_tx.rank, fnd_tx.name, fnd_tx.ncbi_id) == leaf[:3]
            assert fnd_tx.taxa_hier[0] == leaf
            assert (fnd_tx.non_canon_name, fnd_tx.non_canon_name_class) == non_canon

        log_check(caplog, exp_log)


# Third function that checks that taxa_search throws the appropriate errors
@pytest.mark.parametrize(
    "test_input,expected_exception",
    [
        pytest.param(
            dict(nnme="E coli", taxon_hier=27), ValueError, id="hier not a list"
        ),
        pytest.param(
            dict(nnme=27, taxon_hier=[("species", "Escherichia coli")]),
            ValueError,
            id="nnme not a string",
        ),
        pytest.param(
            dict(nnme="E coli", taxon_hier=[]), ValueError, id="empty taxon hierarchy"
        ),
        pytest.param(
            dict(nnme="E coli", taxon_hier=[("species", 27)]),
            ValueError,
            id="non string hierarchy entry",
        ),
    ],
)
def test_taxa_search_errors(fixture_ncbi_validator, test_input, expected_exception):
    """This test checks validator.taxa_search inputs throw errors as expected."""

    with pytest.raises(expected_exception):
        _ = fixture_ncbi_validator.taxa_search(**test_input)


# ------------------------------------------
# Testing NCBITaxa
# ------------------------------------------

# Start with the validate_and_add_taxon function


# First check that expected output is recovered
@pytest.mark.parametrize(
    "test_input,counts,tx_index",
    # Basic case to begin with
    [
        pytest.param(
            dict(
                m_name="E coli",
                taxon_hier=[("genus", "Escherichia"), ("species", "Escherichia coli")],
            ),
            (dict(ntax=1, nnms=1, nhier=7)),
            [("E coli", 562, 561, "Escherichia coli", "species", "accepted")],
            id="simple species",
        ),
        pytest.param(
            dict(
                m_name="E coli",
                taxon_hier=[("genus", "Escheria"), ("species", "Escherichia coli")],
            ),
            dict(ntax=0, nnms=0, nhier=0),
            [],
            id="simple species incongruent hier",
        ),
        pytest.param(
            dict(
                m_name="C vulpes",
                taxon_hier=[("genus", "Vulpes"), ("species", "Canis vulpes")],
            ),
            dict(ntax=2, nnms=1, nhier=8),
            [
                ("C vulpes", 9627, 9625, "Canis vulpes", "species", "synonym"),
                ("C vulpes", 9627, 9625, "Vulpes vulpes", "species", "accepted"),
            ],
            id="synonym with congruent hierarchy",
        ),
        pytest.param(
            dict(
                m_name="C vulpes",
                taxon_hier=[("genus", "Cani"), ("species", "Canis vulpes")],
            ),
            dict(ntax=0, nnms=0, nhier=0),
            [],
            id="synonym with incongruent hierarchy",
        ),
        pytest.param(
            dict(
                m_name="C marina",
                taxon_hier=[("species", "Cytophaga marina")],
            ),
            dict(ntax=2, nnms=1, nhier=7),
            [
                ("C marina", 107401, 104267, "Cytophaga marina", "species", "synonym"),
                (
                    "C marina",
                    107401,
                    104267,
                    "Tenacibaculum maritimum",
                    "species",
                    "accepted",
                ),
            ],
            id="synonym with no hierarchy",
        ),
        pytest.param(
            dict(
                m_name="Bacteria",
                taxon_hier=[("kingdom", "Bacteria")],
            ),
            dict(ntax=0, nnms=0, nhier=0),
            [],
            id="Bacteria as kingdom fails",
        ),
        pytest.param(
            dict(
                m_name="Bacteria",
                taxon_hier=[("superkingdom", "Bacteria")],
            ),
            dict(ntax=1, nnms=1, nhier=1),
            [("Bacteria", 2, None, "Bacteria", "superkingdom", "accepted")],
            id="Bacteria as superkingdom",
        ),
        pytest.param(
            dict(
                m_name="Eukaryotes",
                taxon_hier=[("superkingdom", "Eukaryota")],
            ),
            dict(ntax=1, nnms=1, nhier=1),
            [("Eukaryotes", 2759, None, "Eukaryota", "superkingdom", "accepted")],
            id="Eukaryota",
        ),
        pytest.param(
            dict(
                m_name="Fungi",
                taxon_hier=[("superkingdom", "Eukaryota"), ("kingdom", "Fungi")],
            ),
            dict(ntax=1, nnms=1, nhier=2),
            [("Fungi", 4751, 2759, "Fungi", "kingdom", "accepted")],
            id="Fungi",
        ),
        pytest.param(
            dict(
                m_name="Fungi",
                taxon_hier=[("kingdom", "Fungi"), ("superkingdom", "Eukaryota")],
            ),
            dict(ntax=1, nnms=1, nhier=2),
            [("Fungi", 4751, 2759, "Fungi", "kingdom", "accepted")],
            id="Fungi reversed order",
        ),
        pytest.param(
            dict(
                m_name="New Fungi",
                taxon_hier=[
                    ("superkingdom", "Eukaryota"),
                    ("kingdom", "Fungi"),
                    ("phylum", "new fungal phylum"),
                ],
                new=True,
            ),
            dict(ntax=1, nnms=1, nhier=2),
            [("New Fungi", -1, 4751, "new fungal phylum", "phylum", "user")],
            id="Fungi new phylum",
        ),
    ],
)
def test_validate_and_add_taxon(fixture_resources, test_input, counts, tx_index):
    """Checks that NCBI taxon validation stores the expected information."""

    ncbi_instance = taxa.NCBITaxa(fixture_resources)
    ncbi_instance.validate_and_add_taxon(**test_input)

    assert len(ncbi_instance.taxon_index) == counts["ntax"]
    assert len(ncbi_instance.taxon_names) == counts["nnms"]
    assert len(ncbi_instance.hierarchy) == counts["nhier"]

    # Check that taxon info recorded is as expected
    for idx_entry in tx_index:
        # Check that provided taxon name is used
        assert idx_entry[0] in ncbi_instance.taxon_names
        assert idx_entry in ncbi_instance.taxon_index


@pytest.mark.parametrize(
    "test_input,raises,expected_log_entries",
    [
        pytest.param(
            ["E coli", [("species", "Escherichia coli")]],
            does_not_raise(),
            ((INFO, "Match found for E coli"),),
            id="good species",
        ),
        pytest.param(
            ["Eukaryota", [("superkingdom", "Eukaryota")]],
            does_not_raise(),
            ((INFO, "Match found for Eukaryota"),),
            id="good superkingdom",
        ),
        pytest.param(
            ["Bacteria", [("superkingdom", "Bacteria")]],
            does_not_raise(),
            ((INFO, "Match found for Bacteria"),),
            id="bacteria as superkingdom",
        ),
        pytest.param(
            ["Bacteria", [("kingdom", "Bacteria")]],
            does_not_raise(),
            (
                (ERROR, "Taxa Bacteria not found with name Bacteria and rank kingdom"),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
            id="bacteria as kingdom fails",
        ),
        pytest.param(
            [" E coli", [("species", "Escherichia coli")]],
            does_not_raise(),
            (
                (ERROR, "Worksheet name has whitespace padding: ' E coli'"),
                (INFO, "Match found for E coli"),
            ),
            id="whitespace padding error",
        ),
        pytest.param(
            ["    ", [("species", "Escherichia coli")]],
            does_not_raise(),
            ((ERROR, "Worksheet name missing, whitespace only or not text"),),
            id="whitespace only name",
        ),
        pytest.param(
            [None, [("species", "Escherichia coli")]],
            does_not_raise(),
            ((ERROR, "Worksheet name missing, whitespace only or not text"),),
            id="no name empty cell",
        ),
        pytest.param(
            ["E coli", "Escherichia coli"],
            pytest.raises(ValueError),
            ((CRITICAL, "Taxon hierarchy should be a list"),),
            id="non list hierarchy",
        ),
        pytest.param(
            ["E coli", []],
            does_not_raise(),
            ((ERROR, "No taxonomy provided"),),
            id="empty list hierarchy",
        ),
        pytest.param(
            [" E coli", []],
            does_not_raise(),
            (
                (ERROR, "Worksheet name has whitespace padding: ' E coli'"),
                (ERROR, "No taxonomy provided"),
            ),
            id="multiple errors",
        ),
        pytest.param(
            ["E coli", ["notatuple"], None],
            pytest.raises(ValueError),
            (
                (
                    CRITICAL,
                    "Not all taxa hierachy entries are 2 tuples of non-empty strings",
                ),
            ),
            id="non tuple hierarchy element",
        ),
        pytest.param(
            ["E coli", [("  ", "Escherichia coli")]],
            pytest.raises(ValueError),
            (
                (
                    CRITICAL,
                    "Not all taxa hierachy entries are 2 tuples of non-empty strings",
                ),
            ),
            id="only whitespace in hierarchy tuple element",
        ),
        pytest.param(
            ["E coli", [(26, "Escherichia coli")]],
            pytest.raises(ValueError),
            (
                (
                    CRITICAL,
                    "Not all taxa hierachy entries are 2 tuples of non-empty strings",
                ),
            ),
            id="non string hierarchy tuple element",
        ),
        pytest.param(
            ["E coli", [(" species", "Escherichia coli")]],
            does_not_raise(),
            (
                (
                    ERROR,
                    "Hierarchy contains whitespace: "
                    "rank ' species', name 'Escherichia coli'",
                ),
                (INFO, "Match found for E coli"),
            ),
            id="whitespace padding in hierarchy",
        ),
        # # Taxon hierarchy in wrong order
        # (
        #     ["E coli", {"species": "Escherichia coli", "genus": "Escherichia"}, None],
        #     does_not_raise(),
        #     ((ERROR, "Taxon hierarchy not in correct order"),),
        # ),
        pytest.param(
            ["Morus", [("genus", "Morus")]],
            does_not_raise(),
            (
                (
                    ERROR,
                    "Multiple matches for taxon Morus: "
                    "no additional taxonomy provided.",
                ),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
            id="ambiguous taxon",
        ),
        pytest.param(
            ["Morus", [("family", "Moraceae"), ("genus", "Morus")]],
            does_not_raise(),
            ((INFO, "Match found for Morus"),),
            id="ambiguous taxon resolved",
        ),
        pytest.param(
            [
                "E coli strain",
                [
                    ("species", "Escherichia coli"),
                    ("strain", "Escherichia coli 1-110-08_S1_C1"),
                ],
            ],
            does_not_raise(),
            ((INFO, "Match found for E coli strain"),),
            id="non backbone",
        ),
        pytest.param(
            ["C marina", [("species", "Cytophaga marina")]],
            does_not_raise(),
            (
                (
                    WARNING,
                    "Non-canon usage: Cytophaga marina is synonym for "
                    "Tenacibaculum maritimum",
                ),
                (INFO, "Match found for C marina"),
            ),
            id="synonym",
        ),
        pytest.param(
            [
                "E coli",
                [
                    ("family", "Escherichia coli"),
                ],
            ],
            does_not_raise(),
            (
                (
                    ERROR,
                    "Taxa E coli not found with name Escherichia coli and rank family",
                ),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
            id="incorrect rank assignment",
        ),
        pytest.param(
            [
                "Utter nonsense",
                [("species", "Nonsense species")],
            ],
            does_not_raise(),
            (
                (
                    ERROR,
                    "Taxa Utter nonsense not found with name "
                    "Nonsense species and rank species",
                ),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
            id="invalid taxon",
        ),
        pytest.param(
            ["E coli", [("genus", "Escheria"), ("species", "Escherichia coli")]],
            does_not_raise(),
            (
                (
                    ERROR,
                    "Taxonomy mismatch for Escheria at rank genus: expecting "
                    "Escherichia",
                ),
                (ERROR, "Match found for E coli with incongruent taxonomy"),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
            id="invalid higher taxa",
        ),
    ],
)
def test_validate_and_add_taxon_logging(
    caplog, fixture_resources, test_input, raises, expected_log_entries
):
    """Checks the logging of NCBI taxon validation."""

    ncbi_instance = taxa.NCBITaxa(fixture_resources)

    with raises:
        ncbi_instance.validate_and_add_taxon(*test_input)

    log_check(caplog, expected_log_entries)


# First test whether sensible output is produced
@pytest.mark.parametrize(
    "test_input,n_ind,exp_index,exp_log",
    [
        pytest.param(
            ["E coli", [("genus", "Escherichia"), ("species", "Escherichia coli")]],
            7,
            [
                ("E coli", 562, 561, "Escherichia coli", "species", "accepted"),
                (None, 2, None, "Bacteria", "superkingdom", "accepted"),
                (None, 1224, 2, "Pseudomonadota", "phylum", "accepted"),
                (None, 1236, 1224, "Gammaproteobacteria", "class", "accepted"),
                (None, 91347, 1236, "Enterobacterales", "order", "accepted"),
                (None, 543, 91347, "Enterobacteriaceae", "family", "accepted"),
                (None, 561, 543, "Escherichia", "genus", "accepted"),
            ],
            (
                (INFO, "Match found for E coli"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Pseudomonadota"),
                (INFO, "Added class Gammaproteobacteria"),
                (INFO, "Added order Enterobacterales"),
                (INFO, "Added family Enterobacteriaceae"),
                (INFO, "Added genus Escherichia"),
            ),
            id="index one accepted species",
        ),
        pytest.param(
            ["C marina", [("species", "Cytophaga marina")]],
            8,
            [
                ("C marina", 107401, 104267, "Cytophaga marina", "species", "synonym"),
                (
                    "C marina",
                    107401,
                    104267,
                    "Tenacibaculum maritimum",
                    "species",
                    "accepted",
                ),
                (None, 2, None, "Bacteria", "superkingdom", "accepted"),
                (None, 976, 2, "Bacteroidota", "phylum", "accepted"),
                (None, 117743, 976, "Flavobacteriia", "class", "accepted"),
                (None, 200644, 117743, "Flavobacteriales", "order", "accepted"),
                (None, 49546, 200644, "Flavobacteriaceae", "family", "accepted"),
                (None, 104267, 49546, "Tenacibaculum", "genus", "accepted"),
            ],
            (
                (
                    WARNING,
                    "Non-canon usage: Cytophaga marina is synonym for "
                    "Tenacibaculum maritimum",
                ),
                (INFO, "Match found for C marina"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Bacteroidota"),
                (INFO, "Added class Flavobacteriia"),
                (INFO, "Added order Flavobacteriales"),
                (INFO, "Added family Flavobacteriaceae"),
                (INFO, "Added genus Tenacibaculum"),
            ),
            id="index one synonym",
        ),
        pytest.param(
            ["Bacteria", [("superkingdom", "Bacteria")]],
            1,
            [("Bacteria", 2, None, "Bacteria", "superkingdom", "accepted")],
            (
                (INFO, "Match found for Bacteria"),
                (INFO, "Indexing taxonomic hierarchy"),
            ),
            id="bacteria as superkingdom",
        ),
        pytest.param(
            ["Fungi", [("superkingdom", "Eukaryota"), ("kingdom", "Fungi")]],
            2,
            [
                ("Fungi", 4751, 2759, "Fungi", "kingdom", "accepted"),
                (None, 2759, None, "Eukaryota", "superkingdom", "accepted"),
            ],
            (
                (INFO, "Match found for Fungi"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Eukaryota"),
            ),
            id="fungi",
        ),
    ],
)
def test_index_higher_taxa(
    caplog, fixture_resources, test_input, n_ind, exp_index, exp_log
):
    """Checks that higher taxonomic ranks for a taxon are stored correctly."""

    ncbi_instance = taxa.NCBITaxa(fixture_resources)
    ncbi_instance.validate_and_add_taxon(*test_input)

    # Then index higher taxa
    ncbi_instance.index_higher_taxa()

    assert len(ncbi_instance.taxon_index) == n_ind

    for item in exp_index:
        assert item in ncbi_instance.taxon_index

    log_check(caplog, exp_log)


@pytest.mark.parametrize(
    argnames=["mock_output", "expected_log_entries"],
    argvalues=[
        pytest.param(
            DotMap({"data_columns": []}),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (ERROR, "No data or only headers in Taxa worksheet"),
            ),
            id="No data in sheet",
        ),
        pytest.param(
            DotMap({"data_columns": ["some_columns"], "bad_headers": "duplicated"}),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (ERROR, "Cannot parse taxa with duplicated headers"),
            ),
            id="Duplicated headers",
        ),
        pytest.param(
            DotMap({"data_columns": [tuple()], "headers": ["genus"]}),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (ERROR, "NCBI taxa sheet is missing the name field"),
            ),
            id="Missing name field",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple(), tuple(), tuple(), tuple()],
                    "headers": [
                        "name",
                        "genus",
                        "species",
                        "random extra header",
                    ],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "Additional fields provided:"),
                (INFO, "2 NCBI rank fields found:"),
                (INFO, "No taxon rows found"),
            ),
            id="Extra header",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple(), tuple()],
                    "headers": ["name", "family", "species"],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "2 NCBI rank fields found:"),
                (ERROR, "A genus field is required with species"),
            ),
            id="species without genus",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [tuple(), tuple(), tuple(), tuple()],
                    "headers": ["name", "family", "genus", "subspecies"],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "3 NCBI rank fields found:"),
                (ERROR, "A species field is required with subspecies"),
            ),
            id="subspecies without species",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("E coli",),
                        ("Enterobacteriaceae",),
                        ("Escherichia",),
                        ("coli",),
                    ],
                    "headers": ["name", "family", "genus", "species"],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "3 NCBI rank fields found:"),
                (INFO, "Validating row 1: E coli"),
                (INFO, "Match found for E coli"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Pseudomonadota"),
                (INFO, "Added class Gammaproteobacteria"),
                (INFO, "Added order Enterobacterales"),
                (INFO, "Added family Enterobacteriaceae"),
                (INFO, "Added genus Escherichia"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="good taxon",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("E coli",),
                        ("f__Enterobacteriaceae",),
                        ("g__Escherichia",),
                        ("s__coli",),
                    ],
                    "headers": ["name", "family", "genus", "species"],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "3 NCBI rank fields found:"),
                (INFO, "Validating row 1: E coli"),
                (INFO, "Match found for E coli"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Pseudomonadota"),
                (INFO, "Added class Gammaproteobacteria"),
                (INFO, "Added order Enterobacterales"),
                (INFO, "Added family Enterobacteriaceae"),
                (INFO, "Added genus Escherichia"),
                (INFO, "1 taxa loaded correctly"),
            ),
            id="good taxon with NCBI prefixes",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("E coli",),
                        ("f__Enterobacteriaceae",),
                        ("k__Escherichia",),
                        ("s__coli",),
                    ],
                    "headers": ["name", "family", "genus", "species"],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "3 NCBI rank fields found:"),
                (INFO, "Validating row 1: E coli"),
                (ERROR, "Prefix of taxon k__Escherichia inconsistent with rank genus"),
                (INFO, "Match found for E coli"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Pseudomonadota"),
                (INFO, "Added class Gammaproteobacteria"),
                (INFO, "Added order Enterobacterales"),
                (INFO, "Added family Enterobacteriaceae"),
                (INFO, "Added genus Escherichia"),
                (INFO, "NCBITaxa contains 1 errors"),
            ),
            id="taxon with bad NCBI prefixes",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        (562,),
                        ("Enterobacteriaceae",),
                        ("Escherichia",),
                        ("coli",),
                    ],
                    "headers": ["name", "family", "genus", "species"],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "3 NCBI rank fields found:"),
                (INFO, "Validating row 1: 562"),
                (ERROR, "Worksheet name is not a string: 562"),
                (INFO, "Match found for 562"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Pseudomonadota"),
                (INFO, "Added class Gammaproteobacteria"),
                (INFO, "Added order Enterobacterales"),
                (INFO, "Added family Enterobacteriaceae"),
                (INFO, "Added genus Escherichia"),
                (INFO, "NCBITaxa contains 1 errors"),
            ),
            id="numeric name",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("E coli ",),
                        ("Enterobacteriaceae",),
                        ("Escherichia",),
                        ("coli",),
                    ],
                    "headers": ["name", "family", "genus", "species"],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "3 NCBI rank fields found:"),
                (INFO, "Validating row 1:"),
                (ERROR, "Worksheet name has whitespace padding: 'E coli '"),
                (INFO, "Match found for E coli"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Pseudomonadota"),
                (INFO, "Added class Gammaproteobacteria"),
                (INFO, "Added order Enterobacterales"),
                (INFO, "Added family Enterobacteriaceae"),
                (INFO, "Added genus Escherichia"),
                (INFO, "NCBITaxa contains 1 errors"),
            ),
            id="padded name",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("E coli",),
                        (123,),
                        ("Escherichia",),
                        ("coli",),
                    ],
                    "headers": ["name", "family", "genus", "species"],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "3 NCBI rank fields found:"),
                (INFO, "Validating row 1:"),
                (ERROR, "Rank family has non-string or empty string value: 123"),
                (INFO, "NCBITaxa contains 1 errors"),
            ),
            id="non string taxon",
        ),
        pytest.param(
            DotMap(
                {
                    "data_columns": [
                        ("E coli",),
                        ("Enterobacteriaceae",),
                        ("Escherichia ",),
                        ("coli",),
                    ],
                    "headers": ["name", "family", "genus", "species"],
                }
            ),
            (
                (INFO, "Loading NCBITaxa worksheet"),
                (INFO, "Reading NCBI taxa data"),
                (INFO, "3 NCBI rank fields found:"),
                (INFO, "Validating row 1: E coli"),
                (ERROR, "Rank genus has whitespace padding: 'Escherichia '"),
                (INFO, "Match found for E coli"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Pseudomonadota"),
                (INFO, "Added class Gammaproteobacteria"),
                (INFO, "Added order Enterobacterales"),
                (INFO, "Added family Enterobacteriaceae"),
                (INFO, "Added genus Escherichia"),
                (INFO, "NCBITaxa contains 1 errors"),
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

    # Create GBIFTaxa class
    tx = taxa.NCBITaxa(fixture_resources)

    # Setup mocking
    mock_get = mocker.patch("safedata_validator.taxa.GetDataFrame")
    mock_get.return_value = mock_output

    tx.load("meaningless_string")
    # This is needed to ensure that the logging depth is not altered by this test
    FORMATTER.pop()

    log_check(caplog, expected_log_entries)


# Finally check load function (starting by doing an error overview)
@pytest.mark.parametrize(
    "example_ncbi_files, n_errors, n_taxa, t_taxa",
    [
        pytest.param("good", 0, 11, 30, id="good"),
        pytest.param("weird", 0, 5, 19, id="weird"),
        pytest.param("bad", 12, 0, 0, id="bad"),
    ],
    indirect=["example_ncbi_files"],  # take actual params from fixture
)
def test_taxa_load(fixture_resources, example_ncbi_files, n_errors, n_taxa, t_taxa):
    """This tests the ensemble loading of (ncbi) taxa from a file.

    It uses indirect parametrisation to access the fixtures containing the sample excel
    files.
    """

    tx = taxa.NCBITaxa(fixture_resources)
    tx.load(example_ncbi_files["NCBITaxa"])

    assert tx.n_errors == n_errors
    # Compare both named taxa and total taxa
    assert len(tx.taxon_names) == n_taxa
    assert len(tx.taxon_index) == t_taxa

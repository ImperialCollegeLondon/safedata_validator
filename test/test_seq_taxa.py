"""Tests checking that the SeqTaxa specific classes work as intended."""

import pytest


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

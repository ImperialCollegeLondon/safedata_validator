import copy
import os
from logging import ERROR, INFO, WARNING

import pytest

from safedata_validator import taxa
from safedata_validator.resources import Resources

from .conftest import log_check


# Check that SDV_NO_REMOTE has not been set to an invalid parameter
def test_sdv_no_remote_set_correctly():
    assert os.getenv("SDV_NO_REMOTE") in [
        None,
        "0",
        "1",
        "false",
        "true",
        "False",
        "True",
    ]


# ------------------------------------------
# Testing NCBITaxon
# ------------------------------------------


@pytest.mark.parametrize(
    "test_input,expected_exception",
    [
        (dict(), TypeError),  # no parameters
        # non string name
        (
            dict(
                name=1,
                rank="genus",
                ncbi_id=37577,
                taxa_hier={"genus": ("Morus", 37577)},
            ),
            TypeError,
        ),
        # non-numeric ncbi_id
        (
            dict(
                name="Morus",
                rank="genus",
                ncbi_id="error",
                taxa_hier={"genus": ("Morus", 37577)},
            ),
            TypeError,
        ),
        # non-integer ncbi_id
        (
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=3.141,
                taxa_hier={"genus": ("Morus", 37577)},
            ),
            TypeError,
        ),
        # string instead of dictionary
        (dict(name="Morus", rank="genus", ncbi_id=37577, taxa_hier="test"), TypeError),
        # empty dictionary
        (dict(name="Morus", rank="genus", ncbi_id=37577, taxa_hier={}), ValueError),
        # non-string key
        (
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=37577,
                taxa_hier={1: ("Morus", 37577)},
            ),
            ValueError,
        ),
        (
            dict(
                name="Morus", rank="genus", ncbi_id=37577, taxa_hier={"genus": 27}
            ),  # non-tuple value
            ValueError,
        ),
        # non-tuple value
        (
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=37577,
                taxa_hier={"genus": (37577, "Morus")},
            ),
            ValueError,
        ),
        # 3 elements in tuple
        (
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=37577,
                taxa_hier={"genus": ("Morus", 37577, "extra")},
            ),
            ValueError,
        ),
        # supplied rank doesn't match
        (
            dict(
                name="Morus",
                rank="species",
                ncbi_id=37577,
                taxa_hier={"genus": ("Morus", 37577)},
            ),
            ValueError,
        ),
        # supplied name doesn't match
        (
            dict(
                name="Bombus bombus",
                rank="genus",
                ncbi_id=37577,
                taxa_hier={"genus": ("Morus", 37577)},
            ),
            ValueError,
        ),
        # supplied NCBI ID doesn't match
        (
            dict(
                name="Morus",
                rank="genus",
                ncbi_id=27,
                taxa_hier={"genus": ("Morus", 37577)},
            ),
            ValueError,
        ),
    ],
)
def test_taxon_init_errors(test_input, expected_exception):
    """This test checks NCBI taxon inputs expected errors"""
    with pytest.raises(expected_exception):
        _ = taxa.NCBITaxon(**test_input)


# ------------------------------------------
# Testing taxa_strip
# ------------------------------------------
# Only need to test that output is sensible here
@pytest.mark.parametrize(
    "test_input,expected",
    [
        (dict(name="Bacteria", rank="Kingdom"), ("Bacteria", True)),
        (dict(name="k__Bacteria", rank="Kingdom"), ("Bacteria", True)),
        (dict(name="k__Bacteria", rank="Phylum"), ("Bacteria", False)),
        (dict(name="p__Acidobacteria", rank="Phylum"), ("Acidobacteria", True)),
        (dict(name="s__", rank="Species"), ("", True)),
        (dict(name="s__", rank="Species"), ("", True)),
        (dict(name=None, rank="Order"), (None, True)),
    ],
)
def test_taxa_strip(test_input, expected):
    """This test checks the function that strips taxa strings down to remove k__
    type notation is functioning properly. This function also checks that the
    supplied rank matches the rank implied by the notation.
    """

    s_taxa, match = taxa.taxa_strip(**test_input)

    assert s_taxa == expected[0]
    assert match == expected[1]


# ------------------------------------------
# Testing species_binomial
# ------------------------------------------
@pytest.mark.parametrize(
    "test_input,expected",
    [
        (dict(genus="Escherichia", species="coli"), "Escherichia coli"),
        (dict(genus="Escherichia", species="Escherichia coli"), "Escherichia coli"),
        (dict(genus="Gorilla", species="gorilla"), "Gorilla gorilla"),
        (
            dict(genus="Candidatus Koribacter", species="Candidatus versatilis"),
            "Candidatus Koribacter versatilis",
        ),
        (
            dict(genus="Candidatus Koribacter", species="versatilis"),
            "Candidatus Koribacter versatilis",
        ),
        (dict(genus="Over long genus name", species="vulpes"), None),
        (dict(genus="Canis", species="Vulpes vulpes"), None),
    ],
)
def test_species_binomial(test_input, expected):
    """This test checks the function that constructs species binomials from species
    and genus names. We test that it can catch when the species name is already
    a binomial, and that it catches Candidatus names and handles them properly.
    """

    s_bi = taxa.species_binomial(**test_input)

    assert s_bi == expected


# Now test that the function logs errors correctly
@pytest.mark.parametrize(
    "test_input,expected_log_entries",
    [
        (dict(genus="Escherichia", species="coli"), ()),  # Fine so empty
        (
            dict(genus="Over long genus name", species="vulpes"),  # Over long name
            ((ERROR, "Genus name (Over long genus name) appears to be too long"),),
        ),
        (
            dict(genus="Canis", species="Vulpes vulpes"),  # Genus name not in binomial
            (
                (
                    ERROR,
                    "Species name (Vulpes vulpes) appears to be binomial but does not"
                    " contain genus name (Canis)",
                ),
            ),
        ),
    ],
)
def test_validate_species_binomial(caplog, test_input, expected_log_entries):
    """This test checks that the function to construct species binomials logs the
    correct errors and warnings.
    """

    taxa.species_binomial(**test_input)

    log_check(caplog, expected_log_entries)


# ------------------------------------------
# Testing subspecies_trinomial
# ------------------------------------------
@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            dict(species="Vulpes vulpes", subspecies="japonica"),
            "Vulpes vulpes japonica",
        ),
        (
            dict(species="Candidatus Koribacter versatilis", subspecies="Ellin345"),
            "Candidatus Koribacter versatilis Ellin345",
        ),
        (
            dict(
                species="Candidatus Koribacter versatilis",
                subspecies="Candidatus Ellin345",
            ),
            "Candidatus Koribacter versatilis Ellin345",
        ),
        (
            dict(species="Vulpes vulpes", subspecies="Vulpes vulpes schrenckii"),
            "Vulpes vulpes schrenckii",
        ),
        (dict(species="Canis vulpes", subspecies="Vulpes vulpes schrenckii"), None),
        (dict(species="Over long name", subspecies="schrenckii"), None),
        (dict(species="Vulpes", subspecies="Vulpes vulpes schrenckii"), None),
    ],
)
def test_subspecies_trinomial(test_input, expected):
    """This test checks the function that constructs subspecies trinomials from
    subspecies and species names. We test that it can catch when the subspecies
    name is already a trinomial, and that it catches Candidatus names and handles
    them properly.
    """

    s_bi = taxa.subspecies_trinomial(**test_input)

    assert s_bi == expected


# Now test that the function logs errors correctly
@pytest.mark.parametrize(
    "test_input,expected_log_entries",
    [
        (
            dict(
                species="Vulpes vulpes", subspecies="Vulpes vulpes japonica"
            ),  # Fine so empty
            (),
        ),
        # Species name and genus name don't match
        (
            dict(species="Canis vulpes", subspecies="Vulpes vulpes schrenckii"),
            (
                (
                    ERROR,
                    "Subspecies name (Vulpes vulpes schrenckii) appears to be trinomial"
                    "but does not contain species name (Canis vulpes)",
                ),
            ),
        ),
        (
            dict(species="Over long name", subspecies="schrenckii"),  # Over long name
            ((ERROR, "Species name (Over long name) too long"),),
        ),
        (
            dict(species="Vulpes", subspecies="Vulpes vulpes schrenckii"),  # Too short
            ((ERROR, "Species name (Vulpes) too short"),),
        ),
    ],
)
def test_validate_subspecies_trinomial(caplog, test_input, expected_log_entries):
    """This test checks that the function to construct subspecies trinomials logs
    the correct errors and warnings.
    """

    taxa.subspecies_trinomial(**test_input)

    log_check(caplog, expected_log_entries)


# ------------------------------------------
# Testing taxon validators
# ------------------------------------------

# First test the search function
@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            dict(nnme="E coli", ncbi_id=562),
            ("Escherichia coli", "species", 562, False, None, 561),
        ),
        (
            dict(nnme="E coli strain", ncbi_id=1444049),
            ("Escherichia coli 1-110-08_S1_C1", "strain", 1444049, False, None, 562),
        ),
        (
            dict(nnme="Streptophytina", ncbi_id=131221),
            ("Streptophytina", "subphylum", 131221, False, None, 35493),
        ),
        (
            dict(nnme="Opisthokonta", ncbi_id=33154),
            ("Opisthokonta", "clade", 33154, False, None, 2759),
        ),
        (
            dict(nnme="Cytophaga marina", ncbi_id=1000),
            ("Tenacibaculum maritimum", "species", 107401, True, None, 104267),
        ),
    ],
)
def test_id_lookup(fixture_ncbi_validators, test_input, expected):
    """This test checks the expected results of looking up a specific NCBI
    taxonomy ID against the actual returned taxonomy details.
    """

    fnd_tx = fixture_ncbi_validators.id_lookup(**test_input)

    assert fnd_tx.name == expected[0]
    assert fnd_tx.rank == expected[1]
    assert fnd_tx.ncbi_id == expected[2]
    assert fnd_tx.superseed == expected[3]

    # Find last dictionary key
    f_key = list(fnd_tx.taxa_hier.keys())[-1]

    assert f_key == expected[1]

    assert fnd_tx.taxa_hier[f_key] == (expected[0], expected[2], expected[5])
    assert fnd_tx.orig == expected[4]


# Now test that the search function logs errors correctly
@pytest.mark.parametrize(
    "test_input,expected_log_entries",
    [
        (dict(nnme="E coli", ncbi_id=562), ()),  # Fine so empty
        (
            dict(nnme="Streptophytina", ncbi_id=131221),  # Non-backbone rank
            ((WARNING, "Streptophytina of non-backbone rank: subphylum"),),
        ),
        (
            dict(nnme="C marina", ncbi_id=1000),  # Non-backbone rank
            ((WARNING, "NCBI ID 1000 has been superseded by ID 107401"),),
        ),
        (
            dict(nnme="Cells", ncbi_id=131567),  # No backbone ranks whatsoever
            ((ERROR, "Taxon hierarchy for Cells contains no backbone ranks"),),
        ),
    ],
)
def test_validate_id_lookup(
    caplog, test_input, expected_log_entries, fixture_ncbi_validators
):
    """This test checks that the function to search for a taxon by NCBI ID logs
    the correct errors and warnings.
    """

    fixture_ncbi_validators.id_lookup(**test_input)

    log_check(caplog, expected_log_entries)


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
def test_id_lookup_errors(fixture_ncbi_validators, test_input, expected_exception):
    """This test checks that validator.id_lookup inputs throw errors as expected"""

    with pytest.raises(expected_exception):
        _ = fixture_ncbi_validators.id_lookup(**test_input)


# Then do the same for the taxa search function
@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            dict(
                nnme="E coli",
                taxah={"genus": "Escherichia", "species": "Escherichia coli"},
            ),
            ("Escherichia coli", "species", 562, False, None, 561),
        ),
        (
            dict(
                nnme="Entero",
                taxah={"order": "Enterobacterales", "family": "Enterobacteriaceae"},
            ),
            ("Enterobacteriaceae", "family", 543, False, None, 91347),
        ),
        (
            dict(
                nnme="E coli strain",
                taxah={
                    "species": "Escherichia coli",
                    "strain": "Escherichia coli 1-110-08_S1_C1",
                },
            ),
            ("Escherichia coli 1-110-08_S1_C1", "strain", 1444049, False, None, 562),
        ),
        (
            dict(
                nnme="Strepto",
                taxah={"phylum": "Streptophyta", "subphylum": "Streptophytina"},
            ),
            ("Streptophytina", "subphylum", 131221, False, None, 35493),
        ),
        (
            dict(
                nnme="Opistho",
                taxah={"superkingdom": "Eukaryota", "clade": "Opisthokonta"},
            ),
            ("Opisthokonta", "clade", 33154, False, None, 2759),
        ),
        (
            dict(
                nnme="Vulpes vulpes",
                taxah={"genus": "Vulpes", "species": "Vulpes vulpes"},
            ),
            ("Vulpes vulpes", "species", 9627, False, None, 9625),
        ),
        (
            dict(nnme="M morus", taxah={"family": "Moraceae", "genus": "Morus"}),
            ("Morus", "genus", 3497, False, None, 3487),
        ),
        (
            dict(nnme="S morus", taxah={"family": "Sulidae", "genus": "Morus"}),
            ("Morus", "genus", 37577, False, None, 30446),
        ),
        (
            dict(nnme="C morus", taxah={"phylum": "Chordata", "genus": "Morus"}),
            ("Morus", "genus", 37577, False, None, 30446),
        ),
        (
            dict(
                nnme="T maritimum",
                taxah={"genus": "Tenacibaculum", "species": "Tenacibaculum maritimum"},
            ),
            ("Tenacibaculum maritimum", "species", 107401, False, None, 104267),
        ),
        (
            dict(
                nnme="C marina",
                taxah={"genus": "Cytophaga", "species": "Cytophaga marina"},
            ),
            ("Tenacibaculum maritimum", "species", 107401, True, None, 104267),
        ),
        (
            dict(nnme="Bacteria", taxah={"superkingdom": "Bacteria"}),
            ("Bacteria", "superkingdom", 2, False, None, None),
        ),
        (
            dict(
                nnme="Unknown strain",
                taxah={"species": "Escherichia coli", "strain": "Nonsense strain"},
            ),
            ("Escherichia coli", "species", 562, False, "strain", 561),
        ),
    ],
)
def test_taxa_search(fixture_ncbi_validators, test_input, expected):
    """This test checks the results of searching for a taxa in the NCBI taxonomy
    database against what was expected to be found.
    """
    fnd_tx = fixture_ncbi_validators.taxa_search(**test_input)

    assert fnd_tx.name == expected[0]
    assert fnd_tx.rank == expected[1]
    assert fnd_tx.ncbi_id == expected[2]
    assert fnd_tx.superseed == expected[3]

    # Find last dictionary key
    f_key = list(fnd_tx.taxa_hier.keys())[-1]

    assert f_key == expected[1]

    assert fnd_tx.taxa_hier[f_key] == (expected[0], expected[2], expected[5])
    assert fnd_tx.orig == expected[4]


# Now test that the search function logs errors correctly
@pytest.mark.parametrize(
    "test_input,expected_log_entries",
    # Fine so empty
    [
        (
            dict(
                nnme="E coli",
                taxah={"genus": "Escherichia", "species": "Escherichia coli"},
            ),
            (),
        ),
        (
            dict(nnme="Nonsense", taxah={"species": "Nonsense garbage"}),  # Nonsense
            ((ERROR, "Taxa Nonsense cannot be found"),),
        ),
        (
            dict(nnme="Morus", taxah={"genus": "Morus"}),  # ambiguous Morus
            (
                (
                    ERROR,
                    "Taxa Morus cannot be found using only one taxonomic level,"
                    " more should be provided",
                ),
            ),
        ),
        # Nonsense Morus
        (
            dict(nnme="N Morus", taxah={"family": "Nonsense", "genus": "Morus"}),
            ((ERROR, "Provided parent taxa for N Morus not found"),),
        ),
        # Aegle Morus
        (
            dict(nnme="A Morus", taxah={"family": "Aegle", "genus": "Morus"}),
            ((ERROR, "More than one possible parent taxa for A Morus found"),),
        ),
        # Eukaryote Morus
        (
            dict(nnme="E Morus", taxah={"superkingdom": "Eukaryota", "genus": "Morus"}),
            (
                (
                    ERROR,
                    "Parent taxa for E Morus refers to multiple possible child taxa",
                ),
            ),
        ),
        # Carnivora Morus
        (
            dict(nnme="C Morus", taxah={"superkingdom": "Carnivora", "genus": "Morus"}),
            ((ERROR, "Parent taxa not actually a valid parent of C Morus"),),
        ),
        # Cytophaga marina
        (
            dict(
                nnme="C marina",
                taxah={"genus": "Cytophaga", "species": "Cytophaga marina"},
            ),
            (
                (
                    WARNING,
                    "Cytophaga marina not accepted usage should be Tenacibaculum"
                    " maritimum instead",
                ),
            ),
        ),
        # E coli strain
        (
            dict(
                nnme="E coli strain",
                taxah={"strain": "Escherichia coli 1-110-08_S1_C1"},
            ),
            (
                (
                    WARNING,
                    "No backbone ranks provided in E coli strain's taxa hierarchy",
                ),
                (WARNING, "E coli strain of non-backbone rank: strain"),
            ),
        ),
        # Unknown E coli strain
        (
            dict(
                nnme="Unknown strain",
                taxah={"species": "Escherichia coli", "strain": "Nonsense strain"},
            ),
            (
                (
                    WARNING,
                    "Nonsense strain not registered with NCBI, but higher level taxon "
                    "Escherichia coli is",
                ),
            ),
        ),
        # Ambiguous species
        (
            dict(
                nnme="Ambiguous taxa",
                taxah={"genus": "Morus", "species": "Unknown species"},
            ),
            (
                (
                    ERROR,
                    "Taxa Ambiguous taxa cannot be found and its higher taxonomic "
                    "hierarchy is ambiguous",
                ),
            ),
        ),
        # Nonsense taxonomy
        (
            dict(
                nnme="Utter nonsense",
                taxah={"species": "Nonsense species", "strain": "Nonsense strain"},
            ),
            (
                (
                    ERROR,
                    "Taxa Utter nonsense cannot be found and neither can its higher "
                    "taxonomic hierarchy",
                ),
            ),
        ),
        # No higher taxonomy provided
        (
            dict(nnme="Nonsense", taxah={"species": "Nonsense species"}),
            (
                (
                    ERROR,
                    "Taxa Nonsense cannot be found and its higher taxonomic hierarchy "
                    "is absent",
                ),
            ),
        ),
    ],
)
def test_validate_taxa_search(
    caplog, test_input, expected_log_entries, fixture_ncbi_validators
):
    """This test checks that the function that searches the NCBI database to find
    information on particular taxa logs the correct errors and warnings.
    """

    fixture_ncbi_validators.taxa_search(**test_input)

    log_check(caplog, expected_log_entries)


# Third function that checks that taxa_search throws the appropriate errors
@pytest.mark.parametrize(
    "test_input,expected_exception",
    [
        (None, TypeError),  # no parameters
        (dict(nnme="E coli", taxah=27), TypeError),  # integer instead of dictionary
        (
            dict(
                nnme=27, taxah={"species": "Escherichia coli"}
            ),  # named using a number
            TypeError,
        ),
        (
            dict(
                nnme="E coli", taxah={27: "Escherichia coli"}
            ),  # dictionary key integer
            ValueError,
        ),
        (
            dict(nnme="E coli", taxah={"species": 27}),  # dictionary value not string
            ValueError,
        ),
    ],
)
def test_taxa_search_errors(fixture_ncbi_validators, test_input, expected_exception):
    """This test checks validator.taxa_search inputs throw errors as expected"""

    with pytest.raises(expected_exception):
        _ = fixture_ncbi_validators.taxa_search(**test_input)


# ------------------------------------------
# Testing NCBITaxa
# ------------------------------------------

# Start with the validate_and_add_taxon function

# First check that expected output is recovered
@pytest.mark.parametrize(
    "test_input,expected",
    # Basic case to begin with
    [
        (
            ["E coli", {"genus": "Escherichia", "species": "Escherichia coli"}, 562],
            (
                1,
                1,
                7,
                "E coli",
                ["E coli"],
                [562],
                [561],
                ["Escherichia coli"],
                ["species"],
                ["accepted"],
            ),
        ),
        # Incorrect genus but should be found fine
        (
            ["E coli", {"genus": "Escheria", "species": "Escherichia coli"}, 562],
            (
                1,
                1,
                7,
                "E coli",
                ["E coli"],
                [562],
                [561],
                ["Escherichia coli"],
                ["species"],
                ["accepted"],
            ),
        ),
        # Superseded species name
        (
            ["C vulpes", {"genus": "Canis", "species": "Canis vulpes"}, None],
            (
                2,
                1,
                8,
                "C vulpes",
                ["C vulpes", "C vulpes"],
                [9627, 9627],
                [9625, 9625],
                ["Canis vulpes", "Vulpes vulpes"],
                ["species", "species"],
                ["merged", "accepted"],
            ),
        ),
        # Superseded species name + ID
        (
            ["C marina", {"species": "Cytophaga marina"}, 1000],
            (
                2,
                1,
                7,
                "C marina",
                ["C marina", "C marina"],
                [1000, 107401],
                [104267, 104267],
                ["Cytophaga marina", "Tenacibaculum maritimum"],
                ["species", "species"],
                ["merged", "accepted"],
            ),
        ),
        # Superseded ID
        (
            ["T maritimum", {"species": "Tenacibaculum maritimum"}, 1000],
            (
                2,
                1,
                7,
                "T maritimum",
                ["T maritimum", "T maritimum"],
                [1000, 107401],
                [104267, 104267],
                ["Tenacibaculum maritimum", "Tenacibaculum maritimum"],
                ["species", "species"],
                ["merged", "accepted"],
            ),
        ),
        # Superseded name + correct ID
        (
            ["C marina", {"species": "Cytophaga marina"}, 107401],
            (
                2,
                1,
                7,
                "C marina",
                ["C marina", "C marina"],
                [107401, 107401],
                [104267, 104267],
                ["Cytophaga marina", "Tenacibaculum maritimum"],
                ["species", "species"],
                ["merged", "accepted"],
            ),
        ),
        # Bacteria
        (
            ["Bacteria", {"kingdom": "Bacteria"}, 2],
            (
                1,
                1,
                1,
                "Bacteria",
                ["Bacteria"],
                [2],
                [None],
                ["Bacteria"],
                ["superkingdom"],
                ["accepted"],
            ),
        ),
        # Eukaryota
        (
            ["Eukaryotes", {"superkingdom": "Eukaryota"}, 2759],
            (
                1,
                1,
                1,
                "Eukaryotes",
                ["Eukaryotes"],
                [2759],
                [None],
                ["Eukaryota"],
                ["superkingdom"],
                ["accepted"],
            ),
        ),
        # Fungi
        (
            ["Fungi", {"superkingdom": "Eukaryota", "kingdom": "Fungi"}, 4751],
            (
                1,
                1,
                2,
                "Fungi",
                ["Fungi"],
                [4751],
                [2759],
                ["Fungi"],
                ["kingdom"],
                ["accepted"],
            ),
        ),
        # Unknown strain
        (
            [
                "Unknown strain",
                {"species": "Escherichia coli", "strain": "NBAvgdft"},
                None,
            ],
            (
                1,
                1,
                7,
                "Unknown strain",
                ["Unknown strain"],
                [-1],
                [562],
                ["NBAvgdft"],
                ["strain"],
                ["user"],
            ),
        ),
    ],
)
def test_validate_and_add_taxon(ncbi_resources_local_and_remote, test_input, expected):
    """This test checks the function to validate a taxon against the NCBI
    database actually stores the expected information.
    """

    ncbi_instance = taxa.NCBITaxa(ncbi_resources_local_and_remote)
    ncbi_instance.validate_and_add_taxon(test_input)

    assert len(ncbi_instance.taxon_index) == expected[0]  # Number of taxa added
    assert len(ncbi_instance.taxon_names) == expected[1]  # Number of taxon names
    assert len(ncbi_instance.hierarchy) == expected[2]  # Size of hierarchy

    # Check that provided taxon name is used
    assert list(ncbi_instance.taxon_names)[0] == expected[3]
    # Check that taxon info recorded is as expected
    assert [item[0] for item in ncbi_instance.taxon_index] == expected[4]
    assert [item[1] for item in ncbi_instance.taxon_index] == expected[5]
    assert [item[2] for item in ncbi_instance.taxon_index] == expected[6]
    assert [item[3] for item in ncbi_instance.taxon_index] == expected[7]
    assert [item[4] for item in ncbi_instance.taxon_index] == expected[8]
    assert [item[5] for item in ncbi_instance.taxon_index] == expected[9]


# Now test that the search function logs errors correctly
@pytest.mark.parametrize(
    "test_input,expected_log_entries",
    # Fine so empty
    [
        (
            ["E coli", {"species": "Escherichia coli"}, None],
            ((INFO, "Taxon (E coli) found in NCBI database"),),
        ),
        # Same but with valid code provided
        (
            ["E coli", {"species": "Escherichia coli"}, 562],
            ((INFO, "Taxon (E coli) found in NCBI database"),),
        ),
        # whitespace padding error
        (
            [" E coli", {"species": "Escherichia coli"}, None],
            (
                (ERROR, "Worksheet name has whitespace padding: ' E coli'"),
                (INFO, "Taxon (E coli) found in NCBI database"),
            ),
        ),
        # String of just whitespace provided as name
        (
            [" ", {"species": "Escherichia coli"}, None],
            ((ERROR, "Worksheet name missing, whitespace only or not text"),),
        ),
        # No name error
        (
            [None, {"species": "Escherichia coli"}, None],
            ((ERROR, "Worksheet name missing, whitespace only or not text"),),
        ),
        # Blank string provided as name
        (
            ["", {"species": "Escherichia coli"}, None],
            ((ERROR, "Worksheet name missing, whitespace only or not text"),),
        ),
        # Floats that can be converted to integers are allowed
        (
            ["E coli", {"species": "Escherichia coli"}, 562.0],
            ((INFO, "Taxon (E coli) found in NCBI database"),),
        ),
        # A true float results in multiple errors
        (
            ["E coli", {"species": "Escherichia coli"}, 562.5],
            (
                (ERROR, "NCBI ID contains value that is not an integer"),
                (ERROR, "Improper NCBI ID provided, cannot be validated"),
            ),
        ),
        # As does a string
        (
            ["E coli", {"species": "Escherichia coli"}, "ID"],
            (
                (ERROR, "NCBI ID contains value that is not an integer"),
                (ERROR, "Improper NCBI ID provided, cannot be validated"),
            ),
        ),
        # This checks that multiple errors can fire at once
        (
            ["E coli", {}, 562.5],
            (
                (ERROR, "NCBI ID contains value that is not an integer"),
                (ERROR, "Taxa hierarchy should be a (not empty) dictionary"),
                (ERROR, "Improper NCBI ID provided, cannot be validated"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Taxa hierarchy provided as a string
        (
            ["E coli", "Escherichia coli", None],
            (
                (ERROR, "Taxa hierarchy should be a (not empty) dictionary"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Taxa hierarchy dictionary empty
        (
            ["E coli", {}, None],
            (
                (ERROR, "Taxa hierarchy should be a (not empty) dictionary"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Example of multiple errors
        (
            [" E coli", {}, None],
            (
                (ERROR, "Worksheet name has whitespace padding: ' E coli'"),
                (ERROR, "Taxa hierarchy should be a (not empty) dictionary"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Missing dictionary key
        (
            ["E coli", {" ": "Escherichia coli"}, None],
            (
                (ERROR, "Empty dictionary key used"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Non-string dictionary key
        (
            ["E coli", {26: "Escherichia coli"}, None],
            (
                (ERROR, "Non-string dictionary key used: 26"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Padded dictionary key
        (
            ["E coli", {" species": "Escherichia coli"}, None],
            (
                (ERROR, "Dictionary key has whitespace padding: ' species'"),
                (INFO, "Taxon (E coli) found in NCBI database"),
            ),
        ),
        # Missing dictionary value
        (
            ["E coli", {"species": ""}, None],
            (
                (ERROR, "Empty dictionary value used"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Non-string dictionary value
        (
            ["E coli", {"species": 26}, None],
            (
                (ERROR, "Non-string dictionary value used: 26"),
                (ERROR, "Taxon details not properly formatted, cannot validate"),
            ),
        ),
        # Padding on dictionary value
        (
            ["E coli", {"species": " Escherichia coli"}, None],
            (
                (ERROR, "Dictionary value has whitespace padding: ' Escherichia coli'"),
                (INFO, "Taxon (E coli) found in NCBI database"),
            ),
        ),
        # Taxon hierarchy in wrong order
        (
            ["E coli", {"species": "Escherichia coli", "genus": "Escherichia"}, None],
            ((ERROR, "Taxon hierarchy not in correct order"),),
        ),
        # Right order but non-backbone rank
        (
            [
                "Strepto",
                {"phylum": "Streptophyta", "subphylum": "Streptophytina"},
                None,
            ],
            (
                (WARNING, "Strepto of non-backbone rank: subphylum"),
                (INFO, "Taxon (Strepto) found in NCBI database"),
            ),
        ),
        # Nonsense taxon provided
        (
            ["N garbage", {"species": "Nonsense garbage"}, None],
            (
                (ERROR, "Taxa N garbage cannot be found"),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
        ),
        # Ambiguous taxon provided
        (
            ["Morus", {"genus": "Morus"}, None],
            (
                (
                    ERROR,
                    "Taxa Morus cannot be found using only one taxonomic level, more "
                    "should be provided",
                ),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
        ),
        # Ambiguous taxon resolved
        (
            ["M Morus", {"family": "Moraceae", "genus": "Morus"}, None],
            ((INFO, "Taxon (M Morus) found in NCBI database"),),
        ),
        # E coli with incorrect NCBI ID
        (
            ["E coli", {"species": "Escherichia coli"}, 333],
            (
                (
                    ERROR,
                    "The NCBI ID supplied for E coli does not match hierarchy: expected"
                    " 562 got 333",
                ),
            ),
        ),
        # Non-backbone case with code
        (
            [
                "E coli strain",
                {
                    "species": "Escherichia coli",
                    "strain": "Escherichia coli 1-110-08_S1_C1",
                },
                1444049,
            ],
            (
                (WARNING, "E coli strain of non-backbone rank: strain"),
                (WARNING, "E coli strain of non-backbone rank: strain"),
                (INFO, "Taxon (E coli strain) found in NCBI database"),
            ),
        ),
        # Superseded taxa ID and name
        (
            ["C marina", {"species": "Cytophaga marina"}, 1000],
            (
                (
                    WARNING,
                    "Cytophaga marina not accepted usage should be Tenacibaculum "
                    "maritimum instead",
                ),
                (WARNING, "NCBI ID 1000 has been superseded by ID 107401"),
                (
                    WARNING,
                    "Taxonomic classification superseded for C marina, using new "
                    "taxonomic classification",
                ),
            ),
        ),
        # Just superseded taxa ID
        (
            ["T maritimum", {"species": "Tenacibaculum maritimum"}, 1000],
            (
                (WARNING, "NCBI ID 1000 has been superseded by ID 107401"),
                (
                    WARNING,
                    "NCBI taxa ID superseded for T maritimum, using new taxa ID",
                ),
            ),
        ),
        # Just superseded name
        (
            ["C marina", {"species": "Cytophaga marina"}, 107401],
            (
                (
                    WARNING,
                    "Cytophaga marina not accepted usage should be Tenacibaculum "
                    "maritimum instead",
                ),
                (
                    WARNING,
                    "Taxonomic classification superseded for C marina, using new "
                    "taxonomic classification",
                ),
            ),
        ),
        # E coli recorded as a family rather than a species
        (
            ["E coli", {"family": "Escherichia coli"}, None],
            (
                (ERROR, "Escherichia coli is a species not a family"),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
        ),
        # E coli recorded as a subspecies rather than a species
        (
            ["E coli", {"subspecies": "Escherichia coli"}, None],
            (
                (ERROR, "Escherichia coli is a species not a subspecies"),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
        ),
        # Same idea for a non-backbone case
        (
            ["Streptophytina", {"phylum": "Streptophytina"}, None],
            (
                (WARNING, "Streptophytina of non-backbone rank: subphylum"),
                (ERROR, "Streptophytina is a subphylum not a phylum"),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
        ),
        # Superkingdom not included in GBIF case
        (
            ["Eukaryota", {"superkingdom": "Eukaryota"}, 2759],
            ((INFO, "Taxon (Eukaryota) found in NCBI database"),),
        ),
        # Can actually deal with the bacterial case
        (
            ["Bacteria", {"superkingdom": "Bacteria"}, 2],
            ((INFO, "Taxon (Bacteria) found in NCBI database"),),
        ),
        # Unknown E coli strain
        (
            [
                "Unknown strain",
                {"species": "Escherichia coli", "strain": "Nonsense strain"},
                None,
            ],
            (
                (
                    WARNING,
                    "Nonsense strain not registered with NCBI, but higher level taxon "
                    "Escherichia coli is",
                ),
                (INFO, "Higher taxon for (Unknown strain) resolved in NCBI"),
            ),
        ),
        # Nonsense taxonomy
        (
            [
                "Utter nonsense",
                {"species": "Nonsense species", "strain": "Nonsense strain"},
                None,
            ],
            (
                (
                    ERROR,
                    "Taxa Utter nonsense cannot be found and neither can its higher "
                    "taxonomic hierarchy",
                ),
                (ERROR, "Search based on taxon hierarchy failed"),
            ),
        ),
        # Valid species name, but incorrect genus
        (
            ["E coli", {"genus": "Escheria", "species": "Escherichia coli"}, None],
            (
                (
                    WARNING,
                    "Hierarchy mismatch for E coli its genus should be Escherichia not"
                    " Escheria",
                ),
                (INFO, "Taxon (E coli) found in NCBI database"),
            ),
        ),
        # Valid non-backbone rank but higher taxa wrong
        (
            ["Strepto", {"phylum": "Strephyta", "subphylum": "Streptophytina"}, None],
            (
                (WARNING, "Strepto of non-backbone rank: subphylum"),
                (
                    WARNING,
                    "Hierarchy mismatch for Strepto its phylum should be Streptophyta "
                    "not Strephyta",
                ),
                (INFO, "Taxon (Strepto) found in NCBI database"),
            ),
        ),
    ],
)
def test_validate_and_add_taxon_validate(
    caplog, test_input, expected_log_entries, ncbi_resources_local_and_remote
):
    """This test checks the function to validate a taxon against the NCBI
    database logs the correct information.
    """

    test_input = copy.deepcopy(test_input)
    ncbi_instance = taxa.NCBITaxa(ncbi_resources_local_and_remote)
    ncbi_instance.validate_and_add_taxon(test_input)

    log_check(caplog, expected_log_entries)


# Third function that checks that validate_and_add_taxon throws the appropriate errors
@pytest.mark.parametrize(
    "test_input,expected_exception",
    [
        (None, TypeError),  # no parameters
        ([], ValueError),  # empty list
        (["taxon name", 252], ValueError),  # too few elements
        # Bad code
        (["E coli", {"species": "Escherichia coli"}, -1], ValueError),
        # Bad code
        (["E coli", {"species": "Escherichia coli"}, 100000000000000], taxa.NCBIError),
    ],
)
def test_validate_and_add_taxon_errors(
    ncbi_resources_local_and_remote, test_input, expected_exception
):
    """This test checks validator.validate_and_add_taxon inputs throw errors as
    expected"""

    ncbi_instance = taxa.NCBITaxa(ncbi_resources_local_and_remote)

    with pytest.raises(expected_exception):
        _ = ncbi_instance.validate_and_add_taxon(test_input)


# Then do tests on the index_higher_taxa function

# First test whether sensible output is produced
@pytest.mark.parametrize(
    "test_input,expected",
    # Basic case to begin with
    [
        (
            ["E coli", {"genus": "Escherichia", "species": "Escherichia coli"}, 562],
            (
                7,
                1,
                7,
                "E coli",
                ["E coli", None, None, None, None, None, None],
                [562, 2, 1224, 1236, 91347, 543, 561],
                [561, None, 2, 1224, 1236, 91347, 543],
                [
                    "Escherichia coli",
                    "Bacteria",
                    "Proteobacteria",
                    "Gammaproteobacteria",
                    "Enterobacterales",
                    "Enterobacteriaceae",
                    "Escherichia",
                ],
                [
                    "species",
                    "superkingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                ],
                [
                    "accepted",
                    "accepted",
                    "accepted",
                    "accepted",
                    "accepted",
                    "accepted",
                    "accepted",
                ],
            ),
        ),
        # Superseded ID
        (
            ["T maritimum", {"species": "Tenacibaculum maritimum"}, 1000],
            (
                8,
                1,
                7,
                "T maritimum",
                ["T maritimum", "T maritimum", None, None, None, None, None, None],
                [1000, 107401, 2, 976, 117743, 200644, 49546, 104267],
                [104267, 104267, None, 2, 976, 117743, 200644, 49546],
                [
                    "Tenacibaculum maritimum",
                    "Tenacibaculum maritimum",
                    "Bacteria",
                    "Bacteroidetes",
                    "Flavobacteriia",
                    "Flavobacteriales",
                    "Flavobacteriaceae",
                    "Tenacibaculum",
                ],
                [
                    "species",
                    "species",
                    "superkingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                ],
                [
                    "merged",
                    "accepted",
                    "accepted",
                    "accepted",
                    "accepted",
                    "accepted",
                    "accepted",
                    "accepted",
                ],
            ),
        ),
        # Bacteria
        (
            ["Bacteria", {"kingdom": "Bacteria"}, 2],
            (
                1,
                1,
                1,
                "Bacteria",
                ["Bacteria"],
                [2],
                [None],
                ["Bacteria"],
                ["superkingdom"],
                ["accepted"],
            ),
        ),
        # Eukaryota
        (
            ["Eukaryotes", {"superkingdom": "Eukaryota"}, 2759],
            (
                1,
                1,
                1,
                "Eukaryotes",
                ["Eukaryotes"],
                [2759],
                [None],
                ["Eukaryota"],
                ["superkingdom"],
                ["accepted"],
            ),
        ),
        # Fungi
        (
            ["Fungi", {"superkingdom": "Eukaryota", "kingdom": "Fungi"}, 4751],
            (
                2,
                1,
                2,
                "Fungi",
                ["Fungi", None],
                [4751, 2759],
                [2759, None],
                ["Fungi", "Eukaryota"],
                ["kingdom", "superkingdom"],
                ["accepted", "accepted"],
            ),
        ),
    ],
)
def test_index_higher_taxa(ncbi_resources_local_and_remote, test_input, expected):
    """This test checks the function to store higher taxonomic ranks for a taxon
    actually stores the correct information.
    """

    ncbi_instance = taxa.NCBITaxa(ncbi_resources_local_and_remote)
    ncbi_instance.validate_and_add_taxon(test_input)

    # Then index higher taxa
    ncbi_instance.index_higher_taxa()

    assert len(ncbi_instance.taxon_index) == expected[0]  # Number of taxa added
    assert len(ncbi_instance.taxon_names) == expected[1]  # Number of taxon names
    assert len(ncbi_instance.hierarchy) == expected[2]  # Size of hierarchy

    # Check that provided taxon name is used
    assert list(ncbi_instance.taxon_names)[0] == expected[3]
    # Check that taxon info recorded is as expected
    assert [item[0] for item in ncbi_instance.taxon_index] == expected[4]
    assert [item[1] for item in ncbi_instance.taxon_index] == expected[5]
    assert [item[2] for item in ncbi_instance.taxon_index] == expected[6]
    assert [item[3] for item in ncbi_instance.taxon_index] == expected[7]
    assert [item[4] for item in ncbi_instance.taxon_index] == expected[8]
    assert [item[5] for item in ncbi_instance.taxon_index] == expected[9]


# Now test that the index hierarchy function logs correctly
@pytest.mark.parametrize(
    "test_input,expected_log_entries",
    # Test that E coli works fine
    [
        (
            ["E coli", {"species": "Escherichia coli"}, 562],
            (
                (INFO, "Taxon (E coli) found in NCBI database"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Proteobacteria"),
                (INFO, "Added class Gammaproteobacteria"),
                (INFO, "Added order Enterobacterales"),
                (INFO, "Added family Enterobacteriaceae"),
                (INFO, "Added genus Escherichia"),
            ),
        ),
        # Superseded taxon name used
        (
            ["C marina", {"species": "Cytophaga marina"}, None],
            (
                (
                    WARNING,
                    "Cytophaga marina not accepted usage should be Tenacibaculum "
                    "maritimum instead",
                ),
                (
                    WARNING,
                    "Taxonomic classification superseded for C marina, using new "
                    "taxonomic classification",
                ),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Bacteroidetes"),
                (INFO, "Added class Flavobacteriia"),
                (INFO, "Added order Flavobacteriales"),
                (INFO, "Added family Flavobacteriaceae"),
                (INFO, "Added genus Tenacibaculum"),
            ),
        ),
        # Ambiguous taxon resolved
        (
            ["M Morus", {"family": "Moraceae", "genus": "Morus"}, None],
            (
                (INFO, "Taxon (M Morus) found in NCBI database"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Eukaryota"),
                (INFO, "Added kingdom Viridiplantae"),
                (INFO, "Added phylum Streptophyta"),
                (INFO, "Added class Magnoliopsida"),
                (INFO, "Added order Rosales"),
                (INFO, "Added family Moraceae"),
            ),
        ),
        # Non-backbone case with code
        (
            [
                "E coli strain",
                {
                    "species": "Escherichia coli",
                    "strain": "Escherichia coli 1-110-08_S1_C1",
                },
                1444049,
            ],
            (
                (WARNING, "E coli strain of non-backbone rank: strain"),
                (WARNING, "E coli strain of non-backbone rank: strain"),
                (INFO, "Taxon (E coli strain) found in NCBI database"),
                (INFO, "Indexing taxonomic hierarchy"),
                (INFO, "Added superkingdom Bacteria"),
                (INFO, "Added phylum Proteobacteria"),
                (INFO, "Added class Gammaproteobacteria"),
                (INFO, "Added order Enterobacterales"),
                (INFO, "Added family Enterobacteriaceae"),
                (INFO, "Added genus Escherichia"),
                (INFO, "Added species Escherichia coli"),
            ),
        ),
    ],
)
def test_validate_index_higher_taxa(
    caplog, ncbi_resources_local_and_remote, test_input, expected_log_entries
):
    """This test checks the function to store higher taxonomic ranks for a taxon
    logs the correct information.
    """

    ncbi_instance = taxa.NCBITaxa(ncbi_resources_local_and_remote)
    ncbi_instance.validate_and_add_taxon(test_input)

    # Then index higher taxa
    ncbi_instance.index_higher_taxa()

    log_check(caplog, expected_log_entries)


# Finally check load function (starting by doing an error overview)
@pytest.mark.parametrize(
    "example_ncbi_files, n_errors, n_taxa, t_taxa",
    [("good", 0, 10, 30), ("weird", 0, 5, 20), ("bad", 12, 5, 0)],
    indirect=["example_ncbi_files"],  # take actual params from fixture
)
def test_taxa_load(
    ncbi_resources_local_and_remote, example_ncbi_files, n_errors, n_taxa, t_taxa
):
    """This tests the ensemble loading of (ncbi) taxa from a file using
    indirect parametrisation to access the fixtures containing the
    sample excel files.
    """

    tx = taxa.NCBITaxa(ncbi_resources_local_and_remote)
    tx.load(example_ncbi_files["NCBITaxa"])

    assert tx.n_errors == n_errors
    # Compare both named taxa and total taxa
    assert len(tx.taxon_names) == n_taxa
    assert len(tx.taxon_index) == t_taxa

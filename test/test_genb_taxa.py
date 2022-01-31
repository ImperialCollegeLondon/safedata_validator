import pytest
from safedata_validator import genb_taxa

@pytest.fixture(scope='module')
def validators():
    """Parameterised fixture to return local GBIF validator
    """
    return genb_taxa.RemoteNCBIValidator()

# ------------------------------------------
# Testing NCBITaxon
# ------------------------------------------

@pytest.mark.parametrize(
    'test_input,expected_exception',
    [(dict(),  # no parameters
      TypeError),
     (dict(name=1, genbank_id=2, taxa_hier={'genus': 'Morus'}),  # non string name
      TypeError),
     # non-numeric genbank_id
     (dict(name='Bombus bombus', genbank_id='error', taxa_hier={'genus': 'Morus'}),
      TypeError),
     # non-integer genbank_id
     (dict(name='Bombus bombus', genbank_id=3.141, taxa_hier={'genus': 'Morus'}),
      ValueError),
     # string instead of dictonary
     (dict(name='Bombus bombus', genbank_id=3, taxa_hier="test"),
      TypeError),
     (dict(name='Bombus bombus', genbank_id=3, taxa_hier={}), # empty dictonary
      ValueError),
     (dict(name='Bombus bombus', genbank_id=3, taxa_hier={1: 'Morus'}), # non-string key
      ValueError),
     (dict(name='Bombus bombus', genbank_id=3, taxa_hier={'genus': 27}), # non-string value
      ValueError)])
def test_taxon_init_errors(test_input, expected_exception):
    """This test checks NCBI taxon inputs throw errors as expected
    """
    with pytest.raises(expected_exception):
        _ = genb_taxa.NCBITaxon(**test_input)

# TESTS TO PUT IN:
# Somewhere need a non-backbone taxa rank check
# HOW SHOULD I ALSO TEST THAT THE OUTPUT OF THE LOGGER IS CORRECT

# ------------------------------------------
# Testing taxon validators
# ------------------------------------------

@pytest.mark.parametrize(
    'test_input,expected',
    [(dict(nnme='E coli', genbank_id=562),
      ('E coli', 562, None, False, "species", "Escherichia coli")),
     (dict(nnme='E coli strain', genbank_id=1444049),
      ('E coli strain', 1444049, "strain", False, "species", "Escherichia coli")),
     (dict(nnme='Streptophytina', genbank_id=131221),
      ('Streptophytina', 131221, "subphylum", False, "phylum", "Streptophyta")),
     (dict(nnme='Opisthokonta', genbank_id=33154),
      ('Opisthokonta', 33154, "clade", False, "superkingdom", "Eukaryota")),
     (dict(nnme='Cytophaga marina', genbank_id=1000),
      ('Cytophaga marina', 1000, None, True, "species", "Tenacibaculum maritimum"))
     ])
def test_id_lookup(validators, test_input, expected):
    """This test checks the results of looking up a specific NCBI taxonomy against
    the actual returned taxonomy details. At the moment this for the Remote validator,
    but if a local validator is defined this should also be checked against.

    """

    fnd_tx = validators.id_lookup(**test_input)

    assert fnd_tx.name == expected[0]
    assert fnd_tx.genbank_id == expected[1]
    assert fnd_tx.diverg == expected[2]
    assert fnd_tx.superseed == expected[3]

    # Find last dictonary key
    f_key = list(fnd_tx.taxa_hier.keys())[-1]

    assert f_key == expected[4]
    assert fnd_tx.taxa_hier[f_key] == expected[5]

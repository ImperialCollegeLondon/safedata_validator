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
# HOW SHOULD I ALSO TEST THAT THE OUTPUT OF THE LOGGER IS CORRECT

# ------------------------------------------
# Testing taxon validators
# ------------------------------------------

# First test the search function
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

# Then do the same for the taxa serach function
@pytest.mark.parametrize(
    'test_input,expected',
    [(dict(nnme="E coli",taxa={'genus': 'Escherichia', 'species': 'Escherichia coli'}),
      ("E coli", 562, None, False, "species", "Escherichia coli")),
     (dict(nnme="Entero",taxa={'order': 'Enterobacterales', 'family': 'Enterobacteriaceae'}),
      ("Entero", 543, None, False, "family", "Enterobacteriaceae")),
     (dict(nnme="E coli strain",taxa={'species': 'Escherichia coli', 'strain':
     'Escherichia coli 1-110-08_S1_C1'}),
      ("E coli strain", 1444049, "strain", False, "species", "Escherichia coli")),
     (dict(nnme="Strepto",taxa={'phylum': 'Streptophyta', 'subphylum': 'Streptophytina'}),
      ("Strepto", 131221, "subphylum", False, "phylum", "Streptophyta")),
     (dict(nnme="Opistho",taxa={'superkingdom': 'Eukaryota', 'clade': 'Opisthokonta'}),
      ("Opistho", 33154, "clade", False, "superkingdom", "Eukaryota")),
     (dict(nnme="Vulpes",taxa={'genus': 'Vulpes', 'species': 'Vulpes vulpes'}),
      ("Vulpes", 9627, None, False, "species", "Vulpes vulpes")),
     (dict(nnme="M morus",taxa={'family': 'Moraceae', 'genus': 'Morus'}),
      ("M morus", 3497, None, False, "genus", "Morus")),
     (dict(nnme="S morus",taxa={'family': 'Sulidae', 'genus': 'Morus'}),
      ("S morus", 37577, None, False, "genus", "Morus")),
     (dict(nnme="C morus",taxa={'phylum': 'Chordata', 'genus': 'Morus'}),
      ("C morus", 37577, None, False, "genus", "Morus")),
     (dict(nnme="T maritimum",taxa={'genus': 'Tenacibaculum', 'species':
     'Tenacibaculum maritimum'}),
      ("T maritimum", 107401, None, False, "species", "Tenacibaculum maritimum")),
     (dict(nnme="C marina",taxa={'genus': 'Cytophaga', 'species': 'Cytophaga marina'}),
      ("C marina", 107401, None, True, "species", "Tenacibaculum maritimum"))
    ])
def test_taxa_search(validators, test_input, expected):
    """This test checks the results of searching for a taxa in the NCBI taxonomy
    database against what was expected to be found. At the moment this for the
    Remote validator, but if a local validator is defined this should also be
    checked against.

    """
    fnd_tx = validators.taxa_search(**test_input)

    assert fnd_tx.name == expected[0]
    assert fnd_tx.genbank_id == expected[1]
    assert fnd_tx.diverg == expected[2]
    assert fnd_tx.superseed == expected[3]

    # Find last dictonary key
    f_key = list(fnd_tx.taxa_hier.keys())[-1]

    assert f_key == expected[4]
    assert fnd_tx.taxa_hier[f_key] == expected[5]

# THESE GIVE ERRORS SO MUST BE HANDLED SEPERATLY
# # Morus (NA)
# d7 = {'genus': 'Morus'}
# # Eukaryota morus (NA)
# d11 = {'superkingdom': 'Eukaryota', 'genus': 'Morus'}
# # Nonsense garbage (NA)
# d13 = {'genus': 'Nonsense', 'species': 'Nonsense garbage'}
# d16 => NEITHER BACKBONE RANK => DOES THIS ACTUALLY NEED ITS OWN ERROR

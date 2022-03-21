import pytest
from safedata_validator import genb_taxa
from logging import ERROR, WARNING, INFO
from .conftest import log_check

# TESTS TO PUT IN:
# NEED TO CHECK THAT GENBANKTAXA INSTANCES ARE INITIALISED CORRECTLY
# UNIT TESTS FOR LOAD FUNCTIONS => ADD THIS WHEN I'VE FIXED VALIDATE ETC
# CAN I TEST INDEX HIGHER TAXA???

# MOVE THIS TO CONFTEST EVENTUALLY
@pytest.fixture(scope='module')
def gb_instance():
    """Parameterised fixture to return a GenBankTaxa instance
    """
    return genb_taxa.GenBankTaxa()

# ------------------------------------------
# Testing NCBITaxon
# ------------------------------------------

@pytest.mark.parametrize(
    'test_input,expected_exception',
    [(dict(),  # no parameters
      TypeError),
     (dict(name=1, ncbi_id=2, taxa_hier={'genus': ('Morus', 37577)}),  # non string name
      TypeError),
     # non-numeric ncbi_id
     (dict(name='Bombus bombus', ncbi_id='error',
      taxa_hier={'genus': ('Morus', 37577)}),
      TypeError),
     # non-integer ncbi_id
     (dict(name='Bombus bombus', ncbi_id=3.141, taxa_hier={'genus': ('Morus', 37577)}),
      TypeError),
     # string instead of dictonary
     (dict(name='Bombus bombus', ncbi_id=3, taxa_hier="test"),
      TypeError),
     (dict(name='Bombus bombus', ncbi_id=3, taxa_hier={}), # empty dictonary
      ValueError),
     # non-string key
     (dict(name='Bombus bombus', ncbi_id=3, taxa_hier={1: ('Morus', 37577)}),
      ValueError),
     (dict(name='Bombus bombus', ncbi_id=3, taxa_hier={'genus': 27}), # non-tuple value
      ValueError),
     # non-tuple value
     (dict(name='Bombus bombus', ncbi_id=3, taxa_hier={'genus': (37577, 'Morus')}),
      ValueError),
     # 3 elements in tuple
     (dict(name='Bombus bombus', ncbi_id=3,
      taxa_hier={'genus': ('Morus', 37577, "extra")}),
      ValueError)])
def test_taxon_init_errors(test_input, expected_exception):
    """This test checks NCBI taxon inputs throw errors as expected
    """
    with pytest.raises(expected_exception):
        _ = genb_taxa.NCBITaxon(**test_input)

# ------------------------------------------
# Testing taxa_strip
# ------------------------------------------
# Only need to test that output is sensible here
@pytest.mark.parametrize(
    'test_input,expected',
    [(dict(name='Bacteria', rank='Kingdom'),
      ('Bacteria', True)),
     (dict(name='k__Bacteria', rank='Kingdom'),
      ('Bacteria', True)),
     (dict(name='k__Bacteria', rank='Phylum'),
      ('Bacteria', False)),
     (dict(name='p__Acidobacteria', rank='Phylum'),
      ('Acidobacteria', True)),
     (dict(name='s__', rank='Species'),
      ('', True)),
     (dict(name='s__', rank='Species'),
      ('', True)),
     (dict(name=None, rank='Order'),
      (None, True)),
     ])
def test_taxa_strip(test_input, expected):
    """This test checks the function that strips taxa strings down to remove k__
    type notation is functioning properly. This function also checks that the
    supplied rank matches the rank implied by the notation, we also test this
    behaviour.
    """

    s_taxa, match = genb_taxa.taxa_strip(**test_input)

    assert s_taxa == expected[0]
    assert match == expected[1]

# ------------------------------------------
# Testing species_binomial
# ------------------------------------------
@pytest.mark.parametrize(
    'test_input,expected',
    [(dict(genus='Escherichia', species='coli'),
      'Escherichia coli'),
     (dict(genus='Escherichia', species='Escherichia coli'),
       'Escherichia coli'),
     (dict(genus='Gorilla', species='gorilla'),
      'Gorilla gorilla'),
     (dict(genus='Candidatus Koribacter', species='Candidatus versatilis'),
      'Candidatus Koribacter versatilis'),
     (dict(genus='Candidatus Koribacter', species='versatilis'),
      'Candidatus Koribacter versatilis'),
     (dict(genus='Over long genus name', species='vulpes'),
      None),
     (dict(genus='Canis', species='Vulpes vulpes'),
      None),
     ])
def test_species_binomial(test_input, expected):
    """This test checks the function that strips constructs species binomals from
    species and genus names. We test that it can catch when the species name is
    already a binomial, and that it can catch Candidatus names and handle them
    properly.
    """

    s_bi = genb_taxa.species_binomial(**test_input)

    assert s_bi == expected

# Now test that the function logs errors correctly
@pytest.mark.parametrize(
    'test_input,expected_log_entries',
    [(dict(genus='Escherichia', species='coli'), # Fine so empty
      ()),
     (dict(genus='Over long genus name', species='vulpes'), # Over long name
      ((ERROR, 'Genus name (Over long genus name) appears to be too long'),
      )),
     (dict(genus='Canis', species='Vulpes vulpes'), # Genus name not in binomial
      ((ERROR, 'Species name (Vulpes vulpes) appears to be binomal but does not'
               ' contain genus name (Canis)'),
      ))
     ])
def test_validate_species_binomial(caplog, test_input, expected_log_entries):
    """This test checks that the function to construct species binomials logs the
     correct errors and warnings.
    """

    s_bi = genb_taxa.species_binomial(**test_input)

    log_check(caplog, expected_log_entries)

# ------------------------------------------
# Testing subspecies_trinomial
# ------------------------------------------
@pytest.mark.parametrize(
    'test_input,expected',
    [(dict(species='Vulpes vulpes', subspecies='japonica'),
      'Vulpes vulpes japonica'),
     (dict(species='Candidatus Koribacter versatilis', subspecies='Ellin345'),
       'Candidatus Koribacter versatilis Ellin345'),
     (dict(species='Candidatus Koribacter versatilis', subspecies='Candidatus Ellin345'),
      'Candidatus Koribacter versatilis Ellin345'),
     (dict(species='Vulpes vulpes', subspecies='Vulpes vulpes schrenckii'),
      'Vulpes vulpes schrenckii'),
     (dict(species='Canis vulpes', subspecies='Vulpes vulpes schrenckii'),
      None),
     (dict(species='Over long name', subspecies='schrenckii'),
      None),
     (dict(species='Vulpes', subspecies='Vulpes vulpes schrenckii'),
      None),
     ])
def test_subspecies_trinomial(test_input, expected):
    """This test checks the function that strips constructs species binomals from
    species and genus names. We test that it can catch when the species name is
    already a binomial, and that it can catch Candidatus names and handle them
    properly.
    """

    s_bi = genb_taxa.subspecies_trinomial(**test_input)

    assert s_bi == expected

# Now test that the function logs errors correctly
@pytest.mark.parametrize(
    'test_input,expected_log_entries',
    [(dict(species='Vulpes vulpes', subspecies='Vulpes vulpes japonica'), # Fine so empty
      ()),
     # Species name and genus name don't match
     (dict(species='Canis vulpes', subspecies='Vulpes vulpes schrenckii'),
      ((ERROR, 'Subspecies name (Vulpes vulpes schrenckii) appears to be trinomal'
                   f'but does not contain species name (Canis vulpes)'),
      )),
     (dict(species='Over long name', subspecies='schrenckii'), # Over long name
      ((ERROR, 'Species name (Over long name) too long'),
      )),
     (dict(species='Vulpes', subspecies='Vulpes vulpes schrenckii'), # Too short
      ((ERROR, 'Species name (Vulpes) too short'),
      ))
     ])
def test_validate_subspecies_trinomial(caplog, test_input, expected_log_entries):
    """This test checks that the function to construct species binomials logs the
     correct errors and warnings.
    """

    s_bi = genb_taxa.subspecies_trinomial(**test_input)

    log_check(caplog, expected_log_entries)

# ------------------------------------------
# Testing taxon validators
# ------------------------------------------

# First test the search function
@pytest.mark.parametrize(
    'test_input,expected',
    [(dict(nnme='E coli', ncbi_id=562),
      ('E coli', 562, False, "species", "Escherichia coli", 562)),
     (dict(nnme='E coli strain', ncbi_id=1444049),
      ('E coli strain', 1444049, False, "strain", "Escherichia coli 1-110-08_S1_C1", 1444049)),
     (dict(nnme='Streptophytina', ncbi_id=131221),
      ('Streptophytina', 131221, False, "subphylum", "Streptophytina", 131221)),
     (dict(nnme='Opisthokonta', ncbi_id=33154),
      ('Opisthokonta', 33154, False, "clade", "Opisthokonta", 33154)),
     (dict(nnme='Cytophaga marina', ncbi_id=1000),
      ('Cytophaga marina', 1000, True, "species", "Tenacibaculum maritimum", 107401))
     ])
def test_id_lookup(fixture_ncbi_validators, test_input, expected):
    """This test checks the results of looking up a specific NCBI taxonomy against
    the actual returned taxonomy details. At the moment this for the Remote validator,
    but if a local validator is defined this should also be checked against.

    """

    fnd_tx = fixture_ncbi_validators.id_lookup(**test_input)

    assert fnd_tx.name == expected[0]
    assert fnd_tx.ncbi_id == expected[1]
    assert fnd_tx.superseed == expected[2]

    # Find last dictonary key
    f_key = list(fnd_tx.taxa_hier.keys())[-1]

    assert f_key == expected[3]
    assert fnd_tx.taxa_hier[f_key] == (expected[4], expected[5])

# Now test that the search function logs errors correctly
@pytest.mark.parametrize(
    'test_input,expected_log_entries',
    [(dict(nnme='E coli', ncbi_id=562), # Fine so empty
      ()),
     (dict(nnme='Streptophytina', ncbi_id=131221), # Non-backbone rank
      ((WARNING, 'Streptophytina of non-backbone rank: subphylum'),
      )),
     (dict(nnme='C marina', ncbi_id=1000), # Non-backbone rank
      ((WARNING, 'NCBI ID 1000 has been superseeded by ID 107401'),
      ))
     ])
def test_validate_id_lookup(caplog, test_input, expected_log_entries, fixture_ncbi_validators):
    """This test checks that the function to generate NCBITaxon objects by searching
    an ID logs the correct errors and warnings. At the moment this for the Remote
    validator, but if a local validator is defined this should also be checked against.

    """

    fnd_tx = fixture_ncbi_validators.id_lookup(**test_input)

    log_check(caplog, expected_log_entries)

# Third function that checks that id_lookup throws the appropriate errors
@pytest.mark.parametrize(
    'test_input,expected_exception',
    [(None,  # no parameters
      TypeError),
     (dict(nnme='E coli', ncbi_id='invalid_string'),  # a string
      ValueError),
     (dict(nnme='E coli', ncbi_id=27.5),  # non-integer ID
      ValueError),
     (dict(nnme=27, ncbi_id=27),  # named using a number
      TypeError),
     (dict(nnme='E coli', ncbi_id=-1),  # bad ID
      ValueError),
     (dict(nnme='E coli', ncbi_id=100000000000000),  # bad ID
      genb_taxa.NCBIError)])
def test_id_lookup_errors(fixture_ncbi_validators, test_input, expected_exception):
    """This test checks validator.id_lookup inputs throw errors as expected
    """

    with pytest.raises(expected_exception):
        _ = fixture_ncbi_validators.id_lookup(**test_input)

# Then do the same for the taxa serach function
@pytest.mark.parametrize(
    'test_input,expected',
    [(dict(nnme="E coli",taxah={'genus': 'Escherichia', 'species': 'Escherichia coli'}),
      ("E coli", 562, False, "species", "Escherichia coli", 562)),
     (dict(nnme="Entero",taxah={'order': 'Enterobacterales', 'family': 'Enterobacteriaceae'}),
      ("Entero", 543, False, "family", "Enterobacteriaceae", 543)),
     (dict(nnme="E coli strain",taxah={'species': 'Escherichia coli', 'strain':
     'Escherichia coli 1-110-08_S1_C1'}),
      ("E coli strain", 1444049, False, "strain", "Escherichia coli 1-110-08_S1_C1", 1444049)),
     (dict(nnme="Strepto",taxah={'phylum': 'Streptophyta', 'subphylum': 'Streptophytina'}),
      ("Strepto", 131221, False, "subphylum", "Streptophytina", 131221)),
     (dict(nnme="Opistho",taxah={'superkingdom': 'Eukaryota', 'clade': 'Opisthokonta'}),
      ("Opistho", 33154, False, "clade", "Opisthokonta", 33154)),
     (dict(nnme="Vulpes",taxah={'genus': 'Vulpes', 'species': 'Vulpes vulpes'}),
      ("Vulpes", 9627, False, "species", "Vulpes vulpes", 9627)),
     (dict(nnme="M morus",taxah={'family': 'Moraceae', 'genus': 'Morus'}),
      ("M morus", 3497, False, "genus", "Morus", 3497)),
     (dict(nnme="S morus",taxah={'family': 'Sulidae', 'genus': 'Morus'}),
      ("S morus", 37577, False, "genus", "Morus", 37577)),
     (dict(nnme="C morus",taxah={'phylum': 'Chordata', 'genus': 'Morus'}),
      ("C morus", 37577, False, "genus", "Morus", 37577)),
     (dict(nnme="T maritimum",taxah={'genus': 'Tenacibaculum', 'species':
     'Tenacibaculum maritimum'}),
      ("T maritimum", 107401, False, "species", "Tenacibaculum maritimum", 107401)),
     (dict(nnme="C marina",taxah={'genus': 'Cytophaga', 'species': 'Cytophaga marina'}),
      ("C marina", 107401, True, "species", "Tenacibaculum maritimum", 107401)),
     (dict(nnme="Bacteria",taxah={'superkingdom': 'Bacteria'}),
      ("Bacteria", 2, False, "superkingdom", "Bacteria", 2))
    ])
def test_taxa_search(fixture_ncbi_validators, test_input, expected):
    """This test checks the results of searching for a taxa in the NCBI taxonomy
    database against what was expected to be found. At the moment this for the
    Remote validator, but if a local validator is defined this should also be
    checked against.

    """
    fnd_tx = fixture_ncbi_validators.taxa_search(**test_input)

    assert fnd_tx.name == expected[0]
    assert fnd_tx.ncbi_id == expected[1]
    assert fnd_tx.superseed == expected[2]

    # Find last dictonary key
    f_key = list(fnd_tx.taxa_hier.keys())[-1]

    assert f_key == expected[3]
    assert fnd_tx.taxa_hier[f_key] == (expected[4], expected[5])

# Now test that the search function logs errors correctly
@pytest.mark.parametrize(
    'test_input,expected_log_entries',
    # Fine so empty
    [(dict(nnme='E coli', taxah={'genus': 'Escherichia', 'species': 'Escherichia coli'}),
      ()),
     (dict(nnme='Nonsense', taxah={'species': 'Nonsense garbage'}), # Nonsense
      ((ERROR, 'Taxa Nonsense cannot be found'),
      )),
     (dict(nnme='Morus', taxah={'genus': 'Morus'}), # ambigious Morus
      ((ERROR, 'Taxa Morus cannot be found using only one taxonomic level,'
      ' more should be provided'),
      )),
     # Nonsense Morus
     (dict(nnme='N Morus', taxah={'family': 'Nonsense', 'genus': 'Morus'}),
      ((ERROR, 'Provided parent taxa for N Morus not found'),
      )),
     # Aegle Morus
     (dict(nnme='A Morus', taxah={'family': 'Aegle', 'genus': 'Morus'}),
      ((ERROR, 'More than one possible parent taxa for A Morus found'),
      )),
     # Eukaryote Morus
     (dict(nnme='E Morus', taxah={'superkingdom': 'Eukaryota', 'genus': 'Morus'}),
      ((ERROR, 'Parent taxa for E Morus refers to multiple possible child taxa'),
      )),
     # Carnivora Morus
     (dict(nnme='C Morus', taxah={'superkingdom': 'Carnivora', 'genus': 'Morus'}),
      ((ERROR, 'Parent taxa not actually a valid parent of C Morus'),
      )),
     # Cytophaga marina
     (dict(nnme='C marina', taxah={'genus': 'Cytophaga', 'species': 'Cytophaga marina'}),
      ((WARNING, 'Cytophaga marina not accepted usage should be Tenacibaculum'
      ' maritimum instead'),
      )),
     # E coli strain
     (dict(nnme='E coli strain', taxah={'strain': 'Escherichia coli 1-110-08_S1_C1'}),
      ((WARNING, "No backbone ranks provided in E coli strain's taxa hierarchy"),
       (WARNING, "E coli strain of non-backbone rank: strain"),
      )),
     ])
def test_validate_taxa_search(caplog, test_input, expected_log_entries,
                              fixture_ncbi_validators):
    """This test checks that the function that searches the NCBI database to find
    information on particular taxa logs the correct errors and warnings. At the
    moment this for the Remote validator, but if a local validator is defined
    this should also be checked against.

    """

    fnd_tx = fixture_ncbi_validators.taxa_search(**test_input)

    log_check(caplog, expected_log_entries)

# Third function that checks that taxa_search throws the appropriate errors
@pytest.mark.parametrize(
    'test_input,expected_exception',
    [(None,  # no parameters
      TypeError),
     (dict(nnme='E coli', taxah=27),  # integer instead of dictonary
      TypeError),
     (dict(nnme=27, taxah={'species': 'Escherichia coli'}),  # named using a number
      TypeError),
     (dict(nnme='E coli', taxah={27: 'Escherichia coli'}),  # dictonary key integer
      ValueError),
     (dict(nnme='E coli', taxah={'species': 27}),  # dictonary value not string
      ValueError),
    ])

def test_taxa_search_errors(fixture_ncbi_validators, test_input, expected_exception):
    """This test checks validator.taxa_search inputs throw errors as expected
    """

    with pytest.raises(expected_exception):
        _ = fixture_ncbi_validators.taxa_search(**test_input)

# ------------------------------------------
# Testing GenBankTaxa
# ------------------------------------------

# Start with the validate_and_add_taxon function

# SECTION TESTING EXPECTED OUTPUT

# Now test that the search function logs errors correctly
@pytest.mark.parametrize(
    'test_input,expected_log_entries',
    # Fine so empty
    [(['E coli', {'species': 'Escherichia coli'}, None],
      ((INFO, "Taxon (E coli) found in NCBI database"),
      )),
     # Same but with valid code provided
     (['E coli', {'species': 'Escherichia coli'}, 562],
      ((INFO, "Taxon (E coli) found in NCBI database"),
      )),
     # whitespace padding error
     ([' E coli', {'species': 'Escherichia coli'}, None],
      ((ERROR, "Worksheet name has whitespace padding: ' E coli'"),
       (INFO, "Taxon (E coli) found in NCBI database"),
      )),
     # String of just whitespace provided as name
     ([' ', {'species': 'Escherichia coli'}, None],
      ((ERROR, "Worksheet name missing, whitespace only or not text"),
      )),
     # No name error
     ([None, {'species': 'Escherichia coli'}, None],
      ((ERROR, "Worksheet name missing, whitespace only or not text"),
      )),
     # Blank string provided as name
     (['', {'species': 'Escherichia coli'}, None],
      ((ERROR, "Worksheet name missing, whitespace only or not text"),
      )),
     # Floats that can be cocnverted to integers are allowed
     (['E coli', {'species': 'Escherichia coli'}, 562.0],
      ((INFO, "Taxon (E coli) found in NCBI database"),
      )),
     # A true float results in multiple errors
     (['E coli', {'species': 'Escherichia coli'}, 562.5],
      ((ERROR, "NCBI ID contains value that is not an integer"),
       (ERROR, "Improper NCBI ID provided, cannot be validated"),
      )),
     # As does a string
     (['E coli', {'species': 'Escherichia coli'}, "ID"],
      ((ERROR, "NCBI ID contains value that is not an integer"),
       (ERROR, "Improper NCBI ID provided, cannot be validated"),
      )),
     # This checks that multiple errors can fire at once
     (['E coli', {}, 562.5],
      ((ERROR, "NCBI ID contains value that is not an integer"),
       (ERROR, "Taxa hierarchy should be a (not empty) dictonary"),
       (ERROR, "Improper NCBI ID provided, cannot be validated"),
       (ERROR, "Taxon details not properly formatted, cannot validate"),
      )),
     # Taxa hierarchy provided as a string
     (['E coli', 'Escherichia coli', None],
      ((ERROR, "Taxa hierarchy should be a (not empty) dictonary"),
       (ERROR, "Taxon details not properly formatted, cannot validate"),
      )),
     # Taxa hierarchy dictonary empty
     (['E coli', {}, None],
      ((ERROR, "Taxa hierarchy should be a (not empty) dictonary"),
       (ERROR, "Taxon details not properly formatted, cannot validate"),
      )),
     # Example of multiple errors
     ([' E coli', {}, None],
      ((ERROR, "Worksheet name has whitespace padding: ' E coli'"),
       (ERROR, "Taxa hierarchy should be a (not empty) dictonary"),
       (ERROR, "Taxon details not properly formatted, cannot validate"),
      )),
     # Missing dictonary key
     (['E coli', {' ': 'Escherichia coli'}, None],
      ((ERROR, "Empty dictonary key used"),
       (ERROR, "Taxon details not properly formatted, cannot validate"),
      )),
     # Non-string dictonary key
     (['E coli', {26: 'Escherichia coli'}, None],
      ((ERROR, "Non-string dictonary key used: 26"),
       (ERROR, "Taxon details not properly formatted, cannot validate"),
      )),
     # Padded dictonary key
     (['E coli', {' species': 'Escherichia coli'}, None],
      ((ERROR, "Dictonary key has whitespace padding: ' species'"),
       (INFO, "Taxon (E coli) found in NCBI database"),
      )),
     # Missing dictonary value
     (['E coli', {'species': ''}, None],
      ((ERROR, "Empty dictonary value used"),
       (ERROR, "Taxon details not properly formatted, cannot validate"),
      )),
     # Non-string dictonary value
     (['E coli', {'species': 26}, None],
      ((ERROR, "Non-string dictonary value used: 26"),
       (ERROR, "Taxon details not properly formatted, cannot validate"),
      )),
     # Padding on dictonary value
     (['E coli', {'species': ' Escherichia coli'}, None],
      ((ERROR, "Dictonary value has whitespace padding: ' Escherichia coli'"),
       (INFO, "Taxon (E coli) found in NCBI database"),
      )),
     # Taxon hierarchy in wrong order
     (['E coli', {'species': 'Escherichia coli', 'genus': 'Escherichia'}, None],
      ((ERROR, "Taxon hierarchy not in correct order"),
      )),
     # Right order but non-backbone rank
     (['Strepto', {'phylum': 'Streptophyta', 'subphylum': 'Streptophytina'}, None],
      ((WARNING, "Strepto not of backbone rank, instead resolved to phylum level"),
      )),
     # Superseeded taxon name used
     (['C marina', {'species': 'Cytophaga marina'}, None],
      ((WARNING, "Cytophaga marina not accepted usage should be Tenacibaculum maritimum instead"),
       (WARNING, "Taxonomic classification superseeded for C marina, using new taxonomic classification"),
      )),
     # Nonsense taxon provided
     (['N garbage', {'species': 'Nonsense garbage'}, None],
      ((ERROR, "Taxa N garbage cannot be found"),
       (ERROR, "Search based on taxon hierarchy failed"),
      )),
     # Ambigious taxon provided
     (['Morus', {'genus': 'Morus'}, None],
      ((ERROR, "Taxa Morus cannot be found using only one taxonomic level, more should be provided"),
       (ERROR, "Search based on taxon hierarchy failed"),
      )),
     # Ambigious taxon resolved
     (['M Morus', {'family': 'Moraceae', 'genus': 'Morus'}, None],
      ((INFO, "Taxon (M Morus) found in NCBI database"),
      )),
     # E coli with incorrect NCBI ID
     (['E coli', {'species': 'Escherichia coli'}, 333],
      ((ERROR, "The NCBI ID supplied for E coli does not match hierarchy: expected 562 got 333"),
      )),
     # Non-backbone case with code
     (['E coli strain', {'species': 'Escherichia coli', 'strain': 'Escherichia coli 1-110-08_S1_C1'},
       1444049],
      ((WARNING, "E coli strain not of backbone rank, instead resolved to species level"),
       (WARNING, "E coli strain not of backbone rank, instead resolved to species level"),
      )),
     # Superseeded taxa ID and name
     (['C marina', {'species': 'Cytophaga marina'}, 1000],
      ((WARNING, "Cytophaga marina not accepted usage should be Tenacibaculum maritimum instead"),
       (WARNING, "NCBI ID 1000 has been superseeded by ID 107401"),
       (WARNING, "Taxonomic classification superseeded for C marina, using new taxonomic classification"),
      )),
     # Just superseeded taxa ID
     (['T maritimum', {'species': 'Tenacibaculum maritimum'}, 1000],
      ((WARNING, "NCBI ID 1000 has been superseeded by ID 107401"),
       (WARNING, "NCBI taxa ID superseeded for T maritimum, using new taxa ID"),
      )),
     # Just superseeded name
     (['C marina', {'species': 'Cytophaga marina'}, 107401],
      ((WARNING, "Cytophaga marina not accepted usage should be Tenacibaculum maritimum instead"),
       (WARNING, "Taxonomic classification superseeded for C marina, using new taxonomic classification"),
      )),
     # E coli recorded as a family rather than a species
     (['E coli', {'family': 'Escherichia coli'}, None],
      ((ERROR, "Escherichia coli is a species not a family"),
       (ERROR, "Search based on taxon hierarchy failed"),
      )),
     # E coli recorded as a subspecies rather than a species
     (['E coli', {'subspecies': 'Escherichia coli'}, None],
      ((ERROR, "Escherichia coli is a species not a subspecies"),
       (ERROR, "Search based on taxon hierarchy failed"),
      )),
     # Same idea for a non-backbone case
     (['Streptophytina', {'phylum': 'Streptophytina'}, None],
      ((WARNING, "Streptophytina not of backbone rank, instead resolved to phylum level"),
       (ERROR, "Streptophytina is a subphylum not a phylum"),
       (ERROR, "Search based on taxon hierarchy failed"),
      )),
     # Superkingdom not included in GBIF case
     (['Eukaryota', {'superkingdom': 'Eukaryota'}, 2759],
      ((INFO, "Taxon (Eukaryota) found in NCBI database"),
      )),
     # Can actually deal with the bacterial case
     (['Bacteria', {'superkingdom': 'Bacteria'}, 2],
      ((INFO, "Taxon (Bacteria) found in NCBI database"),
      )),
     ])
def test_validate_and_add_taxon(caplog, test_input, expected_log_entries):
    """This test checks that the function that searches the NCBI database to find
    information on particular taxa logs the correct errors and warnings. At the
    moment this for the Remote validator, but if a local validator is defined
    this should also be checked against.

    """

    # ONLY MAKING A LOCAL VERSION FOR NOW
    gb_instance = genb_taxa.GenBankTaxa()
    fnd_tx = gb_instance.validate_and_add_taxon(test_input)

    log_check(caplog, expected_log_entries)

# Third function that checks that validate_and_add_taxon throws the appropriate errors
@pytest.mark.parametrize(
    'test_input,expected_exception',
    [(None,  # no parameters
      TypeError),
     ([],  # empty list
      ValueError),
     (['taxon name', 252],  # too few elements
      ValueError),
     # Bad code
     (['E coli', {'species': 'Escherichia coli'}, -1],
      ValueError),
     # Bad code
     (['E coli', {'species': 'Escherichia coli'}, 100000000000000],
      genb_taxa.NCBIError),
    ])
def test_validate_and_add_taxon_errors(test_input, expected_exception):
    """This test checks validator.validate_and_add_taxon inputs throw errors as expected
    """

    # ONLY MAKING A LOCAL VERSION FOR NOW
    gb_instance = genb_taxa.GenBankTaxa()

    with pytest.raises(expected_exception):
        _ = gb_instance.validate_and_add_taxon(test_input)

import pytest
from safedata_validator import genb_taxa

# TESTS TO PUT IN:
# NEED TO CHECK THAT GENBANKTAXA INSTANCES ARE INITIALISED CORRECTLY
# NEED TO CHECK THAT REPEATED ENTRIES ARE ACTUALLY SKIPPED

# MOVE THESE TO CONFTEST
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
     (dict(name=1, genbank_id=2, taxa_hier={'genus': 'Morus'}),  # non string name
      TypeError),
     # non-numeric genbank_id
     (dict(name='Bombus bombus', genbank_id='error', taxa_hier={'genus': 'Morus'}),
      TypeError),
     # non-integer genbank_id
     (dict(name='Bombus bombus', genbank_id=3.141, taxa_hier={'genus': 'Morus'}),
      TypeError),
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
def test_id_lookup(fixture_ncbi_validators, test_input, expected):
    """This test checks the results of looking up a specific NCBI taxonomy against
    the actual returned taxonomy details. At the moment this for the Remote validator,
    but if a local validator is defined this should also be checked against.

    """

    fnd_tx = fixture_ncbi_validators.id_lookup(**test_input)

    assert fnd_tx.name == expected[0]
    assert fnd_tx.genbank_id == expected[1]
    assert fnd_tx.diverg == expected[2]
    assert fnd_tx.superseed == expected[3]

    # Find last dictonary key
    f_key = list(fnd_tx.taxa_hier.keys())[-1]

    assert f_key == expected[4]
    assert fnd_tx.taxa_hier[f_key] == expected[5]

# Now test that the search function logs errors correctly
@pytest.mark.parametrize(
    'test_input,expected_log_entries',
    [(dict(nnme='E coli', genbank_id=562), # Fine so empty
      ()),
     (dict(nnme='Streptophytina', genbank_id=131221), # Non-backbone rank
      (('WARNING', 'Streptophytina not of backbone rank, instead resolved to phylum level'),
      )),
     (dict(nnme='C marina', genbank_id=1000), # Non-backbone rank
      (('WARNING', 'NCBI ID 1000 has been superseeded by ID 107401'),
      ))
     ])
def test_validate_id_lookup(caplog, test_input, expected_log_entries, fixture_ncbi_validators):
    """This test checks that the function to generate NCBITaxon objects by searching
    an ID logs the correct errors and warnings. At the moment this for the Remote
    validator, but if a local validator is defined this should also be checked against.

    """

    fnd_tx = fixture_ncbi_validators.id_lookup(**test_input)

    if len(expected_log_entries) != len(caplog.records):
        pytest.fail('Incorrect number of log records emitted')

    # Note that this implicitly asserts that the order of logging messages is correct too
    for idx, log_rec in enumerate(caplog.records):
        exp_rec = expected_log_entries[idx]
        assert log_rec.levelname == exp_rec[0]
        assert exp_rec[1] in log_rec.message

# Third function that checks that id_lookup throws the appropriate errors
@pytest.mark.parametrize(
    'test_input,expected_exception',
    [(None,  # no parameters
      TypeError),
     (dict(nnme='E coli', genbank_id='invalid_string'),  # a string
      ValueError),
     (dict(nnme='E coli', genbank_id=27.5),  # non-integer ID
      ValueError),
     (dict(nnme=27, genbank_id=27),  # named using a number
      TypeError),
     (dict(nnme='E coli', genbank_id=-1),  # bad ID
      ValueError),
     (dict(nnme='E coli', genbank_id=100000000000000),  # bad ID
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
      ("E coli", 562, None, False, "species", "Escherichia coli")),
     (dict(nnme="Entero",taxah={'order': 'Enterobacterales', 'family': 'Enterobacteriaceae'}),
      ("Entero", 543, None, False, "family", "Enterobacteriaceae")),
     (dict(nnme="E coli strain",taxah={'species': 'Escherichia coli', 'strain':
     'Escherichia coli 1-110-08_S1_C1'}),
      ("E coli strain", 1444049, "strain", False, "species", "Escherichia coli")),
     (dict(nnme="Strepto",taxah={'phylum': 'Streptophyta', 'subphylum': 'Streptophytina'}),
      ("Strepto", 131221, "subphylum", False, "phylum", "Streptophyta")),
     (dict(nnme="Opistho",taxah={'superkingdom': 'Eukaryota', 'clade': 'Opisthokonta'}),
      ("Opistho", 33154, "clade", False, "superkingdom", "Eukaryota")),
     (dict(nnme="Vulpes",taxah={'genus': 'Vulpes', 'species': 'Vulpes vulpes'}),
      ("Vulpes", 9627, None, False, "species", "Vulpes vulpes")),
     (dict(nnme="M morus",taxah={'family': 'Moraceae', 'genus': 'Morus'}),
      ("M morus", 3497, None, False, "genus", "Morus")),
     (dict(nnme="S morus",taxah={'family': 'Sulidae', 'genus': 'Morus'}),
      ("S morus", 37577, None, False, "genus", "Morus")),
     (dict(nnme="C morus",taxah={'phylum': 'Chordata', 'genus': 'Morus'}),
      ("C morus", 37577, None, False, "genus", "Morus")),
     (dict(nnme="T maritimum",taxah={'genus': 'Tenacibaculum', 'species':
     'Tenacibaculum maritimum'}),
      ("T maritimum", 107401, None, False, "species", "Tenacibaculum maritimum")),
     (dict(nnme="C marina",taxah={'genus': 'Cytophaga', 'species': 'Cytophaga marina'}),
      ("C marina", 107401, None, True, "species", "Tenacibaculum maritimum")),
     (dict(nnme="Bacteria",taxah={'superkingdom': 'Bacteria'}),
      ("Bacteria", 2, None, False, "superkingdom", "Bacteria"))
    ])
def test_taxa_search(fixture_ncbi_validators, test_input, expected):
    """This test checks the results of searching for a taxa in the NCBI taxonomy
    database against what was expected to be found. At the moment this for the
    Remote validator, but if a local validator is defined this should also be
    checked against.

    """
    fnd_tx = fixture_ncbi_validators.taxa_search(**test_input)

    assert fnd_tx.name == expected[0]
    assert fnd_tx.genbank_id == expected[1]
    assert fnd_tx.diverg == expected[2]
    assert fnd_tx.superseed == expected[3]

    # Find last dictonary key
    f_key = list(fnd_tx.taxa_hier.keys())[-1]

    assert f_key == expected[4]
    assert fnd_tx.taxa_hier[f_key] == expected[5]

# Now test that the search function logs errors correctly
@pytest.mark.parametrize(
    'test_input,expected_log_entries',
    # Fine so empty
    [(dict(nnme='E coli', taxah={'genus': 'Escherichia', 'species': 'Escherichia coli'}),
      ()),
     (dict(nnme='Nonsense', taxah={'species': 'Nonsense garbage'}), # Nonsense
      (('ERROR', 'Taxa Nonsense cannot be found'),
      )),
     (dict(nnme='Morus', taxah={'genus': 'Morus'}), # ambigious Morus
      (('ERROR', 'Taxa Morus cannot be found using only one taxonomic level,'
      ' more should be provided'),
      )),
     # Nonsense Morus
     (dict(nnme='N Morus', taxah={'family': 'Nonsense', 'genus': 'Morus'}),
      (('ERROR', 'Provided parent taxa for N Morus not found'),
      )),
     # Aegle Morus
     (dict(nnme='A Morus', taxah={'family': 'Aegle', 'genus': 'Morus'}),
      (('ERROR', 'More than one possible parent taxa for A Morus found'),
      )),
     # Eukaryote Morus
     (dict(nnme='E Morus', taxah={'superkingdom': 'Eukaryota', 'genus': 'Morus'}),
      (('ERROR', 'Parent taxa for E Morus refers to multiple possible child taxa'),
      )),
     # Carnivora Morus
     (dict(nnme='C Morus', taxah={'superkingdom': 'Carnivora', 'genus': 'Morus'}),
      (('ERROR', 'Parent taxa not actually a valid parent of C Morus'),
      )),
     # Cytophaga marina
     (dict(nnme='C marina', taxah={'genus': 'Cytophaga', 'species': 'Cytophaga marina'}),
      (('WARNING', 'Cytophaga marina not accepted usage should be Tenacibaculum'
      ' maritimum instead'),
      )),
     # Cytophaga marina
     (dict(nnme='E coli strain', taxah={'strain': 'Escherichia coli 1-110-08_S1_C1'}),
      (('WARNING', "No backbone ranks provided in E coli strain's taxa hierarchy"),
       ('WARNING', "E coli strain not of backbone rank, instead resolved to species level"),
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

    if len(expected_log_entries) != len(caplog.records):
        pytest.fail('Incorrect number of log records emitted')

    # Note that this implicitly asserts that the order of logging messages is correct too
    for idx, log_rec in enumerate(caplog.records):
        exp_rec = expected_log_entries[idx]
        assert log_rec.levelname == exp_rec[0]
        assert exp_rec[1] in log_rec.message

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

# Then do the same for the validate_and_add_taxon function

# SECTION TESTING EXPECTED OUTPUT

# Now test that the search function logs errors correctly
@pytest.mark.parametrize(
    'test_input,expected_log_entries',
    # Fine so empty
    [(['E coli', {'species': 'Escherichia coli'}, None],
      (('INFO', "Taxon (E coli) found in NCBI database"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Escherichia coli) accepted in GBIF"),
      )),
     # Same but with valid code provided
     (['E coli', {'species': 'Escherichia coli'}, 562],
      (('INFO', "Taxon (E coli) found in NCBI database"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Escherichia coli) accepted in GBIF"),
      )),
     # whitespace padding error
     ([' E coli', {'species': 'Escherichia coli'}, None],
      (('ERROR', "Worksheet name has whitespace padding: ' E coli'"),
       ('INFO', "Taxon (E coli) found in NCBI database"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Escherichia coli) accepted in GBIF"),
      )),
     # String of just whitespace provided as name
     ([' ', {'species': 'Escherichia coli'}, None],
      (('ERROR', "Worksheet name missing, whitespace only or not text"),
      )),
     # No name error
     ([None, {'species': 'Escherichia coli'}, None],
      (('ERROR', "Worksheet name missing, whitespace only or not text"),
      )),
     # Blank string provided as name
     (['', {'species': 'Escherichia coli'}, None],
      (('ERROR', "Worksheet name missing, whitespace only or not text"),
      )),
     # Floats that can be cocnverted to integers are allowed
     (['E coli', {'species': 'Escherichia coli'}, 562.0],
      (('INFO', "Taxon (E coli) found in NCBI database"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Escherichia coli) accepted in GBIF"),
      )),
     # A true float results in multiple errors
     (['E coli', {'species': 'Escherichia coli'}, 562.5],
      (('ERROR', "NCBI ID contains value that is not an integer"),
       ('ERROR', "Improper NCBI ID provided, cannot be validated"),
      )),
     # As does a string
     (['E coli', {'species': 'Escherichia coli'}, "ID"],
      (('ERROR', "NCBI ID contains value that is not an integer"),
       ('ERROR', "Improper NCBI ID provided, cannot be validated"),
      )),
     # This checks that multiple errors can fire at once
     (['E coli', {}, 562.5],
      (('ERROR', "NCBI ID contains value that is not an integer"),
       ('ERROR', "Taxa hierarchy should be a (not empty) dictonary"),
       ('ERROR', "Improper NCBI ID provided, cannot be validated"),
       ('ERROR', "Taxon details not properly formatted, cannot validate"),
      )),
     # Taxa hierarchy provided as a string
     (['E coli', 'Escherichia coli', None],
      (('ERROR', "Taxa hierarchy should be a (not empty) dictonary"),
       ('ERROR', "Taxon details not properly formatted, cannot validate"),
      )),
     # Taxa hierarchy dictonary empty
     (['E coli', {}, None],
      (('ERROR', "Taxa hierarchy should be a (not empty) dictonary"),
       ('ERROR', "Taxon details not properly formatted, cannot validate"),
      )),
     # Example of multiple errors
     ([' E coli', {}, None],
      (('ERROR', "Worksheet name has whitespace padding: ' E coli'"),
       ('ERROR', "Taxa hierarchy should be a (not empty) dictonary"),
       ('ERROR', "Taxon details not properly formatted, cannot validate"),
      )),
     # Missing dictonary key
     (['E coli', {' ': 'Escherichia coli'}, None],
      (('ERROR', "Empty dictonary key used"),
       ('ERROR', "Taxon details not properly formatted, cannot validate"),
      )),
     # Non-string dictonary key
     (['E coli', {26: 'Escherichia coli'}, None],
      (('ERROR', "Non-string dictonary key used: 26"),
       ('ERROR', "Taxon details not properly formatted, cannot validate"),
      )),
     # Padded dictonary key
     (['E coli', {' species': 'Escherichia coli'}, None],
      (('ERROR', "Dictonary key has whitespace padding: ' species'"),
       ('INFO', "Taxon (E coli) found in NCBI database"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Escherichia coli) accepted in GBIF"),
      )),
     # Missing dictonary value
     (['E coli', {'species': ''}, None],
      (('ERROR', "Empty dictonary value used"),
       ('ERROR', "Taxon details not properly formatted, cannot validate"),
      )),
     # Non-string dictonary value
     (['E coli', {'species': 26}, None],
      (('ERROR', "Non-string dictonary value used: 26"),
       ('ERROR', "Taxon details not properly formatted, cannot validate"),
      )),
     # Padding on dictonary value
     (['E coli', {'species': ' Escherichia coli'}, None],
      (('ERROR', "Dictonary value has whitespace padding: ' Escherichia coli'"),
       ('INFO', "Taxon (E coli) found in NCBI database"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Escherichia coli) accepted in GBIF"),
      )),
     # Taxon hierarchy in wrong order
     (['E coli', {'species': 'Escherichia coli', 'genus': 'Escherichia'}, None],
      (('ERROR', "Taxon hierarchy not in correct order"),
      )),
     # Right order but non-backbone rank
     (['Strepto', {'phylum': 'Streptophyta', 'subphylum': 'Streptophytina'}, None],
      (('WARNING', "Strepto not of backbone rank, instead resolved to phylum level"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('ERROR', "Taxon (Streptophyta) No match found"),
      )),
     # Superseeded taxon name used
     (['C marina', {'species': 'Cytophaga marina'}, None],
      (('WARNING', "Cytophaga marina not accepted usage should be Tenacibaculum maritimum instead"),
       ('WARNING', "Taxonomic classification superseeded for C marina, using new taxonomic classification"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Tenacibaculum maritimum) accepted in GBIF"),
      )),
     # Nonsense taxon provided
     (['N garbage', {'species': 'Nonsense garbage'}, None],
      (('ERROR', "Taxa N garbage cannot be found"),
       ('ERROR', "Search based on taxon hierarchy failed"),
      )),
     # Ambigious taxon provided
     (['Morus', {'genus': 'Morus'}, None],
      (('ERROR', "Taxa Morus cannot be found using only one taxonomic level, more should be provided"),
       ('ERROR', "Search based on taxon hierarchy failed"),
      )),
     # Ambigious taxon resolved
     (['M Morus', {'family': 'Moraceae', 'genus': 'Morus'}, None],
      (('INFO', "Taxon (M Morus) found in NCBI database"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('ERROR', "Taxon (Morus) Multiple equal matches for Morus"),
      )),
     # E coli with incorrect NCBI ID
     (['E coli', {'species': 'Escherichia coli'}, 333],
      (('ERROR', "The NCBI ID supplied for E coli does not match hierarchy: expected 562 got 333"),
      )),
     # Non-backbone case with code
     (['E coli strain', {'species': 'Escherichia coli', 'strain': 'Escherichia coli 1-110-08_S1_C1'},
       1444049],
      (('WARNING', "E coli strain not of backbone rank, instead resolved to species level"),
       ('WARNING', "E coli strain not of backbone rank, instead resolved to species level"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Escherichia coli) accepted in GBIF"),
      )),
     # Superseeded taxa ID and name
     (['C marina', {'species': 'Cytophaga marina'}, 1000],
      (('WARNING', "Cytophaga marina not accepted usage should be Tenacibaculum maritimum instead"),
       ('WARNING', "NCBI ID 1000 has been superseeded by ID 107401"),
       ('WARNING', "Taxonomic classification superseeded for C marina, using new taxonomic classification"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Tenacibaculum maritimum) accepted in GBIF"),
      )),
     # Just superseeded taxa ID
     (['T maritimum', {'species': 'Tenacibaculum maritimum'}, 1000],
      (('WARNING', "NCBI ID 1000 has been superseeded by ID 107401"),
       ('WARNING', "NCBI taxa ID superseeded for T maritimum, using new taxa ID"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Tenacibaculum maritimum) accepted in GBIF"),
      )),
     # Just superseeded name
     (['C marina', {'species': 'Cytophaga marina'}, 107401],
      (('WARNING', "Cytophaga marina not accepted usage should be Tenacibaculum maritimum instead"),
       ('WARNING', "Taxonomic classification superseeded for C marina, using new taxonomic classification"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('INFO', "Taxon (Tenacibaculum maritimum) accepted in GBIF"),
      )),
     # E coli recorded as a family rather than a species
     (['E coli', {'family': 'Escherichia coli'}, None],
      (('ERROR', "Escherichia coli is a species not a family"),
       ('ERROR', "Search based on taxon hierarchy failed"),
      )),
     # E coli recorded as a subspecies rather than a species
     (['E coli', {'subspecies': 'Escherichia coli'}, None],
      (('ERROR', "Escherichia coli is a species not a subspecies"),
       ('ERROR', "Search based on taxon hierarchy failed"),
      )),
     # Same idea for a non-backbone case
     (['Streptophytina', {'phylum': 'Streptophytina'}, None],
      (('WARNING', "Streptophytina not of backbone rank, instead resolved to phylum level"),
       ('ERROR', "Streptophytina is a subphylum not a phylum"),
       ('ERROR', "Search based on taxon hierarchy failed"),
      )),
     # Superkingdom not included in GBIF case
     (['Eukaryota', {'superkingdom': 'Eukaryota'}, 2759],
      (('INFO', "Taxon (Eukaryota) found in NCBI database"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('WARNING', "Taxon (Eukaryota) is not of a GBIF backbone rank"),
       ('ERROR', "Taxon (Eukaryota) No match found")
      )),
     # Can actually deal with the bacterial case
     (['Bacteria', {'superkingdom': 'Bacteria'}, 2],
      (('INFO', "Taxon (Bacteria) found in NCBI database"),
       ('INFO', "Attempting to validate against GBIF database"),
       ('WARNING', "Taxon (Bacteria) is not of a GBIF backbone rank"),
       ('INFO', "Taxon (Bacteria) defined as kingdom rank in GBIF"),
       ('INFO', "Taxon (Bacteria) accepted in GBIF"),
      )),
     ])
def test_validate_and_add_taxon(caplog, test_input, expected_log_entries,
                                fixture_taxon_validators):
    """This test checks that the function that searches the NCBI database to find
    information on particular taxa logs the correct errors and warnings. At the
    moment this for the Remote validator, but if a local validator is defined
    this should also be checked against.

    """

    # ONLY MAKING A LOCAL VERSION FOR NOW
    gb_instance = genb_taxa.GenBankTaxa()
    fnd_tx = gb_instance.validate_and_add_taxon(fixture_taxon_validators,test_input)

    if len(expected_log_entries) != len(caplog.records):
        pytest.fail('Incorrect number of log records emitted')

    # Note that this implicitly asserts that the order of logging messages is correct too
    for idx, log_rec in enumerate(caplog.records):
        exp_rec = expected_log_entries[idx]
        assert log_rec.levelname == exp_rec[0]
        assert exp_rec[1] in log_rec.message

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
def test_validate_and_add_taxon_errors(fixture_taxon_validators, test_input, expected_exception):
    """This test checks validator.validate_and_add_taxon inputs throw errors as expected
    """

    # ONLY MAKING A LOCAL VERSION FOR NOW
    gb_instance = genb_taxa.GenBankTaxa()

    with pytest.raises(expected_exception):
        _ = gb_instance.validate_and_add_taxon(fixture_taxon_validators,test_input)

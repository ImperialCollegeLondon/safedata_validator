import pytest
from safedata_validator import taxa


@pytest.fixture(scope='module', params=['remote', 'local'])
def validators(resources_with_local_gbif, request):
    """Parameterised fixture to return local and remote taxon validators
    """
    if request.param == 'remote':
        return taxa.RemoteGBIFValidator()

    elif request.param == 'local':
        return taxa.LocalGBIFValidator(resources_with_local_gbif)

# ------------------------------------------
# Testing Taxon
# ------------------------------------------


@pytest.mark.parametrize(
    'test_input,expected_exception',
    [(dict(),  # no parameters
      TypeError),
     (dict(name=1, rank=1),  # non string name and rank
      TypeError),
     (dict(name='Bombus bombus', rank='species', gbif_id='error'),  # non-numeric gbif_id
      TypeError),
     (dict(name='Bombus bombus', rank='species', gbif_id=3.141),  # non-integer gbif_id
      ValueError)])
def test_taxon_init_errors(test_input, expected_exception):
    """This test checks taxon inputs throw errors as expected
    """
    with pytest.raises(expected_exception):
        _ = taxa.Taxon(**test_input)

# ------------------------------------------
# Testing taxon validators
# ------------------------------------------


@pytest.mark.parametrize(
    'test_input,expected',
    [(dict(name='Crematogaster borneensis', rank='species'),
      ('found', True, 1324716, None)),
     (dict(name='Morus', rank='genus', gbif_id=2480962),
      ('found', True, 2480962, None)),
     (dict(name='Alsomitra simplex', rank='species'),
      ('found', False, 5537041, 3623287))
     ])
def test_validator_search(validators, test_input, expected):
    """This test checks inputs against expected outputs for both the local
    and remote validator classes.

    TODO - proof of concept, think about systematic structure and failure
           cases
    """

    tx = taxa.Taxon(**test_input)
    srch_out = validators.search(tx)

    assert srch_out.lookup_status == expected[0]
    assert srch_out.is_canon == expected[1]
    assert srch_out.gbif_id == expected[2]

    if not srch_out.is_canon:
        assert srch_out.canon_usage.gbif_id == expected[3]


@pytest.mark.parametrize(
    'test_input,expected',
    [(1324716, ('found', 'Crematogaster borneensis')),
     (2480962, ('found', 'Morus'))])
def test_validator_gbif_lookup_outputs(validators, test_input, expected):
    """This test checks inputs against expected outputs for both the local
    and remote validator classes.
    TODO - proof of concept, think about systematic structure and failure
    cases
    """

    found = validators.id_lookup(test_input)

    assert found.lookup_status == expected[0]
    assert found.name == expected[1]


@pytest.mark.parametrize(
    'test_input,expected_exception',
    [(None,  # no parameters
      TypeError),
     ('invalid_string',  # a string
      TypeError),
     (-1,  # bad ID
      ValueError),
     (100000000000000,  # bad ID
      taxa.GBIFError)])
def test_validator_gbif_lookup_errors(validators, test_input, expected_exception):
    """This test checks validator.id_lookup inputs throw errors as expected
    """

    with pytest.raises(expected_exception):
        _ = validators.id_lookup(test_input)

# ------------------------------------------
# Testing Taxa _validate_tuple
# - This private method tests the individual tuples of taxon and parent
#   provided in the Taxa sheet, separating loading from testing. The test
#   cases test the various combinations of good and bad inputs.
# ------------------------------------------


@pytest.mark.parametrize(
    'taxon_tuple,expected_log',
    [(('Crematogaster borneensis',
            ('Crematogaster borneensis', 'Species', None, None),
            None),
       'Taxon found in GBIF backbone (accepted)'),
      (('Dolichoderus sp.',
            ('Dolichoderus', 'Genus', None, None),
            None),
       'Taxon found in GBIF backbone (accepted)'),
       (('Ponerinae #1',
            ('Ponerinae', 'Subfamily', None, None),
            ('Formicidae', 'Family', None, None)),
       'subfamily with valid parent information'),
    ])
def test_validate_tuple(caplog, resources_local_and_remote, taxon_tuple, expected_log):

    taxa_instance = taxa.Taxa(resources_local_and_remote)
    taxa_instance.validate_tuple(taxon_tuple)

    assert expected_log in caplog.text


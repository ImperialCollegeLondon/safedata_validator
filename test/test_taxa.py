import pytest
from .safedata_validator import taxa
from .safedata_validator import config


# ------------------------------------------
# Fixtures: a parameterized fixture containing examples of both local and
#     remote GBIF validators.
#
# TODO - compare local and remote outputs:
#  https://stackoverflow.com/questions/56558823/pytest-run-all-tests-twice-and-compare-results-bet-mock-and-real
# ------------------------------------------


@pytest.fixture(scope='module', params=['remote', 'local'])
def validator(request):
    """Parameterised fixture to run tests using both the local and remote validator.
    """
    if request.param == 'remote':
        return taxa.RemoteGBIFValidator()
    else:
        cfg = config.config()
        return taxa.LocalGBIFValidator(cfg.gbif_database)

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
# Testing Validators
# ------------------------------------------


@pytest.mark.parametrize(
    'test_input,expected',
    [(dict(name='Crematogaster borneensis', rank='species'),
      ('found', True, 1324716)),
     (dict(name='Morus', rank='genus', gbif_id=2480962),
      ('found', True, 2480962))])
def test_validator_search(validator, test_input, expected):
    """This test checks inputs against expected outputs for both the local
    and remote validator classes. 
    TODO - proof of concept, think about systematic structure and failure
    cases
    """

    tx = taxa.Taxon(**test_input)
    found = validator.search(tx)
    
    assert found.lookup_status == expected[0]
    assert found.is_canon == expected[1]
    assert found.gbif_id == expected[2]


@pytest.mark.parametrize(
    'test_input,expected',
    [(1324716, ('found', 'Crematogaster borneensis')),
     (2480962, ('found', 'Morus'))])
def test_validator_gbif_lookup_outputs(validator, test_input, expected):
    """This test checks inputs against expected outputs for both the local
    and remote validator classes.
    TODO - proof of concept, think about systematic structure and failure
    cases
    """

    found = validator.id_lookup(test_input)

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
def test_validator_gbif_lookup_errors(validator, test_input, expected_exception):
    """This test checks validator.id_lookup inputs throw errors as expected
    """

    with pytest.raises(expected_exception):
        _ = validator.id_lookup(test_input)

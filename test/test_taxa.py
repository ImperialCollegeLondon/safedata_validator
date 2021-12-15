from logging import log
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
#   cases test the various combinations of good and bad inputs and checks
#   that the set of emitted logging records is as expected.
# ------------------------------------------


@pytest.mark.parametrize(
    'taxon_tuple,expected_log_entries',
      [ # PARENTS: Testing checking of parent taxon code using morphospecies
        # Non backbone parent
        (('Morphospecies', ('mspecies', 'Morphospecies', None, None), 
                           ('Agroecomyrmecini', 'Tribe', None, None)),
         (('ERROR', 'is not of a backbone rank'),
          ('ERROR', 'invalid parent information'),
          )), 
        # Unknown parent
        (('Morphospecies', ('mspecies', 'Morphospecies', None, None), 
                           ('Madeupnotarealfamilyidae', 'Family', None, None)),
         (('ERROR', 'No match found'),
          ('ERROR', 'invalid parent information'),
          )), 
        # Synonymous parent
        (('Morphospecies', ('mspecies', 'Morphospecies', None, None), 
                           ("Bothroponera", "Genus", None, None)),
         (('WARNING', 'considered a synonym'),
          ('INFO', 'has valid parent information'),
          )), 
        # Normal valid parent
        (('Morphospecies', ('mspecies', 'Morphospecies', None, None),
                           ('Formicidae', 'Family', None, None)),
         (('INFO', 'accepted'),
          ('INFO', 'has valid parent information'),
         )),
        # IGNORED: Testing ignored backbone matches
        # Ignored non backbone (no parent)
        (('Morphospecies', ('mspecies', 'Morphospecies', None, 123456789), 
                           None),
         (('INFO', 'No parent taxon provided'),
          ('ERROR', 'Ignore ID can only be used with GBIF backbone taxon ranks'),
          ('ERROR', 'Taxa with Ignore ID must provide parent information'),
         )),
        # Ignored unknown (no parent)
        (('Unknown species', ('Unknown species', 'species', None, 123456789), 
                           None),
         (('INFO', 'No parent taxon provided'),
          ('ERROR', 'Taxon with Ignore ID not found in GBIF backbone'),
          ('ERROR', 'Taxa with Ignore ID must provide parent information'),
         )),
        # - Ignoring good GBIF match with no parent
        (('Crematogaster borneensis', ('Crematogaster borneensis', 'Species', None, 1324716), 
                                      None),
         (('INFO', 'No parent taxon provided'),
          ('INFO', 'Canon GBIF usage ignored'),
          ('ERROR', 'Taxa with Ignore ID must provide parent information'),
         )),
         # - Ignoring GBIF match with bad parent
        (('Crematogaster borneensis', ('Crematogaster borneensis', 'Species', None, 1324716), 
                                      ('Madeupnotarealfamilyidae', 'Family', None, None)),
         (('ERROR', 'No match found'),
          ('INFO', 'Canon GBIF usage ignored'),
          ('ERROR', 'Taxon with Ignore ID has invalid parent information'),
         )),
        # - Ignoring GBIF match on canon taxa with good parent and correct ID
        (('Panthera leo', ('Panthera leo', 'Species', None, 5219404),
                              ('Felis', 'Genus', None, None)),
         (('INFO', 'accepted'),
          ('INFO', 'Canon GBIF usage ignored'),
          ('INFO', 'Taxon with ignored canon usage has valid parent information'),
         )),
        # - Ignoring GBIF match on canon taxa with good parent and bad ID
        (('Panthera leo', ('Panthera leo', 'Species', None, 123456789),
                              ('Felis', 'Genus', None, None)),
         (('INFO', 'accepted'),
          ('ERROR', 'Ignore ID does not match the canon GBIF usage'),
          ('INFO', 'Taxon with ignored canon usage has valid parent information'),
         )),
        # - Ignoring GBIF match on synonym with good parent but non canon ID (synonym)
        (('Felis imperialis', ('Felis imperialis', 'Species', None, 8531298),
                              ('Felis', 'Genus', None, None)),
         (('INFO', 'accepted'),
          ('ERROR', 'Taxon is non-canon and Ignore ID does not match the canon GBIF usage'),
          ('INFO', 'Taxon with ignored canon usage has valid parent information'),
         )),
        # - Ignoring GBIF match on synonym with good parent using correct canon ID
        (('Felis imperialis', ('Felis imperialis', 'Species', None, 5219404),
                              ('Felis', 'Genus', None, None)),
         (('INFO', 'accepted'),
          ('INFO', 'Canon GBIF usage ignored'),
          ('INFO', 'Taxon with ignored canon usage has valid parent information'),
         )),
        # NON-BACKBONE Testing non backbone outcomes
        # No parent on non-backbone
        (('Morphospecies', ('mspecies', 'Morphospecies', None, None),
                           None),
         (('INFO', 'No parent taxon provided'),
          ('ERROR', ' must provide parent information'),
         )),
        # Invalid and valid parent cases tested above already but check main taxon logging
        # Unknown parent
        (('Morphospecies', ('mspecies', 'Morphospecies', None, None), 
                           ('Madeupnotarealfamilyidae', 'Family', None, None)),
         (('ERROR', 'No match found'),
          ('ERROR', 'has invalid parent information'),
         )),
        # Normal valid parent
        (('Morphospecies', ('mspecies', 'Morphospecies', None, None),
                           ('Formicidae', 'Family', None, None)),
         (('INFO', 'accepted'),
          ('INFO', 'has valid parent information'),
         )),
        # BACKBONE: Testing backbone outcomes
        # - Good match - no parent
        (('Crematogaster borneensis', ('Crematogaster borneensis', 'Species', None, None), 
                                      None),
         (('INFO', 'No parent taxon provided'),
          ('INFO', 'Taxon found in GBIF backbone (accepted)'),
         )),
        # - Synonymous match - no parent
        (('Felis imperialis', ('Felis imperialis', 'Species', None, None),
                              None),
         (('INFO', 'No parent taxon provided'),
         ('WARNING', 'Taxon considered a synonym'),
         )),
        # - Good match, valid parent
        (('Crematogaster borneensis', ('Crematogaster borneensis', 'Species', None, None), 
                                      ('Crematogaster', 'Genus', None, None)),
         (('INFO', 'accepted'),
          ('INFO', 'Taxon in GBIF backbone (accepted) with compatible parent information'),
         )),
        # - Good match, incompatible parent
        (('Crematogaster borneensis', ('Crematogaster borneensis', 'Species', None, None), 
                                      ('Pinus', 'Genus', None, None)),
         (('INFO', 'accepted'),
          ('ERROR', 'Taxon in GBIF backbone (accepted) with incompatible parent information'),
         )),
        # - Good match, incompatible parent
        (('Crematogaster borneensis', ('Crematogaster borneensis', 'Species', None, None), 
                                      ('Madeupnotarealgenus', 'Genus', None, None)),
         (('ERROR', 'No match found'),
          ('ERROR', 'Taxon in GBIF backbone (accepted) but with invalid parent information'),
         )),
        #  - No match, no parent
        (("Crematogaster ormei", ("Crematogaster ormei", "Species", None, None), 
                                 None),
        (('INFO', 'No parent taxon provided'),
         ('ERROR', 'Taxon name and rank combination not found'),
         )),
        #  - No match, bad parent
        (("Crematogaster ormei", ("Crematogaster ormei", "Species", None, None), 
                                 ('Madeupnotarealgenus', 'Genus', None, None)),
        (('ERROR', 'No match found'),
         ('ERROR', 'Taxon not found in GBIF and has invalid parent information'),
         )),
        # - No match, valid parent
        (("Crematogaster ormei", ("Crematogaster ormei", "Species", None, None), 
                                 ('Crematogaster', 'Genus', None, None)),
        (('INFO', 'accepted'),
         ('INFO', 'Taxon not found in GBIF but has valid parent information'),
         )),
        # ODDITIES
        # - Subspecies
        (("Water monitor", ("Varanus salvator macromaculatus", "Subspecies", None, None), 
                            None),
         (('INFO', 'No parent taxon provided'),
          ('INFO', 'Taxon found in GBIF backbone (accepted)'),
         )),
        # - Ambiguous
        (("Gannets", ("Morus", "Genus", None, None), 
                      None),
         (('INFO', 'No parent taxon provided'),
          ('ERROR', 'GBIF issue:'),
         )),
        # - Ambiguous resolved with GBIF ID 
        (("Gannets", ("Morus", "Genus", 2480962, None), 
                      None),
         (('INFO', 'No parent taxon provided'),
          ('INFO', 'Taxon found in GBIF backbone (accepted)'),
         )),   
        #  - Ambiguous parent 
        (("Solenopsis #1", ("NA", "Morphospecies", None, None), 
                           ("Solenopsis", "Genus", None, None)),
         (('ERROR', 'Multiple equal matches'),
          ('ERROR', 'invalid parent information'),
         )),
        #  - Ambiguous parent resolved
        (("Solenopsis #1", ("NA", "Morphospecies", None, None), 
                           ("Solenopsis", "Genus", 9107211, None)),
         (('INFO', 'accepted'),
          ('INFO', 'valid parent information'),
         )),
])
def test_validate_tuple(caplog, resources_local_and_remote, taxon_tuple, expected_log_entries):

    taxa_instance = taxa.Taxa(resources_local_and_remote)
    taxa_instance.validate_tuple(taxon_tuple)

    if len(expected_log_entries) != len(caplog.records):
        pytest.fail('Incorrect number of log records emitted')

    # Note that this implicitly asserts that the order of logging messages is correct too
    for idx, log_rec in enumerate(caplog.records):
        exp_rec = expected_log_entries[idx]
        assert log_rec.levelname == exp_rec[0]
        assert exp_rec[1] in log_rec.message

# ("Alsomitra simplex", ("Alsomitra simplex", "Species", None, None), None)
# ("Gannets", ("Morus", "Genus", "2480962", None), None),
# ("Zenicomus photuroides", ("Zenicomus photuroides", "Species", None, None), None),
# ("Cicada sanguinolenta", ("Cicada sanguinolenta", "Species", None, , None), ("Cicada", "Genus", None, None)
# ("Melittia oedippus", ("Melittia oedippus", "Species", None, None, None, None
# ("Melaphorus potteri", ("Melaphorus potteri", "Species", None, None, None, None
# ("Goniopholis tenuidens", ("Goniopholis tenuidens", "Species", None, None, None, None
# ("Solenopsis abdita", ("Solenopsis abdita", "Species", None, None, None, None
# ("Solenopsis #1", ("NA", "Morphospecies", None, "Solenopsis", "Genus", "9107211"
# ("Biarmosuchus tagax", ("Biarmosuchus tagax", "Species", None, None, None, None
# ("Camponotites kraussei", ("Camponotites kraussei", "Species", None, None, None, None
# ("Bothroponera novus", ("Bothroponera novus", "Species", None, ), ("Bothroponera", "Genus", None, None)

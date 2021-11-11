
# flow of process

load workbook 

        if config.gbif_database is None:
            LOGGER.info('Using GBIF online API to validate taxonomy')
            self.validator = RemoteGBIFValidator
        else:
            LOGGER.info('Validating local GBIF database: ' + config.gbif_database)
            self.validator = LocalGBIFValidator(config.gbif_database)

# reorder sys.path to get the local copy not the pip installed
sys.path = sys.path[:6] + sys.path[7:] + [sys.path[6]]

import xlrd
from safedata_validator.config import config
from safedata_validator.taxa import Taxa, LocalGBIFValidator, Taxon

cfg = config()
taxa = Taxa()
val = LocalGBIFValidator(cfg.gbif_database)

workbook = xlrd.open_workbook('/Users/dorme/Research/SAFE/Database/safedata_validator_package/test_files/Test_format_good.xlsx')
taxa.load_taxa(workbook)
taxa.validate(val)

tax = val.id_lookup(5537041)

tx = Taxon('Alsomitra simplex', 'species', 1)
tax2 = val.search(tx)
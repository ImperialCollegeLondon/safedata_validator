import sqlite3
import csv
# This script uses a list of taxa used in testing to cut down the full (~2GB)
# GBIF backbone to a truncated version for inclusion in the test fixtures
# directory for testing local taxon validation. The list cannot be easily
# automated, because e.g. ambiguous taxa detection needs to have copies of
# all possible resolutions. Currently easiest to hand curate this list but
# the details for good matches can be extracted from test files and then
# supplemented. For example

# from safedata_validator.taxa import Taxa
# from safedata_validator.resources import Resources
#
# rs = Resources()
# tx = Taxa(rs)
# wb = openpyxl.load_workbook('Test_format_good.xlsx')
# tx.load(wb['Taxa'])
# [t[1:] for t in tx.taxon_index]  # these values are the row entries in the list

# The details file contains some whole row comments, so these are skipped
fp = open('details_of_test_taxa.csv')
rdr = csv.DictReader(filter(lambda row: row[0] != '#', fp))
data = list(rdr)
fp.close()

# Get unique, non-user GBIF IDs from the data
gbif_ids = set([int(d['gbif_id']) for d in data])
gbif_ids.discard(-1)

# Copy truncated database from local non-git copy.
source_db = sqlite3.connect('../../local_db/gbif_2021-11-26/gbif_backbone.sqlite3')
dest_db = sqlite3.connect('gbif_backbone_truncated.sqlite')

# create table
cur = source_db.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='backbone'")
schema = cur.fetchone()[0]
cur = dest_db.execute(schema)

# copy rows
for each_id in gbif_ids:

    cur = source_db.execute(f"SELECT * FROM backbone WHERE id = {each_id}")
    ins = dest_db.execute('insert into backbone values (' + ','.join(['?'] * 29) + ')' , cur.fetchone())

dest_db.commit()

source_db.close()
dest_db.close()

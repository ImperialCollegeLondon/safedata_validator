"""Script to generate a truncated version of the GBIF taxonomic backbone.

This script uses a list of taxa used in testing to cut down the full (~2GB) GBIF
backbone to a truncated version for inclusion in the test fixtures directory for testing
local taxon validation. The list cannot be easily automated, because e.g. ambiguous taxa
detection needs to have copies of all possible resolutions. Currently easiest to hand
curate this list but the details for good matches can be extracted from test files and
then supplemented. For example:

from safedata_validator.taxa import Taxa
from safedata_validator.resources import Resources

rs = Resources()
tx = Taxa(rs)
wb = openpyxl.load_workbook('Test_format_good.xlsx')
tx.load(wb['Taxa'])
[t[1:] for t in tx.taxon_index]  # these values are the row entries in the list
"""

import csv
import sqlite3

from safedata_validator.resources import Resources

# Get the locally configured full databases
resources = Resources()

# The details file contains some whole row comments, so these are skipped
fp = open("test_gbif_taxa_details.csv")
rdr = csv.DictReader(filter(lambda row: row[0] != "#", fp))
data = list(rdr)
fp.close()

# Get unique, non-user GBIF IDs from the data
gbif_ids = set([int(d["gbif_id"]) for d in data])
gbif_ids.discard(-1)

# Copy truncated database from local non-git copy.
source_db = sqlite3.connect(resources.gbif_database)
dest_db = sqlite3.connect("gbif_backbone_truncated.sqlite")

# create backbone table
cur = source_db.execute(
    "SELECT sql FROM sqlite_master WHERE type='table' AND name='backbone'"
)
schema = cur.fetchone()[0]
cur = dest_db.execute(schema)

# copy rows
for each_id in gbif_ids:

    cur = source_db.execute(f"SELECT * FROM backbone WHERE id = {each_id}")
    ins = dest_db.execute(
        "insert into backbone values (" + ",".join(["?"] * 28) + ")", cur.fetchone()
    )

dest_db.commit()

# timestamp table
cur = source_db.execute(
    "SELECT sql FROM sqlite_master WHERE type='table' AND name='timestamp'"
)
schema = cur.fetchone()[0]
cur = dest_db.execute(schema)

cur = source_db.execute("SELECT * FROM timestamp")
timestamp = cur.fetchone()[0]
ins = dest_db.execute(f'insert into timestamp values ("{timestamp}")')

dest_db.commit()

source_db.close()
dest_db.close()

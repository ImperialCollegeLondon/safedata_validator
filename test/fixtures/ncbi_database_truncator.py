"""Script to generate a truncated version of the NCBI taxonomy backbone.

This script uses a list of taxa used in testing to cut down the full (~0.5GB) NCBI
taxonomy database to a truncated version for inclusion in the test fixtures directory
for testing local taxon validation. The input is a list of valid NCBI taxa IDs,
associated names, nodes, and merged nodes are then all retained in the appropriate
tables. Everything else is chucked out. Only IDs for taxa used are needed, higher taxa
are found automatically and added to the database.
"""

import csv
import sqlite3

from safedata_validator.resources import Resources

# Get the locally configured full databases
resources = Resources()

# The details file contains some whole row comments, so these are skipped
fp = open("test_ncbi_taxa_details.csv")
rdr = csv.DictReader(filter(lambda row: row[0] != "#", fp))
data = list(rdr)
fp.close()

# Copy truncated database from local non-git copy.
source_db = sqlite3.connect(resources.ncbi_database)
dest_db = sqlite3.connect("ncbi_database_truncated.sqlite")

# Get unique, non-user GBIF IDs from the data
ncbi_ids = set([int(d["ncbi_id"]) for d in data])

# Store user provided IDs to iterate over
prov_ids = list(ncbi_ids)

# Then setup loop to find the whole lineage for each id
for ind in range(0, len(prov_ids)):
    # Find details of provided ID
    sql = f"select * from nodes where tax_id = {prov_ids[ind]}"
    taxon_row = source_db.execute(sql).fetchone()

    # Setup while loop to find lineage
    lin_fnd = False
    while not lin_fnd:
        # Find node and name of the parent taxon
        sql = f"select * from nodes where tax_id = {taxon_row[1]}"
        taxon_row = source_db.execute(sql).fetchone()
        # Store all relevant info
        ncbi_ids.add(taxon_row[0])
        # End this when the parent taxon is root (ID=1)
        if taxon_row[1] == 1:
            lin_fnd = True

# create the three tables
tables = ["nodes", "names", "merge"]
for t_name in tables:
    cur = source_db.execute(
        f"SELECT sql FROM sqlite_master WHERE type='table' AND name='{t_name}'"
    )
    schema = cur.fetchone()[0]
    cur = dest_db.execute(schema)

# copy rows
for each_id in ncbi_ids:
    # Search for the node associated with the ID
    cur = source_db.execute(f"SELECT * FROM nodes WHERE tax_id = {each_id}")
    ins = dest_db.execute(
        "insert into nodes values (" + ",".join(["?"] * 4) + ")", cur.fetchone()
    )
    # Search for all names associated with the ID
    cur = source_db.execute(f"SELECT * FROM names WHERE tax_id = {each_id}")
    rqst = cur.fetchall()
    for ind in range(0, len(rqst)):
        ins = dest_db.execute(
            "insert into names values (" + ",".join(["?"] * 4) + ")", rqst[ind]
        )
    # Search for all IDs that have been merged into this ID
    cur = source_db.execute(f"SELECT * FROM merge WHERE new_tax_id = {each_id}")
    rqst = cur.fetchall()
    for ind in range(0, len(rqst)):
        ins = dest_db.execute(
            "insert into merge values (" + ",".join(["?"] * 2) + ")", rqst[ind]
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

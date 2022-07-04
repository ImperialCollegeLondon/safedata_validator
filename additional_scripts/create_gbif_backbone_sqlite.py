"""
This script is used to build an SQLite3 database file containing the contents of
the GBIF backbone taxonomy. It is not currently provided as a part of the package
because the details of file names etc. change and so it is provided as a recipe.
"""

import csv
import sqlite3

db_file = "local_db/gbif_2021-11-26/gbif_backbone.sqlite3"
drop_fields = ["name_published_in", "issues"]
backbone = "local_db/gbif_2021-11-26/simple.txt"
deleted = "local_db/gbif_2021-11-26/simple-deleted.txt"

# Create the schema for the table, using drop fields to remove unwanted fields
# in the schema and the data tuples. The file_schema list describes the full set
# of fields provided by GBIF
file_schema = [
    ("id", "int PRIMARY KEY"),
    ("parent_key", "int"),
    ("basionym_key", "int"),
    ("is_synonym", "boolean"),
    ("status", "text"),
    ("rank", "text"),
    ("nom_status", "text[]"),
    ("constituent_key", "text"),
    ("origin", "text"),
    ("source_taxon_key", "int"),
    ("kingdom_key", "int"),
    ("phylum_key", "int"),
    ("class_key", "int"),
    ("order_key", "int"),
    ("family_key", "int"),
    ("genus_key", "int"),
    ("species_key", "int"),
    ("name_id", "int"),
    ("scientific_name", "text"),
    ("canonical_name", "text"),
    ("genus_or_above", "text"),
    ("specific_epithet", "text"),
    ("infra_specific_epithet", "text"),
    ("notho_type", "text"),
    ("authorship", "text"),
    ("year", "text"),
    ("bracket_authorship", "text"),
    ("bracket_year", "text"),
    ("name_published_in", "text"),
    ("issues", "text[]"),
]

# Get a logical index of which fields are being kept
drop_index = [True if vl[0] in drop_fields else False for vl in file_schema]

# Create the final schema for the backbone table, including a new field to show
# deleted taxa and insert statements for both the backbone and deleted files
output_schema = ", ".join(
    [" ".join(val) for val, drop in zip(file_schema, drop_index) if not drop]
)
output_schema = f"CREATE TABLE backbone ({output_schema}, deleted boolean)"

insert_placeholders = ",".join(["?"] * (len(drop_index) - sum(drop_index)))
insert_statement = f"INSERT INTO backbone VALUES ({insert_placeholders}, False)"

insert_deleted_statement = f"INSERT INTO backbone VALUES ({insert_placeholders}, True)"

# Create the output file and turn off safety features for speed
con = sqlite3.connect(db_file)
con.execute("PRAGMA synchronous = OFF")

# Create the table
con.execute(output_schema)
con.commit()

# Import data from the simple backbone and deleted taxa

# # The approach below is more efficient but makes it impossible to drop fields
# # and substitute \\N to None. Although converting \\N to None can be done later with
# # an update statement, you _cannot_ drop fields in sqlite3, so that has
# # to be done up front.
#
# con.executemany(
#     insert_statement,
#     bb_reader
# )

with open(backbone) as bbn:

    # The files are tab delimited but the quoting is sometimes unclosed,
    # so turning off quoting - includes quotes in the fields where present
    bb_reader = csv.reader(bbn, delimiter="\t", quoting=csv.QUOTE_NONE)

    for row in bb_reader:
        row = [
            None if val == "\\N" else val
            for val, drp in zip(row, drop_index)
            if not drp
        ]

        con.execute(insert_statement, row)

with open(deleted) as dlt:

    # The files are tab delimited but the quoting is sometimes unclosed,
    # so turning off quoting - includes quotes in the fields where present
    dl_reader = csv.reader(dlt, delimiter="\t", quoting=csv.QUOTE_NONE)

    for row in dl_reader:
        row = [
            None if val == "\\N" else val
            for val, drp in zip(row, drop_index)
            if not drp
        ]

        con.execute(insert_deleted_statement, row)

# Create the indices
con.execute("CREATE INDEX backbone_name_rank ON backbone (canonical_name, rank);")
con.execute("CREATE INDEX backbone_id ON backbone (id);")
con.commit()

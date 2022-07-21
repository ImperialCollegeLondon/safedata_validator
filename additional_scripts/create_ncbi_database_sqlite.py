"""
This script is used to build an SQLite3 database file containing the contents of
the NCBI taxonomy. It is not currently provided as a part of the package because
the details of file names etc. change and so it is provided as a recipe.
"""

import csv
import sqlite3

db_file = "local_db/ncbi/ncbi_database.sqlite3"
# Set fields to be dropped for each data table
nodes_drop_fields = [
    "embl_code",
    "division_id",
    "inherited_div_flag",
    "genetic_code_id",
    "inherited_GC_flag",
    "mito_code_id",
    "inherited_MGC_flag",
    "GenBank_hidden_flag",
    "hidden_subtree_root_flag",
]
names_drop_fields = ["unique_name"]
merge_drop_fields: list[str] = []
nodes = "local_db/ncbi/nodes.dmp"
names = "local_db/ncbi/names.dmp"
merge = "local_db/ncbi/merged.dmp"

# Create the schema for the table, using drop fields to remove unwanted fields
# in the schema and the data tuples. The file_schema list describes the full set
# of fields provided by NCBI
nodes_file_schema = [
    ("tax_id", "int PRIMARY KEY"),
    ("parent_tax_id", "int"),
    ("rank", "text"),
    ("embl_code", "text"),
    ("division_id", "int"),
    ("inherited_div_flag", "boolean"),
    ("genetic_code_id", "int"),
    ("inherited_GC_flag", "boolean"),
    ("mito_code_id", "int"),
    ("inherited_MGC_flag", "boolean"),
    ("GenBank_hidden_flag", "boolean"),
    ("hidden_subtree_root_flag", "boolean"),
    ("comments", "text"),
]

names_file_schema = [
    ("tax_id", "int"),
    ("name_txt", "text"),
    ("unique_name", "text"),
    ("name_class", "text"),
]

merge_file_schema = [("old_tax_id", "int"), ("new_tax_id", "int")]

# Get a logical index of which fields are being kept for each table
nodes_drop_index = [
    True if vl[0] in nodes_drop_fields else False for vl in nodes_file_schema
]
names_drop_index = [
    True if vl[0] in names_drop_fields else False for vl in names_file_schema
]
merge_drop_index = [
    True if vl[0] in merge_drop_fields else False for vl in merge_file_schema
]

# Create the final schema for all three tables
nodes_output_schema = ", ".join(
    [
        " ".join(val)
        for val, drop in zip(nodes_file_schema, nodes_drop_index)
        if not drop
    ]
)
nodes_output_schema = f"CREATE TABLE nodes ({nodes_output_schema})"

nodes_insert_placeholders = ",".join(
    ["?"] * (len(nodes_drop_index) - sum(nodes_drop_index))
)
nodes_insert_statement = f"INSERT INTO nodes VALUES ({nodes_insert_placeholders})"

names_output_schema = ", ".join(
    [
        " ".join(val)
        for val, drop in zip(names_file_schema, names_drop_index)
        if not drop
    ]
)
names_output_schema = f"CREATE TABLE names ({names_output_schema})"

names_insert_placeholders = ",".join(
    ["?"] * (len(names_drop_index) - sum(names_drop_index))
)
names_insert_statement = f"INSERT INTO names VALUES ({names_insert_placeholders})"

merge_output_schema = ", ".join(
    [
        " ".join(val)
        for val, drop in zip(merge_file_schema, merge_drop_index)
        if not drop
    ]
)
merge_output_schema = f"CREATE TABLE merge ({merge_output_schema})"

merge_insert_placeholders = ",".join(
    ["?"] * (len(merge_drop_index) - sum(merge_drop_index))
)
merge_insert_statement = f"INSERT INTO merge VALUES ({merge_insert_placeholders})"

# Create the output file and turn off safety features for speed
con = sqlite3.connect(db_file)
con.execute("PRAGMA synchronous = OFF")

# Create the tables
con.execute(nodes_output_schema)
con.execute(names_output_schema)
con.execute(merge_output_schema)
con.commit()

# Import data from the node data
with open(nodes) as bbn:

    # The files are tab delimited but the quoting is sometimes unclosed,
    # so turning off quoting - includes quotes in the fields where present
    bb_reader = csv.reader(bbn, delimiter="|", quoting=csv.QUOTE_NONE)

    for row in bb_reader:
        row = [val.strip() for val, drp in zip(row, nodes_drop_index) if not drp]

        con.execute(nodes_insert_statement, row)

# Import data from the names data
with open(names) as bbn:

    # The files are tab delimited but the quoting is sometimes unclosed,
    # so turning off quoting - includes quotes in the fields where present
    bb_reader = csv.reader(bbn, delimiter="|", quoting=csv.QUOTE_NONE)

    for row in bb_reader:
        row = [val.strip() for val, drp in zip(row, names_drop_index) if not drp]

        con.execute(names_insert_statement, row)

# Import data from the merge data
with open(merge) as bbn:

    # The files are tab delimited but the quoting is sometimes unclosed,
    # so turning off quoting - includes quotes in the fields where present
    bb_reader = csv.reader(bbn, delimiter="|", quoting=csv.QUOTE_NONE)

    for row in bb_reader:
        row = [val.strip() for val, drp in zip(row, merge_drop_index) if not drp]

        con.execute(merge_insert_statement, row)

# Create the indices
con.execute("CREATE INDEX node_id ON nodes (tax_id);")
con.execute("CREATE INDEX all_names ON names (name_txt);")
con.execute("CREATE INDEX id_name_class ON names (tax_id, name_class);")
con.execute("CREATE INDEX merged_id ON merge (old_tax_id);")
con.commit()

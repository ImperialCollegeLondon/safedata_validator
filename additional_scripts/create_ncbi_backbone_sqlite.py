"""
This script is used to build an SQLite3 database file containing the contents of
the NCBI taxonomy. It is not currently provided as a part of the package because
the details of file names etc. change and so it is provided as a recipe.
"""

import csv
import sqlite3

# RECKON I MIGHT NEED TO DEFINE MULTIPLE FILE SCHEMA, SEE HOW I GET ON
# DOES "text[]" MEAN SOMETHING RELATIVE TO "text"?
# DO I NEED TO INSERT NONES LIKE DAVID DID?
# PROBABLY SHOULD EXPAND THIS TO MAKE THREE SEPERATE DATABASES

db_file = 'local_db/ncbi/ncbi_nodes.sqlite3'
drop_fields = ['embl_code', 'division_id', 'inherited_div_flag', 'genetic_code_id',
               'inherited_GC_flag', 'mito_code_id', 'inherited_MGC_flag',
               'GenBank_hidden_flag', 'hidden_subtree_root_flag']
nodes = 'local_db/ncbi/nodes.dmp'
names = 'local_db/ncbi/names.dmp'
merged = 'local_db/ncbi/merged.dmp'

# Create the schema for the table, using drop fields to remove unwanted fields
# in the schema and the data tuples. The file_schema list describes the full set
# of fields provided by NCBI
file_schema = [('tax_id', 'int PRIMARY KEY'),
                ('parent_tax_id', 'int'),
                ('rank', 'text'),
                ('embl_code', 'text'),
                ('division_id', 'int'),
                ('inherited_div_flag', 'boolean'),
                ('genetic_code_id', 'int'),
                ('inherited_GC_flag', 'boolean'),
                ('mito_code_id', 'int'),
                ('inherited_MGC_flag', 'boolean'),
                ('GenBank_hidden_flag', 'boolean'),
                ('hidden_subtree_root_flag', 'boolean'),
                ('comments', 'text')]

# Get a logical index of which fields are being kept
drop_index = [True if vl[0] in drop_fields else False for vl in file_schema]

# Create the final schema for the backbone table, including a new field to show
# deleted taxa and insert statements for both the backbone and deleted files
output_schema = ', '.join([' '.join(val) for val, drop in zip(file_schema, drop_index) if not drop])
output_schema = f"CREATE TABLE nodes ({output_schema})"

insert_placeholders = ','.join(['?'] * (len(drop_index) - sum(drop_index)))
insert_statement  = f"INSERT INTO nodes VALUES ({insert_placeholders})"

# Create the output file and turn off safety features for speed
con = sqlite3.connect(db_file)
con.execute('PRAGMA synchronous = OFF')

# Create the table
con.execute(output_schema)
con.commit()

# Import data from the node data
with open(nodes) as bbn:

    # The files are tab delimited but the quoting is sometimes unclosed,
    # so turning off quoting - includes quotes in the fields where present
    bb_reader = csv.reader(bbn, delimiter='|', quoting=csv.QUOTE_NONE)

    for row in bb_reader:
        row = [None if val == '\\N' else val
               for val, drp in zip(row, drop_index) if not drp]

        con.execute(insert_statement, row)

# Create the indices
# con.execute('CREATE INDEX backbone_name_rank ON backbone (canonical_name, rank);')
con.execute('CREATE INDEX node_id ON nodes (tax_id);')
con.commit()

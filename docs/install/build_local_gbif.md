# Local GBIF backbone database

GBIF provide the backbone taxonomy data as a SQL file defining the table
structure and a tab delimited text file containing the backbone data. To use
this data as a local data resource for taxon validation, the
`safedata_validator` package requires it to be loaded into a SQLite3 database.

You may need to download and install SQLite3. You want the 'bundle' installer
for your operating system.

[https://www.sqlite.org/download.html](https://www.sqlite.org/download.html)

Next, go to this page:

[https://hosted-datasets.gbif.org/datasets/backbone/current/](https://hosted-datasets.gbif.org/datasets/backbone/current/)

The file names change from time to time, but from this directory you will
need two files. Currently these are `simple.txt.gz` and `deleted.txt.gz`.
These files contain the data for a [simplified version](https://hosted-datasets.gbif.org/datasets/backbone/README.html)
of the GBIF backbone, including taxa that have been deleted from the GBIF
backbone.

Those files are both dumps from a PostGRESQL database, and the definition
(schema) for the resulting table can be found
[here](https://raw.githubusercontent.com/gbif/checklistbank/master/checklistbank-mybatis-service/src/main/resources/backbone-ddl.sql)

There are a number of steps needed to convert this data into a SQLite3 database.
For convenience, this conversion has been automated into a Python script
(`scripts/create_gbif_backbone_sqlite.py`). This may need updating if file
names and structures change, but the basic process should work.

* The table contains some very long fields (notably `name_published_in`) that
  include a lot of quotes and add to the file size, but are not used by the
  package. This field is therefore dropped.

* The file contains a lot of `\N` values, which is a PostgreSQL symbol for a
   null field. SQLite3 will treat these as values and so they need to be
   converted to a `null` value.

* The speed of the package is much improved by building covering indices to
   speed up the two kinds of searches used by `safedata_validator`:

      1.  searches on the canonical name  and rank of a taxon and 
      2.  searches on the taxon id.

Once you have an SQLite3 backbone database, you will then need to edit the
`gbif_database` entry in your [configuration file](configuration.md) to provide
the path to your new SQLite file or provide the path as an argument.

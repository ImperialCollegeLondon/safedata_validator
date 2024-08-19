# Local GBIF backbone database

## The data source

GBIF maintain a detailed database that is used to provide the hierarchical backbone
taxonomy underpinning GBIF species observations. They provide multiple versions of the
dataset - at roughly one year to six month intervals - that are identified with a date
timestamp. All versions are freely available from GBIF at the link below, and the
`current` folder is a shortcut to the most recent version.

[https://hosted-datasets.gbif.org/datasets/backbone/current/](https://hosted-datasets.gbif.org/datasets/backbone/current/)

To use the GBIF backbone taxonomy as a local data resource for taxon validation, the
`safedata_validator` package requires a version of the database  to be built into a
SQLite3 database.

## Building the local GBIF database

The `safedata_build_local_gbif` command line tool is used to automatically download and
build the required file. The command line help is shown below - note that a particular
version can be selected by supplying the version date timestamp from the page above.

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_build_local_gbif.txt"
%}
```

You will need to provide an output file path for the database and then use the command:

```sh
safedata_build_local_gbif outfile
```

This should result in the following output:

```txt
- Downloading GBIF data to: /path/to/tempdir
    - Checking for version with timestamp 2023-08-28
    - Downloading simple.txt.gz
100%|████████████████████████████████████████████| 466M/466M [00:05<00:00, 92.2MB/s]
    - Downloading simple-deleted.txt.gz
100%|████████████████████████████████████████████| 90.8M/90.8M [00:00<00:00, 96.2MB/s]
- Building GBIF backbone database in: /path/to/outdir/gbif_backbone_2023-08-28.sqlite
    - Timestamp table created
    - Backbone table created
    - Adding core backbone taxa
7746724it [03:31, 36688.61it/s]
    - Adding deleted taxa
1711901it [00:42, 40470.91it/s]
    - Creating database indexes
    - Removing downloaded files
```

Once you have an SQLite3 backbone database, you will then need to edit the
`gbif_database` entry in your [configuration file](configuration.md) to provide
the path to your new SQLite file.

## Build process overview

From the archive directory for the version, the database is built from two files:
`simple.txt.gz` and `deleted.txt.gz`. These files contain the data for a [simplified
version](https://hosted-datasets.gbif.org/datasets/backbone/README.html) of the GBIF
backbone, including taxa that have been deleted from the GBIF backbone.

Those files are both dumps from a PostGRESQL database, and the definition
(schema) for the resulting table can be found
[here](https://raw.githubusercontent.com/gbif/checklistbank/master/checklistbank-mybatis-service/src/main/resources/backbone-ddl.sql)

There are a number of steps needed to convert this data into a SQLite3 database, but the
basic process is:

* The `simple` file contains some very long fields (notably `name_published_in`) that
  include a lot of quotes and add to the file size, but are not used by the
  package. This field is therefore dropped.

* Both files contain a lot of `\N` values, which is a PostgreSQL symbol for a
   null (empty) field. SQLite3 will treat these as values and so they need to be
   converted to a `null` value.

* The main `backbone` table is created and then data from both files are inserted as the
  data rows.

* A `timestamp` table is inserted to record the timestamp of the database version used
  to build the local file.

* The speed of the package is much improved by building covering indices to
   speed up the two kinds of searches used by `safedata_validator`:

   1. searches on the canonical name  and rank of a taxon and
   2. searches on the taxon id.

# Local NCBI taxonomy database

## The data source

The NCBI maintain a taxonomy database, providing a taxonomic hierarchy for the sequences
stored in the NCBI system. The database is continually updated, but they provide an
`archive` folder containing roughly monthly snapshots of the taxonomy database. We
recommend using these archive snapshots, because they provide a timestamped and
reproducible version of the underlying data. The available versions can be seen here:

[](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/)

To use the NCBI taxonomy as a local data resource for taxon validation, the
`safedata_validator` package requires a version of the database to be built into a
SQLite3 database.

## Building a local copy

The `safedata_build_local_ncbi` command line tool is used to automatically download and
build the required file. The command line help is shown below - note that a particular
version can be selected by supplying the version date timestamp from the page above.

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_build_local_ncbi.txt"
%}
```

You will need to provide an output directory for the database and then use the command:

```sh
safedata_build_local_ncbi outdir
```

This should result in the following output:

```txt
- Downloading NCBI data to: /path/to/tempdir
    - Connecting to FTP server
    - Using most recent archive: 2023-11-01
    - Downloading taxonomy to: /path/to/tempdir/taxdmp_2023-11-01.zip
100%|███████████████████████████████████████████| 60.5M/60.5M [00:03<00:00, 16.9MB/s]
- Building GBIF backbone database in: /path/to/outdir/ncbi_taxonomy_2023-11-01.sqlite
    - Timestamp table created
    - Creating nodes table
    - Populating nodes table from nodes.dmp
2535034it [00:27, 93695.95it/s]
    - Creating names table
    - Populating names table from names.dmp
3914203it [00:31, 123281.67it/s]
    - Creating merge table
    - Populating merge table from merged.dmp
74501it [00:00, 155486.74it/s]
    - Creating unique ranks table
    - Creating database indexes
    - Removing downloaded archive
```

Once you have an SQLite3 backbone database, you will then need to edit the
`gbif_database` entry in your [configuration file](configuration.md) to provide
the path to your new SQLite file.

## Build process overview

* The command line tool downloads the `taxdmp` zipfile archive for the requested
  timestamp.
* This file contains data from a number of tables stored as DMP files, but only the
  `nodes.dmp`, `names.dmp`, and `merged.dmp` are required.  These describe the taxonomic
  nodes, the names assigned to specific nodes, and superseded nodes that have now been
  merged into other nodes, respectively.
* The tool automatically extracts this data from within the zip archive.
* The nodes table contains a lot of un-needed fields - we therefore drop all fields
  other than `tax_id`, `parent_tax_id`, `rank`, and `comments`.
* The names table includes a field (`unique_name`) that gives a unique version of every
  name. Given that `safedata_validator` has built in handling of ambiguous names and
  that these unique names are extremely specific (e.g. `Bacteria <bacteria>`), we do not
  believe this is a useful field for our purpose. It is therefore dropped.
* A single database is then generated with three tables: `nodes`, `names`, `merged`.
* A `timestamp` table is created to hold the version timestamp.
* A `unique_ncbi_ranks` table is created to hold the names of all ranks used in the NCBI
  database.
* The speed of the package is much improved by building covering indices to
   speed up the four kinds of searches used by `safedata_validator`:

  1. searches for specific taxon ids in the nodes table,
  2. searches for taxon names in the names table,
  3. searches for scientific names in the names tables and
  4. searches for specific taxon ids in the merged table.

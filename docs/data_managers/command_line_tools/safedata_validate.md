# Using the SAFE Dataset Checker

The  `safedata_validator` package contains Python code to validate files
containing data using `safedata` formatting  and report on any problems. The code
validates:

1. The data submission formatting of the file.
1. All taxonomic names against the [GBIF taxonomy
   database](../install/taxonomic_validation.md).
1. All location names against a locations gazetteer.

The package can be imported in Python for use within other frameworks. However, the
package also provides a command line tool, that allows it to be  installed and run as a
standalone application, for example by data managers or individual researchers.

## Configuring data resources

The `safedata_validator` package requires external data resources to validate
both dataset locations and taxa. You will need to [create a configuration
file](../install/configuration.md) to set `safedata_validator` up to find those
resources.

Note that a key resource - the GBIF taxonomy database - require a local SQLite3 database
containing the core data from this database. This is a relatively large file (~ 2GB in
total). The package provides a command
([safedata_build_local_gbif](../install/build_local_gbif.md) to download and build this
database, and the path to this file can then be included in the configuration.

Note that you cannot validate DOIs without an internet connection, but this is optional.

### GBIFTaxa

To validate taxonomic information against GBIF, you will need to download a copy of the
**GBIF backbone taxonomy** and build a SQLite3 database from it. The package provides a
template Python script to do this. If you are happy with running Python scripts, then it
is not particularly hard and is described in detail
[here](../install/build_local_gbif.md). The resulting database file is around 1.6GB, so
you'll need file space!

## Using `safedata_validate`

Once you have setup and configured `safedata_validator`, the usage instructions are
below:

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_validate.txt"
%}
```

Essentially, you should now be able to do:

```bash
safedata_validate MyDataset.xlsx
```

The program will then validate the input dataset, printing information about the
validation process and any errors in the dataset as it goes. When a dataset passes
validation, a JSON file containing the metadata for the validated dataset will be
created.

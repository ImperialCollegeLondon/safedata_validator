# Using the SAFE Dataset Checker

The  `safedata_validator` package contains Python code to validate files
containing data using SAFE data formatting  and report on any problems. The code
validates:

  1. The data submission formatting of the file.
  1. All taxonomic names against either the [GBIF taxonomy
     database](../install/gbif_validation.md) or the [NCBI taxonomy
     database](../install/ncbi_validation.md).
  1. All location names against a locations gazetteer.

This package can be used within other frameworks, such as the [SAFE Project
website](https://safeproject.net/datasets/view_datasets). However, the package also
provides a command line tool, that allows it to be  installed and run as a standalone
application, for example by data managers or individual researchers.

## Configuring data resources

The `safedata_validator` package requires external data resources to validate
both dataset locations and taxa. You will need to [create a configuration
file](../install/configuration.md) to set `safedata_validator` up to find those
resources.

Note that two key resources - the GBIF and NCBI taxonomy databases - require local
SQLite3 databases containing the core data from those databases. These are relatively
large files (~ 2GB in total). The package provides two commands
([safedata_build_local_gbif](../install/build_local_gbif.md) and
[safedata_build_local_ncbi](../install/build_local_ncbi.md)) to download and build
these databases, and the path to those files can then be included in the configuration.

Note that you cannot validate DOIs without an internet connection, but
this is optional.

### GBIFTaxa

If you  want to speed up taxon checking and allow offline use then you will need to
download a copy of the **GBIF backbone taxonomy** and build a SQLite3 database from it.
The package provides a template Python script to do this. If you are happy with running
Python scripts, then it is not particularly hard and is described in detail
[here](../install/build_local_gbif.md). The resulting database file is around 1.6GB, so
you'll need file space!

### NCBITaxa

Taxon checking against NCBI can be similarly sped up by downloading a copy of the **NCBI
database** and building a SQLite3 database from it. Using a local database is
substantially faster than using the online NCBI Entrez tools, which has a built-in rate
limitation. Instructions on how to construct the local database are given
[here](../install/build_local_ncbi.md). Again, the resulting database is large (~600 MB)
so you will need to ensure you have sufficient file space!

## Using `safedata_validate`

Once you have setup and configured `safedata_validator`, the usage instructions are
below:

```sh
{!docs/command_line_tools/command_line_usage/safedata_validate.txt!}
```

Essentially, you should now be able to do:

```bash
safedata_validate MyDataset.xlsx
```

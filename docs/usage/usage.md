# Using the SAFE Dataset Checker

!!! Warning
    TODO - Outdated

The  `safedata_validator` package contains Python code to validate files
containing data using SAFE data formatting  and report on any problems. The code
validates:

  1. The data submission formatting of the file.
  1. All taxonomic names against either the [GBIF taxonomy
     database](../install/gbif_validation.md) or the [NCBI taxonomy
     database](../install/ncbi_validation.md).
  1. All location names against the SAFE Gazetteer.

This package is used to validate datasets submitted online to a SAFE data
management website, such as the  [SAFE Project
website](https://safeproject.net/datasets/submit_dataset). However, it can also
be installed and run independently, for example by data managers or individual
researchers.

The code is open source Python and  is maintained on
[GitHub](https://github.com/ImperialCollegeLondon/safedata_validator) but  can
also be installed using  [PyPI](https://pypi.org/project/safedata-validator).
See the [installation notes](../install/install.md) for setup instructions.

The package provides a command line program `safedata_validate`. The usage
instructions are below but you will also need to provide links to some external
data resources used in location and taxon validation.

```sh
{!docs/command_line_usage/safedata_validate.txt!}
```

## Data resources

The `safedata_validator` package requires external data resources to validate
both dataset locations and taxa. The package supports online resources for both
locations and taxa, which is the easiest option for users to get up and running.
The online [GBIF  Search API](https://www.gbif.org/developer/species) is used by
default but the package does need a web service providing valid location data.

For example, the SAFE Project  website provides an API endpoint returning
location data. Using this API and the default online GBIF validation, the
following command will validate `MyData.xlsx`:  

```sh
safedata_validate MyData.xlsx -l https://www.safeproject.net/api/validator_locations
```

This is considerably easier for most users but it can be rather slow and
requires an internet connection. If you are want to improve the speed of
`safedata_validator` for frequent use  or need to be able to use it offline,
then you can provide local copies of the data resources. Note that you cannot
validate DOIs without an internet connection, but this is optional.

The locations of these resources are set by command line arguments shown above
but  can also be set in a [configuration file](usage.md#configuration-file) for
repeated use.

### Locations

Locations are validated against a set of known location names and possible
aliases for those names. The data resource providing this information is set
with the `location` argument. This can either be a link to a web service as
shown above or a static local JSON file to provide faster and offline use:

```sh
safedata_validate MyData.xlsx -l /path/to/validator_locations.json
```

### GBIFTaxa

If you  want to speed up taxon checking and allow offline use then you will need
to download a copy of the backbone taxonomy and build a SQLite3 database from
it. Using a local database is  much faster than using the GBIF API online. This
isn't particularly hard and is described in detail [here](../install/build_local_gbif.md),
but the resulting database is around 1.6GB, so you'll need file space!

Once you have this file, you can use it like this:

```sh
safedata_validate MyData.xlsx -g /path/to/gbif_backbone.sqlite
```

### NCBITaxa

Taxon checking against NCBI can be similarly sped up by downloading a copy of the online database and building a SQLite3 database from it. Using a local database is substantially faster than using the online NCBI Entrez tools. Instructions on how to construct the local database are given [here](../install/build_local_ncbi.md), again the resulting database is large so you'll need to ensure you have sufficient file space!

Once you have this file, you can use it like this:

    TODO - INSERT USAGE HERE

### Fully offline use

If you've done both the above steps then the following example would validate a
file using both local data resources, and won't need the internet at all.

```sh
safedata_validate MyData.xlsx -g /path/to/gbif_backbone.sqlite \
    -l /path/to/validator_locations.json
```

In both cases, validation can now simply use:

```sh
safedata_validate MyData.xlsx
```

If you do provide command line arguments, they will override anything set in the
configuration.

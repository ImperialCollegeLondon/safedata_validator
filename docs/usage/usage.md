# Using the SAFE Dataset Checker

The  `safedata_validator` package contains Python code to validate files containing data using SAFE data formatting  and report on any problems. The code validates:

  1. The data submission formatting of the file.
  1. All taxonomic names against the [GBIF taxonomy database](../install/gbif_validation.md).
  1. All location names against the SAFE Gazetteer.

This package is used to validate datasets submitted online to a SAFE data management website, such as the  [SAFE Project website](https://safeproject.net/datasets/submit_dataset). However, it can also be installed and run independently, for example by data managers or individual researchers.

The code is open source Python and  is maintained on [GitHub](https://github.com/ImperialCollegeLondon/safedata_validator) but  can also be installed using  [PyPI](https://pypi.org/project/safedata-validator). See the [installation notes](../install/install.md) for setup instructions. 

The package provides a command line program `safedata_validate`. The usage instructions are below but you will also need to provide links to some external data resources used in location and taxon validation.

    usage: safedata_validate [-h] [-l LOCATIONS_JSON]
                                   [--gbif_database GBIF_DATABASE]
                                   [--validate_doi]
                                   fname

    This program validates an Excel file formatted as a SAFE dataset. As it runs,
    it outputs a report that highlights any problems with the formatting. Much of
    the validation is to check that the data meets our metadata standards and is
    internally consistent. However, it uses external sources to perform validation
    in three areas.

    1. Taxon validation. The program validates taxonomic names against the GBIF
    taxonomy backbone. By default, it uses the GBIF web API to validate names,
    but can also use a local copy of the backbone provided in a sqlite database:
    this will work offline and is much faster but requires some simple setup.

    2. Location names. The program also validate sampling location names against
    the SAFE gazeteer. By default, this is loaded automatically from the SAFE
    website so requires an internet connection, but a local copy can be provided
    for offline use.

    3. DOI checking. Optionally, the program will validate any DOIs provided as
    having used the database. This requires a web connection and cannot be
    performed offline.

    positional arguments:
      fname                 Path to the Excel file to be validated.

    optional arguments:
      -h, --help            show this help message and exit
      -l LOCATIONS_JSON, --locations_json LOCATIONS_JSON
                            Path to a locally stored json file of valid location
                            names
      -g GBIF_DATABASE, --gbif_database GBIF_DATABASE
                            The path to a local sqlite database containing the
                            GBIF taxonomy backbone.
      --validate_doi        Check the validity of any publication DOIs, provided
                            by the user. Requires a web connection.

# Data resources

The `safedata_validator` package requires external data resources to validate both dataset locations and taxa. The package supports online resources for both locations and taxa, which is the easiest option for users to get up and running. The online [GBIF  Search API](https://www.gbif.org/developer/species) is used by default but the package does need a web service providing valid location data.

For example, the SAFE Project  website provides an API endpoint returning location data. Using this API and the default online GBIF validation, the following command will validate `MyData.xlsx`:  

    safedata_validate MyData.xlsx -l https://www.safeproject.net/api/validator_locations

This is considerably easier for most users but it can be rather slow and requires an internet connection. If you are want to improve the speed of `safedata_validator` for frequent use  or need to be able to use it offline, then you can provide local copies of the data resources. Note that you cannot validate DOIs without an internet connection, but this is optional.

The locations of these resources are set by command line arguments shown above but  can also be set in a [configuration file](usage.md#configuration-file) for repeated use. 

## Locations 

Locations are validated against a set of known location names and possible aliases for those names. The data resource providing this information is set with the `location` argument. This can either be a link to a web service as shown above or a static local JSON file to provide faster and offline use:

    safedata_validate MyData.xlsx -l /path/to/validator_locations.json

## Taxa

If you  want to speed up taxon checking and allow offline use then you will need to download a copy of the backbone taxonomy and build a SQLite3 database from it. Using a local database is  much faster than using the GBIF API online. This isn't particularly hard and is described in detail [here](build_local_gbif.md), but the resulting database is around 1.6GB, so you'll need file space! 

Once you have this file, you can use it like this:

    safedata_validate MyData.xlsx -g /path/to/gbif_backbone.sqlite

## Fully offline use

 If you've done both the above steps then the following example would validate a file using both local data resources, and won't need the internet at all.

    safedata_validate MyData.xlsx -g /path/to/gbif_backbone.sqlite \
        -l /path/to/validator_locations.json



In both cases, validation can now simply use:

    safedata_validate MyData.xlsx

If you do provide command line arguments, they will override anything set in the configuration.

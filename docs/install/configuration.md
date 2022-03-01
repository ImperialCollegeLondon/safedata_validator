# Configuring the `safedata_validator` package

The `safedata_validator` package needs to be configured to use specific
resources for data validation. This is done using a configuration file in the
[INI file format](https://en.wikipedia.org/wiki/INI_file).

## Configuration file format

The file structure for the configuration file is a simple text file containing
the details below:

```
locations = /path/to/locations.json
gbif_database = /path/to/local/backbone.sqlite3
[extents]
temporal_soft_extent = 2002-02-01, 2030-02-01
temporal_hard_extent = 2002-02-01, 2030-02-01
latitudinal_hard_extent = -90, 90
latitudinal_soft_extent = -4, 2
longitudinal_hard_extent = -180, 180
longitudinal_soft_extent = 110, 120
[zenodo]
community_name = safe
use_sandbox = True
zenodo_sandbox_api = https://sandbox.zenodo.org
zenodo_sandbox_token = xyz
zenodo_api = https://api.zenodo.org
zenodo_token = xyz
```

Only the `locations` entry is really required, pointing to a local copy of a
SAFE format gazetteer file in JSON format. An example of this file can be
downloaded from:

[https://www.safeproject.net/api/validator_locations]("https://www.safeproject.net/api/validator_locations")

The `gbif_database` value can be left empty (`gbif_database = `), in which case
the package will use the online GBIF API for taxon validation. Alternatively,
the `gbif_database` entry can be used to point to a local copy of the GBIF
backbone database, which is _very_ much faster but requires the database to be
built (see [here](build_local_gbif))

### Extents

The `safedata_validator` package tracks the geographic and temporal extents of
datasets, which are needed to generate Gemini 2 metadata for a dataset. These
extents can have both hard and soft bounds set: values outside of hard bounds
will cause a validation error, where values outside of soft bounds only result
in a validation warning.

If missing from the configuration file, then the only bounds used are hard bounds
on geographic coordinates:

* `latitudinal_hard_extent`: (-90, 90)
* `longitudinal_hard_extent`: (-180, 180)


### Zenodo

!!! Warning "In development"

    These features are in development and are not yet functional.


The `safedata_validator` package provides functionality for publishing validated
datasets to the Zenodo data repository. In order to do this, the Zenodo
community name and access keys need to be stored in the configuration file.

## Configuration file location

The path to a configuration file can be provided when running validation but, to
avoid having to do this,  `safedata_validator` looks in specific locations for
configuration files. Conventions for config file locations differ across
operating systems and we use the conventions used by the `appdirs` package.

The `safedata_validator` package will look for both user and site configuration
files. Site configurations allow a data manager to set up a specific machine
with data resources for all users. If both are present, the user configuration
is used. The configuration file must be called `safedata_validator.cfg` and the
user and site config folders are:

### Mac OS X

```sh
/Users/username/Library/Application Support/safedata_validator/safedata_validator.cfg
/Library/Application Support/safedata_validator/safedata_validator.cfg
```

### Windows 

The repeated directory name is not an error:

```sh
C:\Users\username\AppData\Local\safedata_validator\safedata_validator\safedata_validator.cfg
```

### Linux

```sh
/home/username/.config/safedata_validator/safedata_validator.cfg
/etc/xdg/safedata_validator/safedata_validator.cfg
```


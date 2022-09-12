# Configuring the `safedata_validator` package

The `safedata_validator` package needs to be configured to use specific
resources for data validation. This is done using a configuration file in the
[INI file format](https://en.wikipedia.org/wiki/INI_file).

## Configuration file format

The file structure for the configuration file is a simple text file containing
the details below:

```ini
locations = /path/to/locations.json
gbif_database = /path/to/local/backbone.sqlite3
ncbi_database = /path/to/local/ncbi_database.sqlite3
[extents]
temporal_soft_extent = 2002-02-02, 2030-01-31
temporal_hard_extent = 2002-02-01, 2030-02-01
latitudinal_hard_extent = -90, 90
latitudinal_soft_extent = -4, 6
longitudinal_hard_extent = -180, 180
longitudinal_soft_extent = 110, 120
[zenodo]
community_name = safe
use_sandbox = true
zenodo_sandbox_api = https://sandbox.zenodo.org
zenodo_api = https://api.zenodo.org
zenodo_sandbox_token = xyz
zenodo_token = abc
contact_name = The SAFE Project
contact_affiliation = Imperial College London
contact_orcid = 0000-0003-3378-2814
[metadata]
api = https://safeproject.net
token = xyz
[ncbi]
api_key = xyz
```

### Locations

Locations are validated against a set of known location names and possible
aliases for those names. The data resource providing this information is set
with the `location` configuration, providing to local JSON formatted file containing
known locations and location aliases. An example of this file format can be
downloaded from:

[https://www.safeproject.net/api/validator_locations](https://www.safeproject.net/api/validator_locations)

### GBIF database

The `gbif_database` value can be left empty (`gbif_database =`), in which case
the package will use the online GBIF API for taxon validation. Alternatively,
the `gbif_database` entry can be used to point to a local copy of the GBIF
backbone database, which is _very_ much faster but [requires the database to be
built](build_local_gbif.md).

### NCBI database

The `ncbi_database` value can also be left empty (`ncbi_database =`), in which case the
package will access the online NCBI taxonomy database using the [Entrez
system](https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html). Alternatively, the
`ncbi_database` entry can be used to point to a local copy of the NCBI database, which
is _very_ much faster but again [requires the database to be built](build_local_ncbi.md).

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

The `safedata_validator` package provides functionality for publishing validated
datasets to the Zenodo data repository. In order to do this, the Zenodo
community name and access keys need to be stored in the configuration file.

On Zenodo datasets are generally uploaded to specific communities. The community to
upload the dataset to is indicated in the configuration file by the `community_name`
option. In our case it is set to upload datasets to the [SAFE project Zenodo
community](https://zenodo.org/communities/safe). Uploading to a community requires that
you hold proper admin credentials for said community. These are provided in the form of
a `zenodo_token`. It is also important to provide the contact details of the community
when interacting with the api. These consist of the full community name, affiliation,
and orchid ID, and are provided by `contact_name`, `contact_affiliation`, and
`contact_orcid`, respectively.

When testing dataset uploads the normal Zenodo api should **not** be used. Instead the
sandbox api should be used. Documentation of this can be found
[here](https://developers.zenodo.org/#testing). The sandbox provides an identical
environment to the real Zenodo site, so that upload function outputs can be checked,
without adding unnecessary (or invalid) datasets to the official repository. Sandbox is
enabled by specifying `use_sandbox = true` in the configuration file. A different token
(`zenodo_sandbox_token`) has to be provided to upload to a community on the sandbox
site.

### Metadata server

Zenodo only allows a fairly limited amount of metadata to be stored for each dataset.
While this is completely adequate to describe the contents of a dataset, more extensive
metadata must be stored elsewhere if detailed searches within datasets are desired. In
our case, the SAFE project [website](https://safeproject.net) is used as a metadata
server, which provides an API allowing searches across all uploaded datasets. A further
package in the SAFE data ecosystem
([`safedata`](https://imperialcollegelondon.github.io/safedata/)) has been built, in
order to simplify the process of querying this API for end users.

### NCBI API

Remote NCBI validation is very slow due to requests being rate limited to a maximum of 3
requests per second. However, if an API key is provided the rate limit increases to 10
requests per second. An API key can be obtained by registering for an NCBI account and
then requesting one, instructions for this can be found
[here](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).
When this key is included in the configuration file `safedata_validator` will make use
of it for all remote validation requests.

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

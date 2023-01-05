# Configuring the `safedata_validator` package

The `safedata_validator` package needs to be configured to use specific
resources for data validation. This is done using a configuration file in the
[INI file format](https://en.wikipedia.org/wiki/INI_file).

## Configuration file format

The file structure for the configuration file is a simple text file containing
the details below:

```ini
gazetteer = /path/to/gazeteer.geojson
location_aliases = /path/to/location_aliases.csv
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
```

### Locations

Locations are validated against a set of known location names and possible aliases for
those names. The data resources providing this information are set with the `gazetteer`
and `location_aliases` configuration options.

The gazetteer file must be a [GeoJSON](https://geojson.org/) file containing a
collection of GIS features providing known locations for a project. The GIS features can
be simple point locations but the file can also include linestrings or polygons to
capture linear features like streams or transects and area features like quadrats. Each
feature will have a set of properties: `safedata_validator` only requires the `location`
property but it is fine for features to have other project specific properties. The
values of `location` **must be unique** within the gazetteer file.

Gazetteer locations must also be **persistent**: this core resource is used to link
location names to geographic data and removing locations from the gazetteer breaks that
link. It is fine to update other properties and the feature geometry.

A simple minimal example of GeoJSON file content is shown below, but this file format is
commonly used in GIS applications, so can be more easily created and edited using GIS
tools like [QGIS](https://qgis.org).

```json
{    "type": "FeatureCollection",
    "features": [
        {
            "type": "Feature",
            "properties": {
                "location": "location_name"
            },
            "geometry": {
                "type": "Point",
                "coordinates": [
                    117.586071,
                    4.710346
                ]
            }
        }
    ]
}
```

The location aliases file provides a way to handle commonly used alternative names for
sampling locations. The file format is a simple CSV file which must contain the field
headers shown below:

```csv
location,alias,zenodo_record_id
location_name,loc,NA
location_name_2,loc2,NA
```

Values added in the `alias` field will be automatically mapped to the given `location`.
Note that `safedata_validator` will _always_ warn about the use of location aliases.
Although the location aliases file must be provided, it is _fine_ to only include the
headers insist that users only use the canonical names from the gazetteer.

The `zenodo_record_id` field can be used to set a location alias for a particular
dataset - this allows new locations in published datasets to be associated with
canonical gazetteer locations. Alternatively, a new version of the dataset can be
published that updates the location data to use the gazetteer locations directly.

### GBIF database

The `gbif_database` entry contains the path to the required local copy of the GBIF
backbone database ([see here](build_local_gbif.md)).

### NCBI database

The `ncbi_database` entry contains the path to the required local copy of the NCBI
database, ([see here](build_local_ncbi.md)).

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

# The gazetteer

Datasets can [use their own unique
locations](../../data_providers/data_format/locations.md), but it is more generally
useful to have a centrally curated set of shared locations. These curated locations are
provided using **gazetteer** and **location alias** files, and known locations in
datasets are then validated against these resources. The data resources providing this
information are set with the `gazetteer` and `location_aliases` [configuration
options](configuration.md). This page exists to explains the purpose and contents of
these two files. See [here](gazetteer_creation.md) for a guide explaining how to create a
gazetteer.

## The gazetteer file

The gazetteer file must be a [GeoJSON](https://geojson.org/) file containing a
collection of GIS features providing known locations for a project.

* The file provides a set of GIS features, each of which is a known location used in
  datasets. These can be simple point features but the gazetteer can include mixed
  feature types to also include linestrings or polygons to capture linear features like
  streams or transects and area features like quadrats.

* Each feature will have a set of properties associated with the feature geometry.
  For `safedata_validator`, we only require the `location` property, which provides
  **a unique name** for each location within the gazetteer file. Although only the
  `location` property is required, it is fine for features to have other project
  specific properties.

* Gazetteer locations must also be **persistent**: this core resource is used to link
  location names to geographic data and removing locations from the gazetteer breaks
  that link. It is fine to update other properties and the feature geometry.

A simple minimal example of GeoJSON file content is shown below, but this file format is
commonly used in GIS applications, so can be more easily created and edited using GIS
tools like [QGIS](https://qgis.org).

```json
{"type": "FeatureCollection",
    "features": [
        {
            "type": "Feature",
            "properties": {
                "location": "Sample_Site_A"
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

## The location aliases file

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
Although the location aliases file must be provided, it is _fine_ to have an empty file
that only includes the column headings: users will then only be able to use the
canonical names from the gazetteer.

The `zenodo_record_id` field can be used to set a location alias for a particular
dataset - this allows new locations in published datasets to be associated with
canonical gazetteer locations. Alternatively, a new version of the dataset can be
published that updates the location data to use the gazetteer locations directly.

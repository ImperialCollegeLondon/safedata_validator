# Gazetteer creation guide

The basic reason that a gazetteer needs to be created is so that locations can be
referred to in datasets by name rather than having to provide the coordinates for every
single location. These locations can be any simple GIS feature, i.e. they can be a
simple point (for a sampling point), a line (for a transect), or a polygon (for large
plots/subplots). This gazetteer should be maintained by the data manager, this is to
prevent data providers from renaming existing points etc. If data providers sample from
genuinely new points they can either request that the data manager adds them to the
gazetteer or alternatively they can just add them to their datasets as [unique
locations](../../data_providers/data_format/locations.md).

The gazetteer file must be a [GeoJSON](https://geojson.org/) file. This file has to be
defined using the [WGS84 coordinate
system](https://en.wikipedia.org/wiki/World_Geodetic_System#WGS_84), with a local
transformation defined so that distances between points can be calculated. Each GIS
feature needs to have a name attribute specified. Other attributes can also be populated
if required, but it is important to note that these **will not** be searchable using the
metadata server. The metadata server can however provide an up to date copy of the
gazetteer that scripts or applications written by the data manager can make use of.

The gazetteer can either be generated using GIS software, or can be generated
programmatically. Guides for both approaches are provided below.

## Creation using GIS software

This can be done using whatever GIS software you like provided that the software can
export the final map as a GeoJSON file. The developers of the `safedata` system use
[QGIS](https://www.qgis.org/), but if you have a preferred GIS software already you
should stick with that. Sampling locations (i.e. points, lines and polygons) can all be
described by hand drawing in the software. However, if possible we would advise against
doing this, as hand drawing locations introduces the possibility of mis-drawing plots
(or transects etc). Instead we would recommend loading in existing shapes files for the
sampling locations, and then using the GIS software to verify that the locations are
correct and to add the relevant attributes. If shape files do not already exist they can
be generated from the `.gpx` files used in GPS units.

## Programmatic generation

If shape files either exist already or GPS points are available (in some format)
describing the locations, then the gazetteer can be generated programmatically (using
code). If you have the relevant files (or details of the points) and have at least some
familiarity with coding we would strongly advise taking this approach. This is to ensure
that the procedure used to create your gazetteer is reproducible, which makes it easier
to track down any errors that get introduced, as well as reducing the chance of entering
erroneous information to begin with.

We provide examples for `R` and `python` below as they are both highly popular
scientific programming languages. However, if you have a different preferred language
it's completely fine to use that instead.

### Generation using `R`

TODO - Option 1 - R + sf package
TODO - Simple example R script adding a point, line, and shape

### Generation using `python`

In python, we recommend that you use the [`shapely`
package](https://pypi.org/project/shapely/) to define your locations, the [`geopandas`
package](https://geopandas.org/en/stable/) to load in existing shape files and the
[`geojson` package](https://pypi.org/project/geojson/) to combine and export the
gazetteer.

The example below demonstrates how to manually create points, lines and polygons
manually using coordinates, how to load in shape files, and finally how to export
everything as combined `GeoJSON` file.

=== "python"

    ```python
    {%
    include "data_managers/install/create_gazetteer.py"
    %}
    ```

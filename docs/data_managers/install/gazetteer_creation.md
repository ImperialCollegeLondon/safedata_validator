# Gazetteer creation guide

TODO - Basic GIS overview, link names to GIS features (points, lines, polygons)

TODO - Explain Data managers need to curate this list (in a gazetteer), so people can
use names

TODO - Explain that Geojson file made by data manager, can use any GIS provided they can
export geojson.

TODO - They need to make a GIS file containing shapes plus names (and other feature
attributes if they want, SAFE example of fractal order)

TODO - Mention that additional attributes aren't currently searchable through the
safedata system, can point an app at the server

TODO - 1st approach: QGIS + hand drawing, combine existing shape files

TODO - Mention somewhere that Geojson file has to be in WGS84 coordinate system, then
setup to use a local transformation for distances

TODO - Second approach: Programmatic generation, explain general advantages
TODO - Option 1 - R + sf package
TODO - Simple example R script adding a point, line, and shape

TODO - Option 2 - Python + shapely
TODO - Simple example python script adding a point, line, and shape

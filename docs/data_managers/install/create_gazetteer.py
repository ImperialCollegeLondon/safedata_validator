"""Example script showing the basics of creating a gazetteer."""

import geojson
import geopandas as gpd
from shapely.geometry import LineString, Point, Polygon, mapping

# The easiest way to create features from coordinates is to create instances of
# Point, LineString and Polygon features using the Shapely package.

# A point sampling location with latitude = 4.250 and longitude = 117.900
sampling_point = Point(4.250, 117.900)

# A line (transect) created by specifying each sampling point along the transect
transect = LineString(
    [(4.250, 117.900), (4.251, 117.902), (4.252, 117.904), (4.253, 117.906)]
)

# A polygon or plot feature is created by specifying a point for each vertex of the
# polygon - the ring of the polygon is automatically closed if the first and last
# coordinates are not identical.
sampling_area = Polygon(
    [(4.250, 117.900), (4.255, 117.900), (4.255, 117.855), (4.250, 117.855)]
)

# Now convert the shapely points into GeoJSON features
feature_sampling_point = geojson.Feature(
    geometry=mapping(sampling_point), properties={"name": "point_A1"}
)
feature_transect = geojson.Feature(
    geometry=mapping(transect), properties={"name": "transect_B57"}
)
feature_sampling_area = geojson.Feature(
    geometry=mapping(sampling_area), properties={"name": "plot_C324"}
)

# Then combine them all the manually defined points into a single feature collection
manual_points = geojson.FeatureCollection(
    [feature_sampling_point, feature_transect, feature_sampling_area]
)

# We can also add features loaded from  an existing shapefile
shape_file = gpd.read_file("example_shape_file.shp")

# Convert the loaded in points to a geojson feature collection
shape_file_points = geojson.loads(shape_file.to_json())

# Finally combine the points from the different sources
all_points = geojson.FeatureCollection(
    manual_points["features"] + shape_file_points["features"]
)
# Then export everything as a geojson
with open("gazetteer.geojson", "w") as gazetteer_file:
    geojson.dump(all_points, gazetteer_file)

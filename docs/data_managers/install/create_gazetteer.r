library(sf)
library(geojsonio)

# Create points, lines, and polygon features
sampling_point <- st_point(c(117.900, 4.250))
transect <- st_linestring(
  cbind(c(117.900, 117.902, 117.904, 117.906), c(4.250, 4.251, 4.252, 4.253))
)
sampling_area <- st_polygon(list(cbind(
  c(117.900, 117.900, 117.855, 117.855, 117.900),
  c(4.250, 4.255, 4.255, 4.250, 4.250)
)))

# Combine into a single geometry collection with name attributes defined
# (coordinate reference system must be 4326 for GeoJSON)
geometry_column <- st_sfc(sampling_point, transect, sampling_area, crs = 4326)
manual_locations <- st_sf(
  name = c("point_A1", "transect_B57", "plot_C324"),
  geometry = geometry_column
)

# Load features from an existing shapefile
shape_file <- st_read("example_shape_file.shp")

# Combine features from different sources
all_locations <- rbind(manual_locations, shape_file)

# Export as GeoJSON
geojson_write(all_locations, file = "gazetteer.geojson", na = "null")

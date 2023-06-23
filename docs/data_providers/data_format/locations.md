# The Location worksheet

All locations used in your data worksheets need to be listed in this worksheet. By
**location**, we mean the common frequently used areas in which research has happened
for the project you are submitting data to. You might have more detail about the precise
place you worked in your dataset - great! - but using these known locations allows us to
get broad spatial data on sampling relatively simply.

So, we expect you'll have a set of location names in your data sheets, all of which
should appear in this worksheet. The worksheet must have the following structure

<!-- markdownlint-disable MD013 -->
| Location name | New | Latitude | Longitude  | Type       | WKT                                          |
| ------------- | --- | -------- | ---------- | ---------- | -------------------------------------------- |
| E_194         | No  |          |            |            |                                              |
| E_195         | No  |          |            |            |                                              |
| My_site_1     | Yes | 4.957721 | 117.776023 | POINT      | NA                                           |
| My_site_2     | Yes | NA       | NA         | POINT      | NA                                           |
| My_site_3     | Yes | NA       | NA         | POINT      | Point(117.7762 4.9576)                       |
| My_transect_1 | Yes | NA       | NA         | Linestring | Linestring(117.7762 4.9576, 117.7862 4.9676) |
<!-- markdownlint-enable MD013 -->

## Known locations

The location names are checked against a curated gazetteer of known locations maintained
for your `safedata` community. This will be constructed and maintained by your data
administrator, and is accessible through the `safedata` R package.

For example, in R:

```r
> library(jsonlite)
> locations <- fromJSON("https://www.safeproject.net/call/json/get_locations")
> str(locations) List of 1 $
locations: chr [1:2691] "SAFE_camp" "Flux_tower" "A_1" "A_2" ...
```

If you only use known locations, then you _only need to provide the  location name column_.

!!! Note "safedata at SAFE"
    Details of how to view and use the SAFE project specific gazetteer can be found
    [here](../../safedata_at_SAFE.md).

## New locations

If your data comes from genuinely new locations or uses a new sampling structure (e.g. a
grid or transect), then you can create new location names and include them in your
locations table. If they become commonly used, your data administrator may consider
adding them to their Gazetteer.

If you include new locations then you will need to include the following columns in your
Locations worksheet:

- **New**: This should simply contain Yes or No to show which rows contain new
  locations. You must enter 'No' for known locations and you cannot create a new
  location with a name that matches an existing location in the Gazetteer.

- **Type**: This is mandatory and is just an indication of the kind of sampling
  location: 'point', 'transect' or 'area' (or the more GIS names of 'point',
  'linestring' or 'polygon').

- **Latitude** and **Longitude**: these should provide GPS coordinates for the new site.
  These must be provided as decimal degrees (not degrees minutes and seconds) and please
  provide 6 decimal places in your coordinates. This level of precision is around ten
  centimetres and, although the GPS from the field is highly unlikely to be accurate to
  this level, we want to record as much sampling precision as possible.

  Note that you can provide a simple latitude and longitude for any location type: you
  might not have the coordinates of the whole transect but can give a start point. Any
  information is better than none but if you don't have any data, then you can enter
  NA.

- **WKT**: This optional field can be used to provide GIS geometry data for the location
  in the [well-known text
  format](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry).
  This is a good way to provide precise GIS data for a location.

Note that these extra columns can be left empty for known locations as in the example
above, but for new locations, please explicitly enter NA in these fields.

## Location aliases and extending the Gazetteer

In addition to the canonical location names, we also support a set of **location
aliases**. Although we prefer the canonical names to be used, aliases can be used
instead: an example is that "463" is accepted as an alias for the sampling point
"OG3_463".

The list of aliases also allows us to adopt new locations into the Gazetteer. If you
have sampled at new  locations that seem likely to be the focus of other research
projects then they may be added to your `safedata` community's Gazetteer so that other
researchers can use them as known locations. We can then use the location aliases to
record that 'New' locations in a published dataset have been adopted: for example that
the new location "river23" in dataset 123 is the same as the Gazetteer entry "RIVER23".

## GIS data

We prefer GIS information as Latitude and Longitude or WKT because we can use it to add
to the spatial index of the datasets, but you can also submit GIS files as [external
data files](other_formats.md) alongside your Excel data. You can even include metadata
about vector data attribute tables as part of the Excel file.

## My data doesn't include any locations

You don't **have** to include the Locations worksheet, although it would be very unusual
to omit it. Possible examples:

- You are working with lab data (and don't need to say where specimens came from in the
  field)
- You are collecting data haphazardly from across the landscape, for example tracking
  animal movements, and the data isn't tied to particular sampling locations. We would
  then want GPS data for each observation!

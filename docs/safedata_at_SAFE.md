# `safedata` at the SAFE project

TODO - Make this into a narrative rather than a bullet point list

- For SAFE project data access queries [contact us here](mailto:data@safeproject.net)
- Embargo period max of 2 years
- TODO - track down Metadata sever info
- Zenodo account (`the_safe_project`)
- This is particularly important for RCUK funded research, who have agreed to let us
    host all of the SAFE data under a common portal on the condition that Research
    Council funding is clearly acknowledged.

## Example bounds

|       |        |  |  |
|-------|--------|--|--|
| West  | 116.75 |  |  |
| East  | 117.82 |  |  |
| South | 4.50   |  |  |
| North | 5.07   |  |  |

The geographic bounds in the example cover Maliau, Danum and the SAFE Project
experimental site and surrounding area. While we would prefer something a bit
more precise, these are sufficient to provide a geographic extent for most work
at SAFE.

## Location links

You can look at the gazetteer webpage to see the available sites and to download
location data:

[https://www.safeproject.net/info/gazetteer](https://www.safeproject.net/info/gazetteer)

You can download the gazetteer locations from the URL above, but if you want to get a
list of valid location names for use in a program or script, then we also provide a web
service that returns a list of valid names as a JSON object:

[https://www.safeproject.net/call/json/get_locations](https://www.safeproject.net/call/json/get_locations)

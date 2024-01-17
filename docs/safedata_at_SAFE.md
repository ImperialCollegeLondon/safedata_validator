# `safedata` at the SAFE project

This document contains the specific details needed to use the `safedata_validator`
system to upload data for the SAFE project.

## Existing datasets from the SAFE project

Existing datasets from the SAFE project are stored on Zenodo as part of the
[`the_safe_project` Zenodo community](https://zenodo.org/communities/safe/).

## Data managers

The data managers for the SAFE project are at present David Orme and Jacob Cook. If you
have any data upload or dataset access queries we you can [contact us
here](mailto:data@safeproject.net).

## Data policies for the SAFE project

The maximum embargo period allowed for datasets submitted to the SAFE project is 2
years. If you have any queries about this please contact us.

Research Councils UK (RCUK) agreed to allow all data from the SAFE project to be hosted
under a single portal on the condition that Research Council funding is clearly
acknowledged. This means that funding details must be included if your research was in
anyway funding by a UK research council.

## SAFE metadata server

The SAFE project has historically used the [SAFE website](https://safeproject.net) as a
metadata server. It is not yet clear whether this will continue to be the case in the
long run or whether the project will transition to using the `safedata_server`
application to provide  a metadata server. If this transition occurs, we will inform
you of the new server to upload metadata to.

## SAFE project Geographic extent

Data collected as part of the SAFE project were gathered from the SAFE Project
experimental site and surrounding area, as well as Maliau, Danum. The bounds given in the
table below capture that entire area:

|       |        |  |  |
|-------|--------|--|--|
| West  | 116.75 |  |  |
| East  | 117.82 |  |  |
| South | 4.50   |  |  |
| North | 5.07   |  |  |

While we would prefer you to provide something a bit more precise, these are sufficient
to provide a geographic extent for most work at SAFE.

## The SAFE project gazetteer

The full set of locations for the SAFE project can be viewed using the gazetteer::

[https://www.safeproject.net/info/gazetteer](https://www.safeproject.net/info/gazetteer)

You can download the gazetteer locations from the URL above, but if you want to get a
list of valid location names for use in a program or script, then we also provide a web
service that returns a list of valid names as a JSON object:

[https://www.safeproject.net/call/json/get_locations](https://www.safeproject.net/call/json/get_locations)

# The Summary worksheet

This worksheet contains a simple set of rows describing the dataset and
identifying the spreadsheets that contain data tables. Each row is labelled on
the left in the first column and then the description data should be typed in
the columns to the right. The description below sets out the possible summary
metadata in blocks.

Some blocks of fields are mandatory (core, authors, worksheets, keywords) but
may include optional fields (such as author ORCID). Other blocks are optional,
but contain fields that may be mandatory if the block is used. We've tried to
make this as clear as possible below!

Some blocks allow multiple records (e.g. authors and data worksheets) with sets
of values in adjacent columns  but other blocks only allow a single record (e.g.
core fields and geographic extents).

<!-- markdownlint-disable MD033 -->
<style>

/*removing table headers and fixing cell widths so everything lines up.*/
thead {
  display: none;
}

table {
  table-layout: fixed;
}

tbody tr td {
  width: 14em;
  min-width: 14em;
  max-width: 14em;
  border: 1px solid lightgrey;
}

tbody tr td:first-child {
  width: 12em;
  min-width: 12em;
  max-width: 12em;
  font-weight: bold;
}

</style>
<!-- markdownlint-enable MD033 -->

## The core fields block

!!! Warning "Mandatory block"
    All fields are mandatory.

This block provides a set of core details for the dataset. You can only provide
a single value for each field.

* **Project ID**: This is only used if the organisation you are uploading to uses
  projects to group datasets. If they are, then you should obtain relevant project IDs
  from the organisation's data manager and add them in your Summary worksheet. This
  field can alternatively be given as `SAFE Project ID`. This alternative name is
  allowed solely in order to make the system backwards compatible with historic
  datasets. Using it is **not** the preferred option, particularly if you are submitting
  your data to someone other than the SAFE project.
* **Title**: This should be a short informative title for the dataset: it will
  be used as the public title for the dataset so make sure it is clear and
  grammatical!
* **Description**: This will be the public description of the dataset. Note that you can
  have paragraphs of text within a single cell in Excel, so please do provide a
  reasonable summary. You will need to use Alt + Enter (or Alt + Shift + Enter on a Mac)
  to insert a carriage return.

|                 |                                   |  |  |
|-----------------|-----------------------------------|--|--|
| Project ID      | 1                                 |  |  |
| Title           | Example data for the safedata system |  |  |
| Description     | This is an example dataset.       |  |  |

## The access block

!!! Warning "Mandatory block"
    **Embargo date** is only required if the **Access status** value is 'Embargo'.
    **Access conditions** is only required if the **Access status** value is
     'Restricted'.

This block provides the access details for the dataset. You can only provide a
single value for each field.

<!-- markdownlint-disable-next-line MD033 -->
* **Access status**<a name="access-status"></a>: You must enter an Access Status
  of  `Open`, `Embargo` or `Restricted`. We prefer as much data as possible to
  be open access: see the [discussion of data availability](../availability.md).
* **Embargo date**: If you choose embargoed access status then you must also
  enter the date when the embargo will end. This must be an Excel date formatted
  value and you cannot embargo data for more than two years. Do not provide
  access conditions: embargoed datasets become freely available when the embargo
  ends.
* **Access conditions**: If you choose restricted access status then you must
  also provide text describing the access conditions. Do not provide an embargo
  date - restricted datasets are permanently restricted.

|                   |             |  |  |
|-------------------|-------------|--|--|
| Access status     | Embargo     |  |  |
| Embargo date      | 03/09/18    |  |  |
| Access conditions |             |  |  |

## The author block

!!! Warning "Mandatory block"
    **Author ORCID** is optional and the remaining fields are mandatory.

These rows provide contact details for the authors of the data. If the datasets
should be credited to more than one author, then provide sets of details in adjacent
columns. If you have an [ORCID](https://orcid.org/), provide it here: this is a
good way to help link all of your academic outputs to you!

Affiliation and email are also optional, but we would **very much prefer
complete author metadata** (name, affiliation, email) for all authors. However,
we realise that sometimes this isn't possible: if you're uploading data
collected by past students who you've lost contact with, then you might not have
these details for **any** author.

Author names **must be** formatted as "last name, first name": "Orme, David" not "David
Orme". Please provide just the numeric part of the ORCIDs, as shown below, not the full
URL `http://orcid.org/0000-0002-7005-1394`.

|                    |                         |  |  |
|--------------------|-------------------------|--|--|
| Author name        | Orme, David             |  |  |
| Author email       | <d.orme@imperial.ac.uk>   |  |  |
| Author affiliation | Imperial College London |  |  |
| Author ORCID       | 0000-0002-7005-1394     |  |  |

!!! Important
    The authors provided here will form part of the permanent citation for the
    published dataset. Authorship on published datasets should be treated in the
    same way as you would consider authorship on papers: you should include not
    only the people responsible for physically collecting the data but also other
    researchers who facilitated the work, such as project supervisors and local
    collaborators.

## The data worksheet block

!!! Warning "Mandatory block"
    The **Worksheet external file** field is only required if a worksheet entry
    describes data held in another file. All other fields are mandatory,

Each **data worksheet** must be described here - do not include the Taxa and
Locations worksheet in this block. As with the authors, you can describe
multiple sheets in adjacent columns.

* The **Worksheet name** row must contain the label of a worksheet in the
   workbook: that is, the **exact text** shown on the worksheet tab at the
   bottom.
* The **Worksheet title** and **Worksheet description** rows are free text to
  provide a longer title and a summary description of the contents of a given
  sheet.
* You only need to use the **Worksheet external file** row if a data worksheet
  describes tabular data held in an external file. The value must then be a
  filename which appears in the [**External file**](#external-files-block)
  block.

<!-- markdownlint-disable MD013 -->
|                         |                        |                     |                                  |
|-------------------------|------------------------|---------------------|----------------------------------|
| Worksheet name          | DF                     | Incidence           | Transects                        |
| Worksheet title         | My shiny dataset       | My incidence matrix | Bait trap transect lines         |
| Worksheet description   | This is a test dataset | A test dataset too  | Attribute table for transect GIS |
| Worksheet external file |                        |                     | BaitTrapTransects.geojson        |
<!-- markdownlint-enable MD013 -->

## Keywords block

!!! Warning "Mandatory block"

This row allows you to enter a set of keywords for the dataset, with one keyword
(or short phrase) **per cell** in the row. Do not use lists of keywords within
cells and provide one set of keywords for the whole dataset, not one per data
worksheet.

|          |           |           |  |
|----------|-----------|-----------|--|
| Keywords | Keyword 1 | Keyword 2 |  |

## External files block

!!! Note "Optional block"
    You only need to provide this information if you are also providing data in
    other file formats. If you do provide this block, all rows are mandatory.

You can include files in other file formats in your data submission as described
[here](other_formats.md). If you do so, then these files must be listed in this
block: we use this information to ensure that all the correct datafiles have
been uploaded to Zenodo and to provide a description in the Zenodo record.

For each file you must provided the exact filename, which **must not** contain
spaces, and a description of the file.

<!-- markdownlint-disable MD013 -->
|                           |                                                         |                                                              |  |
|---------------------------|---------------------------------------------------------|--------------------------------------------------------------|--|
| External file             | bait_trap_images.zip                                      | BaitTrapTransects.geojson                                    |  |
| External file description | Zip file containing 5000 JPEG images of bait trap cards | GeoJSON file containing polylines of the bait trap transects |  |
<!-- markdownlint-enable MD013 -->

## Publication DOI block

!!! Note "Optional block"

This block allows you to provide DOIs for publications using the data here. You
can add multiple DOIs, one per cell in the row. Please format the DOI as a URL
using `https://doi.org/` before the DOI, so
`https://doi.org/10.1098/rstb.2011.0049` not `DOI:10.1098/rstb.2011.0049`. We do
also accept `http://doi.org/`, `http://dx.doi.org/` and `https://dx.doi.org/` as
the root of the URL.

|                 |                                            |  |  |
|-----------------|--------------------------------------------|--|--|
| Publication DOI | `https://doi.org/10.1098/rstb.2011.0049` |  |  |

## Funders block

!!! Note "Optional block"
    Although the funders block is optional, you should provide it in most cases as
    you must provide details of any funding that lead to the collection of the
    data.

    This is particularly important for RCUK funded research, who have agreed to
    let us host all of the SAFE data under a common portal on the condition that
    Research Council funding is clearly acknowledged.

    **Funding body** and **Funding type** are mandatory, but please do provide a
    reference number and a link if possible.

Funding details are provided by completing a block as follows and,  as with
Authors and Worksheets, you can use multiple columns to acknowledge more than
one funder.

|                   |                                                        |  |  |
|-------------------|--------------------------------------------------------|--|--|
| Funding body      | NERC                                                   |  |  |
| Funding type      | Standard grant                                         |  |  |
| Funding reference | NE/K006339/1                                           |  |  |
| Funding link      | `https://gtr.ukri.org/projects?ref=NE%2FK006339%2F1` |  |  |

## Permits block

!!! Note "Optional block"
    If you provide permit details, all the fields are required.

Permits are required for nearly all research conducted at SAFE. Use this block
to record the permits used to collect this data. The permit type value must be
one of `research`, `export` or `ethics`. Again, you can use multiple columns to
record multiple permits.

|                  |                           |  |  |
|------------------|---------------------------|--|--|
| Permit type      | Research                  |  |  |
| Permit authority | Sabah Biodiversity Centre |  |  |
| Permit number    | ABC-123-456               |  |  |

## Gemini metadata

We provide XML metadata for all datasets that is compliant with the UK
[GEMINI](https://www.agi.org.uk/agi-groups/standards-committee/uk-gemini)
standard. This standard includes mandatory **temporal** and **geographic**
extents.

Ordinarily, the dataset checking process will calculate these extents
automatically from the reported locations for the geographic extent and from any
date or datetime fields for the temporal extent. However, if we cannot populate
these extents from the datasets, then you will have to provide additional rows
in your Summary worksheet that provide extent metadata as described below.

!!! Information "Summary extents with location and time data"
    It is _not_ an error to provide extent data in the Summary sheet when there
    is time and location data in the data worksheets. Location and temporal data
    can be incomplete, and so providing wider extents in the Summary is fine.

    If you do this, you _will_ get a warning that both the Summary and data
    extents exist, but this is to check that this is intentional. The Summary
    extents _must_ however be wider than the data extents.

### Temporal extents

!!! Warning "Possibly mandatory block"
    If you do need to provide this block (see [above](#gemini-metadata)) then all
    rows are required.

The start and end date values must be provided as an Excel date formatted cell.

|            |            |  |  |
|------------|------------|--|--|
| Start Date | 01/06/2015 |  |  |
| End Date   | 11/07/2015 |  |  |

### Geographic extents

!!! Warning "Possibly mandatory block"
    If you do need to provide this block (see [above](#gemini-metadata)) then all
    rows are required.

The geographic extents must be provided as decimal degrees (16.75) not degrees,
minutes and seconds (16° 45' 00'") or degrees and decimal minutes (16° 45.00).

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

## Complete example summary table

<!-- markdownlint-disable MD013 -->
|                           |                                                             |                                                              |                                  |
|---------------------------|-------------------------------------------------------------|--------------------------------------------------------------|----------------------------------|
| Project ID                | 1                                                           |                                                              |                                  |
| Title                     | Example data for the safedata system                          |                                                              |                                  |
| Description               | This is an example dataset.                                 |                                                              |                                  |
| Access status             | Embargo                                                     |                                                              |                                  |
| Embargo date              | 03/09/18                                                    |                                                              |                                  |
| Access conditions         |                  |               |
| Author name               | Orme, David                                                 |                                                              |                                  |
| Author email              | <d.orme@imperial.ac.uk>                                       |                                                              |                                  |
| Author affiliation        | Imperial College London                                     |                                                              |                                  |
| Author ORCID              | 0000-0002-7005-1394                                         |                                                              |                                  |
| Worksheet name            | DF                                                          | Incidence                                                    | Transects                        |
| Worksheet title           | My shiny dataset                                            | My incidence matrix                                          | Bait trap transect lines         |
| Worksheet description     | This is a test dataset                                      | A test dataset too                                           | Attribute table for transect GIS |
| Worksheet external file   |                                                             |                                                              | BaitTrapTransects.geojson        |
| Keywords                  | Keyword 1                                                   | Keyword 2                                                    |                                  |
| External file             | bait_trap_images.zip                                          | BaitTrapTransects.geojson                                    |                                  |
| External file description | Zip file containing 5000 JPEG images of bait trap cards     | GeoJSON file containing polylines of the bait trap transects |                                  |
| Publication DOI           | `https://doi.org/10.1098/rstb.2011.0049`                  |                                                              |                                  |
| Funding body              | NERC                                                        |                                                              |                                  |
| Funding type              | Standard grant                                              |                                                              |                                  |
| Funding reference         | NE/K006339/1                                                |                                                              |                                  |
| Funding link              | `https://gtr.ukri.org/projects?ref=NE%2FK006339%2F1`      |                                                              |                                  |
| Permit type               | Research                                                    |                                                              |                                  |
| Permit authority          | Sabah Biodiversity Centre                                   |                                                              |                                  |
| Permit number             | ABC-123-456                                                 |                                                              |                                  |
| Start Date                | 01/06/2015                                                  |                                                              |                                  |
| End Date                  | 11/07/2015                                                  |                                                              |                                  |
| West                      | 116.75                                                      |                                                              |                                  |
| East                      | 117.82                                                      |                                                              |                                  |
| South                     | 4.50                                                        |                                                              |                                  |
| North                     | 5.07                                                        |                                                              |                                  |
<!-- markdownlint-enable MD013 -->
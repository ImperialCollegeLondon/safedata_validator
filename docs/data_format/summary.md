# The Summary worksheet

This worksheet contains a simple set of rows describing the dataset and identifying the
spreadsheets that contain data tables. Each row is labelled on the left in the first column and
then the description data should be typed in the columns to the right.

The following example shows the required rows. You must include **all of these rows** even if
you don't provide that metadata: just leave the row blank. Examples of where this might happen
are embargo date for Open access datasets and Author affiliations, emails and ORCIDs.

|  |  |  |  |
|--------------------|-----------------------------------|---|---|
| SAFE Project ID       | 1 | | |
| Title                 | Example data for the SAFE Project  | | |
| Description           | This is an example dataset.  | | |
| Access status         | Embargo | | |
| Embargo date          | 03/09/18 | | |
| Author name           | Orme, David  | | |
| Author email          | d.orme@imperial.ac.uk  | | |
| Author affiliation    | Imperial College London  | | |
| Author ORCID          | 0000-0002-7005-1394  | |
| Worksheet name        | DF                        | Incidence | Transects |
| Worksheet title       | My shiny dataset          | My incidence matrix  | Bait trap transect lines |
| Worksheet description | This is a test dataset    | A test dataset too | Attribute table for transect GIS |
| Worksheet external file         |                           |           | BaitTrapTransects.geojson |
| External file              | BaitTrapImages.zip | BaitTrapTransects.geojson | |
| External file description       | Zip file containing 5000 JPEG images of bait trap cards | GeoJSON file containing polylines of the bait trap transects | |
| Keywords              | Keyword 1                 | Keyword 2 | |
| Publication DOI       | https://doi.org/10.1098/rstb.2011.0049 |  | |

The first three rows are simple:

 * **SAFE Project ID**: This is the project number from the SAFE project website. When you upload your dataset, you will also be asked to choose a project for your dataset: these two numbers must match. Note that you can only upload data to a project of which you are a member.
* **Title**: This should be a short informative title for the dataset: it will be used as the public title for the dataset so make sure it is clear and grammatical! 
* **Description**: This will be the public description of the dataset. Note that you can have paragraphs of text within a single cell in Excel, so please do provide a reasonable summary. You will need to use Alt + Enter (or Alt + Shift + Enter on a Mac) to insert a carriage return.


# Access status

The **Access status** and **Embargo date** rows allow you to specify when your data will be publically available. We would prefer as much data as possible to be open access but you can embargo your data. You must enter an Access Status of either  `Open` or `Embargo`. If you want to embargo your data, then you must also enter the date when the embargo will end: you cannot embargo data for more than two years. 

# The author block

!!! Important

	The authors provided here will form part of the permanent citation for the published dataset.
	Authorship on published datasets should be treated in the same way as you would consider authorship on
	papers: you should include not only the people responsible for physically collecting the data but also
	other researchers who facilitated the work, such as project supervisors and local collaborators.

These rows provide contact details for the authors of the data. If the datasets should be
credited to more than author, then provide sets of details in adjacent columns. If you have an
[ORCID](https://orcid.org/), provide it here: this is a good way to help link all of your
academic outputs to you!

Affiliation and email are also optional, but we would **very much prefer complete author
metadata** (name, affiliation, email) for all authors. However we realise that sometimes this
isn't possible: if you're uploading data collected by past students who you've lost contact
with, then you might not have these details for **any** author.

Author names **must be** formatted as "last name, first name": "Orme, C David L" not "C David L
Orme".


# The data worksheet block

Each **data worksheet** must be described here - do not include ''Taxa'' and ''Locations''
worksheet in this block. As with the authors, you can describe multiple sheets in adjacent
columns.

 * The ''Worksheet name'' row must contain the label of a worksheet in the workbook: that is, the **exact text** shown on the worksheet tab at the bottom. 
* The ''Worksheet title'' and ''Worksheet description'' rows are free text to provide a longer title and a summary description of the contents of a given sheet. 
* You should only use the ''Worksheet external file'' row if the worksheet describes tabular data held in an external file (see above), which must then appear in ''External file''.

# External files block

As described above, you can include external data files in your data submission. For each file
you must provided the exact filename, which **must not** contain spaces, and a description of
the file. See the [[data_submission_format#non-excel_data_files|section above]] for details.

# Keywords

This row allows you to enter a set of keywords for the dataset, with one keyword (or short phrase) **per cell** in the row. Do not use lists of keywords within cells and provide one set of keywords for the whole dataset, not one per data worksheet.

# Publication DOI

This row allows you to provide DOIs for publication using the data here. You can add multiple DOIs, one
per cell in row. Please format the DOI as a URL using ''https://doi.org/'' before the DOI, so
''https://doi.org/10.1098/rstb.2011.0049'' not ''DOI:10.1098/rstb.2011.0049''. We do also accept
''http://doi.org/'', ''http://dx.doi.org/'' and ''https://dx.doi.org/'' as the root of the URL.

# Funders block

!!! Important

	You must provide details of any funding that lead to the collection of the data. This is particularly
	important for RCUK funded research, who have agreed to let us host all of the SAFE data under a common
	portal on the condition that Research Council funding is clearly acknowledged. 

Funding details are provided by completing a block as follows:

|  |  | 
|---|---|
| Funding body | NERC | 
| Funding type | Standard grant | 
| Funding reference | NE/K006339/1 |
| Funding link | https://gtr.ukri.org/projects?ref=NE%2FK006339%2F1 |

The first two lines are mandatory, but please do provide a reference number and a link if
possible. As with Authors and Worksheets, you can use multiple columns to acknowledge more than
one funder.

# Permits block

Permits are required for nearly all research conducted at SAFE. Use this block to record the permits used to collect this data. The permit type can be one of `research`, `export` and `ethics`. You also have to record the authority that granted the permit and the permit reference number.

|  |  | 
|---|---|
| Permit type | Research | 
| Permit authority | Sabah Biodiversity Centre |
| Permit number | ABC-123-456 |

# Temporal and geographic extents

We provide XML metadata for all datasets that is compliant with the UK
[GEMINI](https://www.agi.org.uk/agi-groups/standards-committee/uk-gemini) standard. This
standard includes mandatory temporal and geographic extents. Our dataset checking process will
ordinarily calculate these extents automatically from the reported locations for the geographic
extent and from any date or datetime fields for the temporal extent.

However, if your dataset does not use locations and/or does not provide date fields within your
data worksheets, then we cannot set these extents. In this case, you will have to provide
additional rows in your Summary worksheet that provide temporal and/or geographic extent
metadata. These extra rows should be labelled in Column A as in the example below and contain a
single value in Column B. Note that geographic extents must be provided as decimal latitude and
longitude values and dates must be Excel formatted dates (not text, for example!).

|  |  | 
|---|---|
| Start Date | 01/06/2015 | 
| End Date | 11/07/2015 | 
| West | 116.75 |
| East | 117.82 |
| South | 4.50 | 
| North | 5.07 |

The geographic bounds in the example cover Maliau, Danum and the SAFE Project experimental site
and surrounding area. While we would prefer something a bit more precise, these are sufficient
to provide a geographic extent for most work at SAFE.
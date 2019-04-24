# Data formatting

Part of the agreement for research projects working at the SAFE Project is that all project
data are submitted to the SAFE Project data repository, so that it is available to future
researchers. In order to make it easy for data to be found and used in the future, we need
researchers to provide some (relatively!) simple information in their datafiles.

The process in overview is:

* Prepare your data following the formatting guide below.
* Go to the SAFE Project website and [Submit your dataset](https://www.safeproject.net/datasets/submit_dataset).
* We will automatically validate the data formatting. 
* If the validation succeeds then we will publish it on Zenodo.
* If the validation fails then you will get an error report so you can fix the problems and resubmit.

Note that you **must not** publish your dataset directly to Zenodo as this skips the validation
step and means that the datasets are not linked together under a common account.

## Authorship and Funding

All datasets submitted to SAFE are published to the [Zenodo data
repository](https://zenodo.org/communities/safe). Zenodo issues DOIs for all data depositions,
making your datasets easily citable.

You will need to provide a list of authors for submitted datasets, which will form part of the
permanent citation for the published dataset. Authorship on published datasets should be
treated in the same way as you would consider authorship on papers: you should include not only
the people responsible for physically collecting the data but also other researchers who
facilitated the work, such as project supervisors and local collaborators. Similarly, you must
acknowledge any funding that supported your research.

## Format overview

In most cases, we expect the data format to be provided as **tabular datasets stored in Excel
spreadsheets**. This accounts for the vast majority of the data files used by researchers. If
you have data that is stored in other sorts of files, because of formatting or file size
issues, then you will still need to fill in an Excel template containing the key metadata for
your dataset. We will then upload your data files manually to Zenodo and link the metadata to
them.

!!! Note
    
    The details described below will be used to automatically publish your data to Zenodo. You
    should choose titles, descriptions and keywords that you would be happy to
    be**permanently** associated with your dataset!

The basic format for a SAFE dataset submission is an Excel Workbook, which must contain the
following three worksheets:

* **Summary**: This contains some simple information about the authors of the dataset, access rights and the individual data tables in the dataset. 
* **Taxa**: This describes all the taxa used in the dataset. 
* **Locations**: This describes all the sampling locations used in the dataset.

After these worksheets come your data tables. You should label these sheets with a sensible
name (not 'Sheet1'!) and each data table must be described in the Summary worksheet. You can
include as many data tables as you like in a single dataset: we don't want you to spend time
rearranging your data and are happy just to take the data in the natural tables you already use.

### File naming

Use a simple short name for your spreadsheet - there will be a lot of information giving more
detail inside. Please **do not use spaces** in your file name - you can use underscores to
separate words.

### Spreadsheet Template and Examples

Click on this link to download the {{ :working_at_safe:template.xlsx |spreadsheet template }}
containing the required worksheets, labels and headers.

You can also look at existing published datasets to see how the format is used:

* [https://safeproject.net/datasets/view_datasets](https://safeproject.net/datasets/view_datasets)
* [https://zenodo.org/communities/safe/](https://zenodo.org/communities/safe/)

### Format checking

We've tried to make the description below as clear as possible but in order to help you prepare
your file:

* It is easier to follow an example than to follow a description, so please use the template and look at the examples. 
* We use a Python program to automatically check the formatting of datasets. When you submit a file, you will get a report back from this program that will highlight any problems with your dataset. If there are problems, fix them and replace the submitted file. Once the file passes through the checker without problem, we will double check the file and then publish your dataset.

### Checking your own data

If you want to check your formatting yourself before submitting it then the code used to check
Excel datasets is freely available online [here](https://github.com/ImperialCollegeLondon/safe_dataset_checker). The link also provides
instructions on how to use the code to check your data. You will need a computer with Python
installed and which is connected to the internet (although the program can be setup to allow
offline use).

### Data availability

When your dataset is published, the metadata will be immediately publicly visible. This
includes details of the data fields, the spatial scope, the date range and the like. If you set
the access status as ''Open'', then the Excel file itself will also be immediately publicly
available.

We would prefer that as much data as possible is submitted with ''Open'' access status, but if
you want to restrict access to the data while you work on papers, then you can use the
''Embargo'' access status and set an embargo date. The metadata will still be visible, so that
researchers can see that the data exists but the data itself will only become available once
the embargo date has passed. You cannot embargo a dataset for more than two years.

Obviously, you can choose to provide embargoed data to other researchers within the embargo
period. If researchers contact the SAFE Project for access to data during the embargo period,
we will **always** pass the request on to you.

## FAQ or 'What does this error mean?'

There is a separate page [[working_at_safe:data_format_faq|here]] for frequently asked
questions about issues with data formatting.

## Non-Excel data files 

If you are providing data in non-Excel formats (e.g. zip file of images, SQL dumpfile, audio
files) then use the Excel template to provide the basic summary metadata. You will need to
provide us with a set of data files, which we will upload to Zenodo. In your summary template
file, use the ''External file'' and ''External file description'' rows to provide a name and
description for **every** external file you want to upload. If you have have many small files,
then zipping them up and providing a broad description is fine.

Some key points:

 * File names **must not** contain spaces and you should use the **exact** file name of the data file in the metadata description. 
* You **must** also complete the Taxa and Locations worksheets if applicable. We will not be able to validate them against the external data but we will be able to add the taxa and locations to our data index. 
* You can mix standard Excel data worksheets and references to external files, as in the example below. In this way you can provide bulk raw data and processed data tables. 
* Please do not use this approach to submit spreadsheet data in a different format (e.g. OpenOffice).

### Tabular data in external files 

Some possible external data files contain tabular data: for example, a zipped file containing a
large CSV dataset, a SQL data dump or a shapefile with an attribute table. If so, then please
do fill in a data worksheet describing the table metadata as described below. The data
worksheet does not need to contain any of the actual data but it will allow us to index the
available fields.

 * Create a worksheet for the table and describe the fields as described below. 
 * In the worksheet summary fields, use the ''Worksheet external file'' row to give the name of the file containing the tabular data. This must appear in the list of external file names and description. 
 * If an external file contains more than one data table, then it is fine to refer to the same file in more than one data worksheet.

## The Summary worksheet

This worksheet contains a simple set of rows describing the dataset and identifying the
spreadsheets that contain data tables. Each row is labelled on the left in the first column and
then the description data should be typed in the columns to the right.

The following example shows the required rows. You must include **all of these rows** even if
you don't provide that metadata: just leave the row blank. Examples of where this might happen
are embargo date for Open access datasets and Author affiliations, emails and ORCIDs.

|  |  |  |  |
|--------------------|-----------------------------------|---|---|
| SAFE Project ID       | 1 | | |
| Access status         | Embargo | | |
| Embargo date          | 03/09/18 | | |
| Title                 | Example data for the SAFE Project  | | |
| Description           | This is an example dataset.  | | |
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

The first rows are simple:

 * **SAFE Project ID**: This is the project number from the SAFE project website. When you upload your dataset, you will also be asked to choose a project for your dataset: these two numbers must match. Note that you can only upload data to a project of which you are a member.
* **Access status** and **Embargo date**: As described above, the access status of the datasets can either be ''Open'' or ''Embargo''. If you want to embargo your data, then provide a date when the embargo will end: you cannot embargo data for more than two years. 
* **Title**: This should be a short informative title for the dataset: it will be used as the public title for the dataset so make sure it is clear and grammatical! 
* **Description**: This will be the public description of the dataset. Note that you can have paragraphs of text within a single cell in Excel, so please do provide a reasonable summary. You will need to use Alt + Enter (or Alt + Shift + Enter on a Mac) to insert a carriage return.

### The author block

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

### The worksheet block

Each **data worksheet** must be described here - do not include ''Taxa'' and ''Locations''
worksheet in this block. As with the authors, you can describe multiple sheets in adjacent
columns.

 * The ''Worksheet name'' row must contain the label of a worksheet in the workbook: that is, the **exact text** shown on the worksheet tab at the bottom. 
* The ''Worksheet title'' and ''Worksheet description'' rows are free text to provide a longer title and a summary description of the contents of a given sheet. 
* You should only use the ''Worksheet external file'' row if the worksheet describes tabular data held in an external file (see above), which must then appear in ''External file''.

### External files

As described above, you can include external data files in your data submission. For each file
you must provided the exact filename, which **must not** contain spaces, and a description of
the file. See the [[data_submission_format#non-excel_data_files|section above]] for details.

### Keywords

Provide keywords for the dataset here, with one keyword (or short phrase) **per cell** in the
row. Do not use lists of keywords within cells and provide one set of keywords for the whole
dataset, not one per data worksheet.

### Publication DOI

Provide DOIs for publication using the data here and you can add multiple DOIs, one per cell in
the row. Please format the DOI as a URL using ''https://doi.org/'' before the DOI, so
''https://doi.org/10.1098/rstb.2011.0049'' not ''DOI:10.1098/rstb.2011.0049''. We do also
accept ''http://doi.org/'', ''http://dx.doi.org/'' and ''https://dx.doi.org/'' as the root of
the URL.

### Funders

You must provide details of any funding that lead to the collection of the data. This is
particularly important for RCUK funded research, who have agreed to let us host all of the SAFE
data under a common portal on the condition that Research Council funding is clearly
acknowledged. To do so, please provide a block as follows:

|  |  | 
|---|---|
| Funding body | NERC | 
| Funding type | Standard grant | 
| Funding reference | NE/K006339/1 |
| Funding link | https://gtr.ukri.org/projects?ref=NE%2FK006339%2F1 |

The first two lines are mandatory, but please do provide a reference number and a link if
possible. As with Authors and Worksheets, you can use multiple columns to acknowledge more than
one funder.

### Temporal and geographic extents

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

## The Taxa worksheet

Many datasets will involve data taken from organisms, whether that is a count of the number of
individuals or measurement of a trait such as body length. In order to help us keep track of
taxa, all datasets using taxa **must** contain a ''Taxa'' spreadsheet, providing taxonomic
information.

Note that you must only provide details for taxa actually used in the data worksheets. This
ensures that the taxonomic index for a dataset is accurate and also double checks that it the
omission of a taxon from the data worksheets is not an error.

### Taxon validation

In order to help keep the taxonomy as clean as possible and to allow us to index the taxonomic
coverage of datasets, we will check all taxon names in Taxa worksheet against the GBIF backbone
taxonomy database. If you want to check your taxon names and ranks, then the search engine is
here:

[https://www.gbif.org/species/search](https://www.gbif.org/species/search)

No online taxonomy is ever going to be 100% up to date (or 100% agree with your taxonomic
usage!) but the GBIF backbone has very good taxonomic coverage and is well curated.

### Taxon table layout

The table format looks like this:

| Name | Taxon name | Taxon type | Taxon ID | Parent name | Parent type | Parent ID |
|---|---|---|---|---|---|---|
| Crematogaster borneensis | Crematogaster borneensis | Species |   |   |   |   |
| Dolichoderus sp. | Dolichoderus | Genus |   |   |   |   |
| Gannets	   | Morus | Genus | 2480962 |   |   |   |
| Morphospecies 1 | NA | Morphospecies |   | Formicidae | Family |   |

The table must contain column headers in the **first row** of the worksheet. The headers must
include:

 * **''Name''**: This column must contain a local name for **all** of the taxa that you are going to use in the rest of the dataset. You cannot have duplicated names! Note that these can be abbreviations or codes: if you want to use ''Crbe'' in your data worksheets, rather than typing out ''Crematogaster borneensis'' every time, then that is fine.

!!! Note

    These are the names that you are going to use in your data worksheet. The 
	other columns are to help us validate the taxonomy of your names.

 * **''Taxon name''**: This column must contain the scientific name of the taxon, which will be used for taxon validation via GBIF. Note that this should not include any taxon authority, so _Panthera tigris_ not _Panthera tigris_ (Linnaeus, 1758).

 * **''Taxon type''**: This column must provide the **taxonomic type** of the named taxon, which is usually the taxonomic **rank**. For example, the taxon _Pongo pygmaeus_ would be of type ''Species'' and the taxon Formicidae would be of type ''Family''. There are two additional possible values which are described in more detail below: ''Morphospecies'' and ''Functional group''.

 * **''Taxon ID''**: This column is optional - you don't have to include it and you don't have to fill in every row if it is included. It is **only needed** if the taxon name and rank are ambiguous.  
     For example, the genus _Morus_ can refer to either mulberries or gannets. In these rare cases, you will need to look up the taxon in GBIF and provide the GBIF ID number. The example in the table allows us to confirm that you mean this _Morus_: [https://www.gbif.org/species/2480962](https://www.gbif.org/species/2480962).

There are three types of taxa, described below, that can't be validated directly against GBIF.
In these cases, the next three columns (**''Parent name''**, **''Parent type''** and **''Parent
ID''**) are used to provide a parent taxon that can be validated and provides a taxonomic hook
to allow us to place the taxon in the backbone taxonomy. They are used in exactly the same way
as the ''Taxon'' columns: the only restriction is that the ''Parent type'' must be one of seven
core ranks used in the GBIF backbone. Again **''Parent ID''** is only needed to resolve ambiguous 
taxon names.

The ''Parent'' columns only need to be filled in for taxa falling in one of the following
groups.

### New and unrecognized taxa

If a taxon is new or not recognized by GBIF (and you're sure you're right!) then provide a
parent name and type to allow us to hook the taxon into the index. For example //Pongo
tapanuliensis// is not currently recognised as a species, so providing //Pongo// as a parent
name of type 'genus' allows us to place the new taxon.

### Morphospecies and Functional groups

For morphospecies and functional groups, the taxon name is the label to be used in the dataset.
Set the Scientific name to be 'NA' - it cannot be blank - and then specify the taxon type as
'Morphospecies' or 'Functional group'.

Now you need to provide a parent taxon and type. The level of taxonomic certainty for
morphospecies and functional groups is quite variable, but we'd like the finest taxonomic level
you can provide. As an example, in the table above, 'Morphospecies #1' is simply identified as
being an ant.

### Less common taxonomic levels

The GBIF backbone taxonomy only includes the following eight major levels: Kingdom, Phylum,
Order, Class, Family, Genus, Species and Subspecies. If you need to use taxa defined at any
intermediate levels, then again provide a parent taxon and type. For example, if you were
counting bees and only identifying to tribe level (Bombini, Euglossini, etc.) then the parent
family Apidae would allow us to hook the taxa into the backbone taxonomy. The subfamily Apinae
would be more precise, but subfamily isn't one of the backbone taxonomic levels.

### My data doesn't contain taxa

Fine. You can omit the Taxa worksheet!

## The Location worksheet

Like the Taxa worksheet, all locations in your data worksheets need to be listed in this
worksheet. By **location**, we mean the common frequently used areas in which research has
happened at SAFE. You might have more detail about the precise place you worked in your dataset
- great! - but using these known locations allows us to get broad spatial data on sampling
relatively simply.

So, we expect you'll have a relatively small set of location names in your data sheets, all of
which should appear in this worksheet: the worksheet should contain a column of location names,
with the column header 'Location name' in the first row.

### Location verification

The location names are checked against the location names known in the SAFE gazetteer. You can
look at the gazetteer webpage to see the available sites and to download location data:

[https://www.safeproject.net/info/gazetteer](https://www.safeproject.net/info/gazetteer)

If you want to get a list of valid location names for use in a program or script, then we
provide a web service that returns a list of valid names as a JSON object:

[https://www.safeproject.net/call/json/get_locations](https://www.safeproject.net/call/json/get_locations)

For example, in R:

```r
> library(jsonlite) 
> locations <- fromJSON("https://www.safeproject.net/call/json/get_locations") 
> str(locations) List of 1 $
locations: chr [1:2691] "SAFE_camp" "Flux_tower" "A_1" "A_2" ... 
```

### New locations

If your data comes from genuinely new locations or uses a sampling structure (e.g. a grid or
transect) that is likely to be used again in the future, then you can create new location names
and include them in your locations table. We will then consider adding them to the Gazetteer.

If you include new locations then you will need to include the following columns in your
''Locations'' worksheet:

 * ''**New**'': This should simply contain ''Yes'' or ''No'' to show which rows contain new
locations. You cannot create a new location with a name that matches an existing location in
the Gazetteer. 
* ''**Latitude**'' and ''**Longitude**'': these should provide GPS coordinates
for the new site. These must be provided as decimal degrees (not degrees minutes and seconds)
and please provide 6 decimal places in your coordinates. This level of precision is around ten
centimetres and, although the GPS from the field is highly unlikely to be accurate to this
level, we want to record as much sampling precision as possible. If you don't have any GPS data
for the new location, please explicitly enter ''NA'' in these fields. 
* ''**Type**'': For most
new locations, this will be ''POINT'', so the latitude and longitude are sufficient. New linear
sampling features (e.g. transects) are ''LINESTRING'' and sampling areas are ''POLYGON''. In
these cases, you will need to email [[mailto:admin@safeproject.net|the SAFE administrators]]
and provide a GIS file containing the spatial information for your new locations.

You only need to provide ''Latitude'', ''Longitude'' and ''Type'' in the rows for new
locations: these rows can be blank for locations that are already in the gazetteer.

### My data doesn't include any locations

You don't **have** to include the Locations worksheet, although it would be very unusual to omit it.
Possible examples:

* You are working with lab data (and don't need to say where specimens came from in the field)
* You are collecting data haphazardly from across the landscape, for example tracking animal movements, and the data isn't tied to particular sampling locations. We would then want GPS data for each observation!

## Data worksheets

Finally, we get to the worksheets containing your actual data!

### Field metadata

The top rows of the worksheet are used to provide metadata descriptors for each of the columns
('fields') in your data worksheet. Each descriptor row has a label, which must appear in Column
A of the worksheet, with the value for each field appearing above that column.

The following are the **mandatory field descriptors**, which are needed for all fields and
which **cannot be blank**.

 * **''field_type''**: This has to be one of the following values indicating the field type
(see the options below). * **''description''**: a short description of the field *
**''field_name''**: the name of the variable. The name format should be suitable for loading
into an analysis package and should not contain spaces: use an underscore (''_'') to put gaps
in names. It definitely must not have white space at the start or end and it should only use
standard ASCII characters. This descriptor **must always be the last descriptor row**,
immediately above the data, so that it can be used as field headers when loading data from the
file for analysis.

There are also some **additional field descriptors**, which are mandatory for some data types
(see the descriptions of the data types below). The options are:

 * **''levels''**: contains the set of level names used in a categorical variable. *
**''method''**: a contain a short description of the method and equipment used to record
numeric, abundance and trait data. * **''units''**: the units of numeric or trait variables. *
**''taxon_name''**: the name of the taxon for which all of the trait or abundance data in a
field is recorded. The taxon name must appear in the Taxa worksheet. * **''taxon_field''**: the
name of a taxon field in the datasheet which shows the taxon for which trait or abundance data
on that row is recorded. * **''interaction_name''**: a set of names giving the interacting taxa
for interaction data . * **''interaction_field''**: a set of field names, where the rows of the
field give the interacting taxa for interaction data.

These descriptors only have to be completed for the appropriate data types: **leave them
blank** for any fields that don't require them. ==== Missing data ====

If your data worksheets contain missing data, **you must enter 'NA' in those cells**, not just
leave them blank. This is to make it absolutely unambiguous that a given value is actually
missing. We know this is picky but it can be absolutely vital: for example, does a blank cell
in an abundance matrix mean that the species wasn't seen (so the cell should be zero) or that
the trap for that species fell over and you don't know if it was recorded (so it should be NA).

==== Row numbers ====

You **must** number the rows in your data worksheet. The row numbers must start at 1 in the
cell directly under the ''field_name'' descriptor, increase by 1 as you move down through the
cells and must continue down to the last row containing data. The row numbers must not extend
below the data: the template numbers rows down to 1000, so delete the numbers for any unused
rows in your data!

===== Field types =====

This section shows the options that can appear in the ''field_type'' descriptor, along with any
further descriptors that might be needed. See the sections below for details on formatting, but
the available types are:

 * **''Date''**, **''Datetime''** and **''Time''**: when were the data collected? *
**''Location''**: where was the data collected? * **''Latitude''**, **''Longitude''**: GPS data
for the exact location. * **''Replicate''**: a record of replication, usually just a repeating
set of numbers. * **''ID''**: a column showing any kind of identification code. *
**''Categorical''**: otherwise known as a factor: a variable that puts data into a fixed set of
groups. * **''Ordered Categorical''**: a factor where there is a logical order to the levels. *
**''Numeric''**: all kinds of numeric data. * **''Taxa''**: what taxa was the data collected
from? * **''Abundance''**: for abundance/density/presence data collected about a taxon. *
**''Categorical Trait''**: for categorical data collected on a taxon. * **''Numeric Trait''**:
for numeric data collected on a taxon. * **''Categorical Interaction''**: for categorical data
on interactions between taxa. * **''Numeric Interaction''**: for numeric data on interactions
between taxa.

==== Date, Datetime and Time ====

We have three kinds of date and time fields!

 * **''Datetime''**: The data in the field includes **both a time and a date** (e.g. 21/05/2016
15:32), which could be a visit time and day to a site or when a camera trap was deployed or
similar data. * **''Date''**: The data in the field **only specifies a date** (e.g.
21/05/2016). * **''Time''**: The data in the field **only specifies a time** (e.g. 15:32).

We don't mind how you provide the date and time information but you do need to be consistent
within a field.

<WRAP tip>Excel cell formatting can make this confusing. Both date and time are stored in Excel
as a single number (N: days since the beginning of January 1900). If N < 1 it represents a time
and if N > 1 it is a date. However, cell formatting can mislead you as to what is actually
stored in the cell.

 * 0.75 is the time ''18:00'', but it could also display in Excel as ''00/01/1900'' or
''00/01/1900 18:00'' if formatted as a date. This is reasonably easy to spot because of the 0th
of January! * 12 is the date ''12/01/1900'' but it could display as the time ''00:00'' or as
''12/01/1900 00:00''. * 12.75 is the datetime ''12/01/1900 18:00'' but could display as the
time ''18:00'' or as ''12/01/1900''.

Note that the value 12 is ambiguous, because Excel doesn't differentiate between integer and
float numbers: it could just refer to the day (the integer 12) or mean exactly midnight on the
day (the float 12.0). This is one reason why we have the three data types!</WRAP>

==== Locations ====

Columns of this type contain location labels showing where the data in the row was recorded.
All of the labels must have been included in the Locations worksheet.

==== Taxa ====

Columns of this type contain taxon names showing the taxon from which other data in the row was
recorded. All of the values in the row must appear in the ''Taxon Names'' column in the Taxa
worksheet.

==== Replicate and ID ====

Both ''Replicate'' and ''ID'' fields could contain almost any values. Replicates are typically
just shown with repeating numbers, but researchers could use other formats. ID can represent
lots of things (for example, PIT tag numbers for individual organism, fine scale spatial
sampling ID, batch number for reagents) and again could have almost any format.

So, both ''ID'' and ''Replicate'' fields are checked for missing data (NAs are permitted) but
no other validation occurs.

==== Categorical data ====

<wrap INFO>Field descriptor ''levels'' required</wrap>

Both categorical and ordered categorical data (also known as a factors) are made up of a set of
**levels** showing the different groups or treatment. The data in the column then shows which
level applies to each row.

In the ''levels'' descriptor, you must provide a complete set of all the levels used in the
column, which will be checked against the data. The level names must be short text labels. **Do
not use integer level names**: they are harder to interpret in statistical analyses and there
is a real risk that they are analysed as a number by mistake.

The format is that the level names are separated using semi-colons ('';''). For example:

 Control;Logged;Burned

We automatically remove spaces around level labels, so these also work and might be easier to
read:

 Control; Logged; Burned Control ; Logged ; Burned

If the levels aren't obvious, we'd also like label descriptions: they come after each label,
separated by colons ('':''). For example:

 Control:sites in reserve forest;Logged:sites in logged forest;Burned:sites in burned forest

**Do not** use colons or semi-colons in your level names or descriptions!

For ''Ordered Categorical'' fields, the order of the entries in the ''levels'' descriptor
should be the logical order of the factor. For example, an ordered disturbance gradient could
be:

 Primary:primary rainforest;Once:once logged rainforest;Twice:twice logged;Salvage:salvage
logged;Oil palm:plantation

==== Numeric data ====

<wrap INFO>Field descriptors ''method'' and ''units'' required</wrap>

This field type should be used to record numeric variables **except numeric variables recorded
from taxa** (see Traits below). The ''method'' descriptor should include information about how
the variable is measured and the ''units'' descriptor must provide the units used.

Not all numeric variables have methods or units: a column of replicate numbers, for example. If
this is the case, enter ''None'' rather than leaving the descriptors blank. (If you prefer to
use ''Dimensionless'' as the unit for dimensionless quantities then that is also fine!)

==== Abundance and trait data ====

Both traits and abundance data tie a value (category or number) to a single taxon. You need to
format your data so that it is clear which taxon each value comes from. There are two possible
formats:

 - All observations in a column are from **a single taxon**: in this case, you can put a valid
taxon name (see Taxa worksheet) in the ''taxon_name'' descriptor for this column.

++++ Example: Observation counts in separate columns for each taxon | | field_type | Abundance
| Abundance | | taxon_name | Tiger leech | Brown leech | | method | Exhaustive search of 50cm
quadrat | Exhaustive search of 50cm quadrat | | description | quadrat count | quadrat count | |
field_name | tiger_count | brown_count | | 1 | 24 | 12 | | 2 | 62 | 3 | ++++

 - Different rows in the column refer to **different taxa**: in this case, you must also have a
Taxa column and the ''taxon_field'' descriptor needs to contain the field name of the
appropriate Taxa column.

++++ Example: Observation counts with different taxa in rows | | field_type | taxa | Abundance
| | taxon_field | | common_name | | method | | Exhaustive search of 50cm quadrat | |
description | Species found | Number found | | field_name | common_name | leech_count | | 1 |
Tiger leech | 24 | | 2 | Brown leech | 12 | | 3 | Tiger leech | 62 | | 4 | Brown leech | 3 |
++++

It is an error to provide both ''taxon_name'' and ''taxon_field'' descriptors for an Abundance
or Trait field. === Abundance ===

<wrap INFO>Field descriptors ''method'' and one of ''taxon_name'' or ''taxon_field''
required</wrap>

Abundance is used here as an umbrella term to cover a wide range of possibilities from casual
observation data ('We saw two clouded leopards on Friday on the road near F100'), through
presence/absence data to precise measurements of abundances or encounter rate.

The ''method'' descriptor needs to provide a detailed description of the sampling method,
including the area surveyed, the length of time spent sampling, the number of samplers and any
equipment. This should be detailed enough to allow the sampling protocol to be replicated. If
other columns provide sampling information, such as survey time or area, then make this clear.

=== Categorical trait ===

<wrap INFO>Field descriptors ''levels'' and one of ''taxon_name'' or ''taxon_field''
required</wrap>

This is just a categorical variable where the groups apply to a taxa. So, we need information
on the levels used, as for a standard categorical variable, and a link to taxonomic information
as described in the examples above.

=== Numeric trait ===

<wrap INFO>Field descriptors ''units'', ''method'' and one of ''taxon_name'' or ''taxon_field''
required</wrap>

This is just a numeric variable where the groups apply to a taxa. So, we need the method and
units for the values, as for a standard numeric variable, and a link to taxonomic information
as described in the examples above.

==== Interaction data ====

Interaction data is essentially just a column of categorical or numeric data that you want to
associate with (at least) two taxon identities, but there are lots of ways that the taxon
identities could be provided.

Interaction data fields do this by using two alternative descriptors to tie the data to taxa:
''interaction_name'' and ''interaction_field''. Each descriptor can provide one or more taxon
names and optionally their role - the formatting is identical to the categorical data
''levels'' descriptors. So, for example an ''interaction_name'' descriptor might be:

 Moon rat:prey;Clouded leopard:predator;

You can use one or both of the descriptors, depending on how your data is laid out. For the
most common case of two interacting taxa, the following three possibilities exist.

1. Both interacting taxa vary from row to row, so taxon names are provided in two fields

++++ Example: Interacting taxa identified in separate columns | | field_type | Taxon | Taxon |
Categorical Interaction | | interaction_field | | | predator;prey | | levels | | |
success;failure | | method | Visual observation | Visual observation | Visual observation | |
description | Predator observed | Prey observed | Outcome of predation event | | field_name |
predator | prey | outcome | | 1 | Clouded leopard | Brown rat | success | | 2 | Flat headed cat
| Moon rat | failure | ++++

2. Alternatively, all of the data might refer to the same two taxa, so the taxon names can be
provided directly.

++++ Example: Interacting taxa identifed in separate columns | | field_type | Categorical
Interaction | | interaction_name | Clouded leopard:predator;Brown rat:prey | | levels |
success;failure | | method | Visual observation | | description | Outcome of predation event |
| field_name | outcome | | 1 | success | | 2 | failure | ++++

3. Finally, one side of the interaction might vary from row to row but the other side is
constant for all rows.

++++ Example: Interacting taxa identified by name and by column | | field_type | Taxon |
Categorical Interaction | | interaction_name | | Clouded leopard:predator; | |
interaction_field | | prey:prey species; | | levels | | success;failure | | method | Visual
observation | Visual observation | | description | Prey observed | Outcome of predation event |
| field_name | prey | outcome | | 1 | Brown rat | success | | 2 | Moon rat | failure | ++++

You must provide at least two taxon names or fields, but you can provide more if you have
tritrophic interactions! Again, you can use any combination of interaction names and fields to
provide your taxon identities.

=== Categorical interactions ===

<wrap INFO>Field descriptors ''levels'' and ''interaction_name'' and/or ''interaction_field''
required</wrap>

=== Numeric interactions ===

<wrap INFO>Field descriptors ''units'', ''method'' and ''interaction_name'' and/or
''interaction_field'' required</wrap>

==== Comments ====

If you have a free text field with notes or comments, then this is the field type to use. We
don't really check anything in comments fields: they're not expected to be complete data and
you can put anything in them.

A word of caution though: it is //highly unlikely// that anyone will ever read your comments
column again. If there is genuinely important information that might apply across multiple
rows, consider coding it as an explicit variable rather than consigning it to a comments field.

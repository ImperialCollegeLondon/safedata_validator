# Data worksheets

<!-- markdownlint-disable MD033 -->
<style>

/*removing table headers and fixing cell widths so everything lines up.*/
thead {
  display: none;
}

table {
  table-layout: fixed;
}

tbody td {
  width: 14em;
  min-width: 14em;
  max-width: 14em;
  border: 1px solid lightgrey;
}

tbody td:first-child {
  width: 12em;
  min-width: 12em;
  max-width: 12em;
  font-weight: bold;
}

</style>
<!-- markdownlint-enable MD033 -->

Finally, we get to the worksheets containing your actual data! Data worksheets are
structured as a "data frame" with each column containing a different field or variable.
The worksheet is organised into two sections:

* field metadata appears at the top of the sheet, and then
* the actual data for each field appear below the metadata.

!!! Note "Naming"
    The Data worksheets can be given whatever names you feel are most appropriate, with
    a few restrictions.

    Firstly, the following names are not allowed as they are reserved for another
    purpose: "Summary", "Taxa", "GBIFTaxa", "NCBITaxa", "Locations"

    Secondly, all worksheet names should be valid `R` names. See
    [here](https://adv-r.hadley.nz/names-values.html#non-syntactic) for details of what
    constitutes a valid name in `R`.

    Finally, the worksheet names must match with the worksheet names that were provided
    in the [Summary sheet](./summary.md#the-data-worksheet-block).

## Field metadata

Each data worksheet will require a particular set of **metadata descriptors** - some are
mandatory and others are only required for some kinds of field. The field metadata
descriptors are arranged in rows, with the label of the descriptor  for each row
appearing in the first column (Column A).

The following are the **mandatory field descriptors**, which are needed for all
fields and which **cannot be blank**.

* **field_type**: This has to be one of the following values indicating the
  field type (see the options [below](#field-types)).
* **description**: a short description of the field .
* **field_name**: the name of the variable. This descriptor **must always be the last
  descriptor row**, immediately above the data, so that it can be used as field headers
  when loading data from the file for analysis. Field names _have_ to follow the
  convention for field names in R (see
  [here](https://adv-r.hadley.nz/names-values.html#non-syntactic)). Broadly, this should
  use standard [ASCII](https://en.wikipedia.org/wiki/ASCII) characters and should not
  contain spaces or punctuation, although underscores and full stops are permitted. In
  addition, names starting with a number (e.g. `1` or `2day`) are also forbidden, along
  with a small set of [reserved
  words](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Reserved.html).

There are also some **additional field descriptors**, which are mandatory for
some data types (see the descriptions of the data types below). The options are:

* **levels**: contains the set of level names used in a categorical variable.
* **method**: contains a short description of the method and equipment used to
  record numeric, abundance and trait data.
* **units**: the units of numeric, abundance or trait variables.
* **taxon_name**: the name of the taxon for which all of the trait or abundance
  data in a field is recorded. The taxon name must appear in a Taxa worksheet.
* **taxon_field**: the name of a taxon field in the datasheet which shows the
  taxon for which trait or abundance data on that row is recorded.
* **interaction_name**: a set of names giving the interacting taxa for
  interaction data .
* **interaction_field**: a set of field names, where the rows of the field give
  the interacting taxa for interaction data.
* **file_container**: a reference to an external data file name that is
  associated with a field.

These descriptors only have to be completed for the appropriate data types:
**leave them blank** for any fields that don't require them. Unlike missing data (see
below), you do not need to use NA in field metadata.

## Missing data

If your data worksheets contain missing data, **you must enter 'NA' in those cells**,
not just leave them blank. This is to make it absolutely unambiguous that a given value
is actually missing. We know this is picky but it can be absolutely vital: for example,
does a blank cell in an abundance matrix mean that the species wasn't seen (so the cell
should be zero) or that the trap for that species fell over and you don't know if it was
recorded (so it should be NA).

## Formula

You must not submit datasets with formula calculations - obviously, feel free to use
calculations when working with the data but please copy and paste 'As Values' before
submitting the dataset.

## Row numbers

You **must** number the rows in your data worksheet. The row numbers must start at 1 in
the cell directly under the `field_name` descriptor, increase by 1 as you move down
through the cells and must continue down to the last row containing data. The row
numbers must not extend below the data: the template numbers rows down to 1000, so
delete the numbers for any unused rows in your data!

## Field types

This section shows the options that can appear in the `field_type` descriptor, along
with any further descriptors that might be needed. See the sections below for details on
formatting, but the available types are:

* **Date**, **Datetime** and **Time**: when were the data collected?
* **Location**: where was the data collected?
* **Latitude**, **Longitude**: GPS data for the exact location.
* **Replicate**: a record of replication, usually just a repeating set of numbers.
* **ID**: a column showing any kind of identification code.
* **Categorical**: otherwise known as a factor: a variable that puts data into
  a fixed set of groups.
* **Ordered Categorical**: a factor where there is a logical order to the levels.
* **Numeric**: all kinds of numeric data.
* **Taxa**: what taxa was the data collected from?
* **Abundance**: for abundance/density/presence data collected about a taxon.
* **Categorical Trait**: for categorical data collected on a taxon.
* **Numeric Trait**: for numeric data collected on a taxon.
* **Categorical Interaction**: for categorical data on interactions between taxa.
* **Numeric Interaction**: for numeric data on interactions between taxa.
* **File**: for references to external files in data tables.
* **Comments**: for general comments about data.

### Date, Datetime and Time

We have three kinds of date and time fields!

* **Datetime**: The data in the field includes **both a time and a date** (e.g.
  21/05/2016 15:32), which could be a visit time and day to a site or when a
  camera trap was deployed or similar data.
* **Date**: The data in the field **only specifies a date** (e.g. 21/05/2016).
* **Time**: The data in the field **only specifies a time** (e.g. 15:32).

We don't mind if dates and times are provided in separate fields or combined: you just
need to be consistent within a field. You can provide date/time information either as
text (which is quite common from data recorders and the like) or as Excel date/time
formats but again you can only use one of these formats within a field.

1. **Text**. If the field contains text values then the checker will try to interpret
  these as date/time values. We only support [ISO8601
  format](https://en.wikipedia.org/wiki/ISO_8601), so the following are valid examples:

    * `2019-04-24 07:00:12`: a datetime value
    * `2019-04-24`: a date value
    * `07:00:12`: a time value

2. **Excel date and time format**. If the field contains Excel date/time cells then we
   will validate the contents of those cells.

!!! Warning

    Excel stores both date and time data as a single number and flags the cell
    as containing a date/time value. Both date and time are stored in Excel as a
    single number (N: days since the beginning of January 1900). If N < 1 it
    represents a time, if N > 1 it is a date and if N > 1 and has a fractional
    part then it is a datetime.

    However, the cell level formatting of dates and times can mislead you as to
    what is actually stored in a cell and this is a common cause of errors.

    * 0.75 is the time `18:00`, but it could also display in Excel as
      `00/01/1900` or `00/01/1900 18:00` if formatted as a date. This is reasonably
      easy to spot because of the 0th of January!

    * 12.75 is the datetime `12/01/1900 18:00` but could display as the time `18:00` or
      as `12/01/1900`.

    * Unfortunately, integer values like 12 are ambiguous, because Excel doesn't
      differentiate between integer and float numbers: it could just refer to the day
      (the integer 12) or mean exactly midnight on the day (the float 12.0). So it 
      could display perfectly sensibly as `12/01/1900 00:00` or `12/01/1900` but could 
      also display more misleadingly as the time `00:00`.

### Locations

Columns of this type contain location labels showing where the data in the row was
recorded. All of the labels must have been included in the Locations worksheet.

### Taxa

Columns of this type contain taxon names showing the taxon from which other data in the
row was recorded. All of the values in the row must appear in the Taxon Names column in
a Taxa worksheet (i.e. either GBIFTaxa or NCBITaxa).

You don't need to complete any other descriptors. In particular, there is no need to
provide a list of accepted values: the entries are validated against the Taxa worksheet.
This is different from Categorical variables, where a list of category levels is
required (see below)

### Replicate and ID

Both Replicate and ID fields could contain almost any values. Replicates are typically
just shown with repeating numbers, but researchers could use other formats. ID can
represent lots of things (for example, PIT tag numbers for individual organism, fine
scale spatial sampling ID, batch number for reagents) and again could have almost any
format.

So, both ID and Replicate fields are checked for missing data (NAs are permitted) but no
other validation occurs.

### Categorical data

!!! Note

    Field descriptor `levels` required

Both categorical and ordered categorical data (also known as a factors) are made up of a
set of **levels** showing the different groups or treatment. The data in the column then
shows which level applies to each row.

In the `levels` descriptor, you must provide a complete set of all the levels used in
the column, which will be checked against the data. The level names must be short text
labels. **Do not use integer level names**: they are harder to interpret in statistical
analyses and there is a real risk that they are analysed as a number by mistake.

The format is that the level names are separated using semi-colons (`;`). For example:

    Control;Logged;Burned

We automatically remove spaces around level labels, so these also work and might be
easier to read:

    Control; Logged; Burned
    Control ; Logged ; Burned

If the levels aren't obvious, we'd also like label descriptions: they come after each
label, separated by colons (`:`). For example:

    Control:sites in reserve forest;
    Logged:sites in logged forest;
    Burned:sites in burned forest

**Do not** use colons or semi-colons in your level names or descriptions!

For Ordered Categorical fields, the order of the entries in the `levels` descriptor
should be the logical order of the factor. For example, an ordered disturbance gradient
could be:

    Primary:primary rainforest;
    Once:once logged rainforest;
    Twice:twice logged;
    Salvage:salvage logged;
    Oil palm:plantation

### Numeric data

!!! Note

    Field descriptors `method` and `units` required

This field type should be used to record numeric variables **except numeric variables
recorded from taxa** (see Traits below). The `method` descriptor should include
information about how the variable is measured and the `units` descriptor must provide
the units used.

Not all numeric variables have methods or units: a column of replicate numbers, for
example. If this is the case, enter None rather than leaving the descriptors blank. (If
you prefer to use Dimensionless as the unit for dimensionless quantities then that is
also fine!)

### Abundance and trait data

Both traits and abundance data tie a value (category or number) to a single taxon. You
need to format your data so that it is clear which taxon each value comes from. There
are two possible formats:

1. All observations in a column are from **a single taxon**: in this case, you can put a
  valid taxon name (see Taxa worksheet) in the `taxon_name` descriptor for this column.

    _Example: Observation counts in separate columns for each taxon_

    <!-- markdownlint-disable MD013 -->

    |             |                                   |                                   |
    | ----------- | --------------------------------- | --------------------------------- |
    | field_type  | Abundance                         | Abundance                         |
    | taxon_name  | Tiger leech                       | Brown leech                       |
    | method      | Exhaustive search of 50cm quadrat | Exhaustive search of 50cm quadrat |
    | description | quadrat count                     | quadrat count                     |
    | units       | individuals                       | individuals                       |
    | field_name  | tiger_count                       | brown_count                       |
    | 1           | 24                                | 12                                |
    | 2           | 62                                | 3                                 |

    <!-- markdownlint-enable MD013 -->

2. Different rows in the column refer to **different taxa**: in this case, you must also
  have a Taxa column and the `taxon_field` descriptor needs to contain the field name of
  the appropriate Taxa column.

    _Example: Observation counts with different taxa in rows_

    |             |               |                                   |
    | ----------- | ------------- | --------------------------------- |
    | field_type  | Taxa          | Abundance                         |
    | taxon_field |               | common_name                       |
    | method      |               | Exhaustive search of 50cm quadrat |
    | description | Species found | Number found                      |
    | units       |               | individuals                       |
    | field_name  | common_name   | leech_count                       |
    | 1           | Tiger leech   | 24                                |
    | 2           | Brown leech   | 12                                |
    | 3           | Tiger leech   | 62                                |
    | 4           | Brown leech   | 3                                 |

It is an error to provide both `taxon_name` and `taxon_field` descriptors for an
Abundance or Trait field.

### Abundance

!!! Note

    Field descriptors `method`, `units` and **one** of `taxon_name` or `taxon_field`
    required

Abundance is used here as an umbrella term to cover a wide range of
possibilities from casual observation data ('We saw two clouded leopards on
Friday on the road near F100'), through presence/absence data to precise
measurements of abundances or encounter rate.

The `method` descriptor needs to provide a detailed description of the sampling
method, including the area surveyed, the length of time spent sampling, the
number of samplers and any equipment. This should be detailed enough to allow
the sampling protocol to be replicated. If other columns provide sampling
information, such as survey time or area, then make this clear.

### Categorical trait

!!! Note

    Field descriptors `levels` and **one** of `taxon_name` or `taxon_field` required

This is just a categorical variable where the groups apply to a taxa. So, we
need information on the levels used, as for a standard categorical variable, and
a link to taxonomic information as described in the examples above. An example
categorical trait would be the observed life stage of the taxon.

### Numeric trait

!!! Note

    Field descriptors `units`, `method` and **one** of `taxon_name` or `taxon_field` 
    required

A 'trait' in safedata is simply an observation of a property of organisms: that could be
an observation on a single individual, aggregate measures across a sample of individuals
or estimates of population or species level values. A numeric trait is just a numeric
variable where the groups apply to a taxa. So, we need the method and units for the
values, as for a standard numeric variable, and a link to taxonomic information as
described in the examples above. An example numeric trait would be mass.

### Interaction data

Interaction data is data that involves more than one taxa (e.g. observations of clouded
leopard hunting success rate with various prey). This data is essentially just a column
of categorical or numeric data that you want to associate with (at least) two taxon
identities, but there are lots of ways that the taxon identities could be provided.

Interaction data fields do this by using two alternative descriptors to tie the
data to taxa: `interaction_name` and `interaction_field`. Each descriptor can
provide one or more taxon names and optionally their role - the formatting is
identical to the categorical data `levels` descriptors. So, for example an
`interaction_name` descriptor might be:

    Moon rat:prey;Clouded leopard:predator;

You can use one or both of the descriptors, depending on how your data is laid
out. For the most common case of two interacting taxa, the following three
possibilities exist.

1. Both interacting taxa vary from row to row, so taxon names are provided in two fields.

    _Example: Interacting taxa identified in separate columns_

    <!-- markdownlint-disable MD013 -->

    |                   |                    |                    |                            |
    | ----------------- | ------------------ | ------------------ | -------------------------- |
    | field_type        | Taxon              | Taxon              | Categorical Interaction    |
    | interaction_field |                    |                    | predator;prey              |
    | levels            |                    |                    | success;failure            |
    | method            | Visual observation | Visual observation | Visual observation         |
    | description       | Predator observed  | Prey observed      | Outcome of predation event |
    | field_name        | predator           | prey               | outcome                    |
    | 1                 | Clouded leopard    | Brown rat          | success                    |
    | 2                 | Flat headed cat    | Moon rat           | failure                    |

    <!-- markdownlint-enable MD013 -->

2. Alternatively, all of the data might refer to the same two taxa, so the taxon names
    can be provided directly.

    _Example: Interacting taxa constant_

    |                  |                                         |
    | ---------------- | --------------------------------------- |
    | field_type       | Categorical Interaction                 |
    | interaction_name | Clouded leopard:predator;Brown rat:prey |
    | levels           | success;failure                         |
    | method           | Visual observation                      |
    | description      | Outcome of predation event              |
    | field_name       | outcome                                 |
    | 1                | success                                 |
    | 2                | failure                                 |

3. Finally, one side of the interaction might vary from row to row but the other side is
    constant for all rows.

    _Example: Interacting taxa identified by name and by column_

    |                   |                    |                            |
    | ----------------- | ------------------ | -------------------------- |
    | field_type        | Taxon              | Categorical Interaction    |
    | interaction_name  |                    | Clouded leopard:predator;  |
    | interaction_field |                    | prey:prey species;         |
    | levels            |                    | success;failure            |
    | method            | Visual observation | Visual observation         |
    | description       | Prey observed      | Outcome of predation event |
    | field_name        | prey               | outcome                    |
    | 1                 | Brown rat          | success                    |
    | 2                 | Moon rat           | failure                    |

You must provide at least two taxon names or fields, but you can provide more if
you have tritrophic interactions! Again, you can use any combination of
interaction names and fields to provide your taxon identities.

### Categorical interactions

!!! Note

    Field descriptors `levels` and `interaction_name` and/or `interaction_field` required

### Numeric interactions

!!! Note

    Field descriptors `units`, `method` and `interaction_name` and/or `interaction_field`
    required

### File

!!! Note

    Field descriptor `file_container` may be required

This type is used to provide references to information stored in external files.
It can be used in two ways:

1. The values in the data are direct references to external files provided in the
  Summary worksheet. The values in the field are checked against the list of external
  file names and they must all appear there.  For example, if the Summary worksheet
  includes an external file row with `My_raster_1.tiff` and `My_raster_2,tiff`:

    |             |                          |                             |
    | ----------- | ------------------------ | --------------------------- |
    | field_type  | Numeric                  | File                        |
    | description | Altitude in metres       | DEM  file used for altitude |
    | method      | Extracted from DEM tiffs |                             |
    | taxon_name  |                          |                             |
    | units       | metres                   |                             |
    | field_name  | Altitude                 | DEM                         |
    | 1           | 100                      | My_raster_1.tiff            |
    | 2           | 200                      | My_raster_1.tiff            |
    | 3           | 300                      | My_raster_1.tiff            |
    | 4           | 400                      | My_raster_1.tiff            |
    | 5           | 500                      | My_raster_2.tiff            |
    | 6           | 600                      | My_raster_2.tiff            |
    | 7           | 700                      | My_raster_2.tiff            |
    | 8           | 800                      | My_raster_2.tiff            |

2. The values in the data are files contained within an external file. In this case, the
  descriptor `file_container` is used to check the external file is present, but the
  values in the field are not checked. For example:

    |                |                          |                    |
    | -------------- | ------------------------ | ------------------ |
    | field_type     | Numeric                  | File               |
    | description    | Altitude in metres       | Image of Quadrat   |
    | method         | Extracted from DEM tiffs |                    |
    | taxon_name     |                          |                    |
    | units          | metres                   |                    |
    | file_container |                          | My_archive.zip     |
    | field_name     | Altitude                 | Quadrat_image      |
    | 1              | 100                      | Site_quadrat_1.jpg |
    | 2              | 200                      | Site_quadrat_2.jpg |
    | 3              | 300                      | Site_quadrat_3.jpg |
    | 4              | 400                      | Site_quadrat_4.jpg |
    | 5              | 500                      | Site_quadrat_5.jpg |
    | 6              | 600                      | Site_quadrat_6.jpg |
    | 7              | 700                      | Site_quadrat_7.jpg |
    | 8              | 800                      | Site_quadrat_8.jpg |

### Comments

If you have a free text field with notes or comments, then this is the field
type to use. We don't validate anything in comments fields: they're not
expected to be complete data and you can put anything in them.

A word of caution though: it is **highly unlikely** that anyone will ever read
your comments column again. If there is genuinely important information that
might apply across multiple rows, consider coding it as an explicit variable
rather than consigning it to a comments field.

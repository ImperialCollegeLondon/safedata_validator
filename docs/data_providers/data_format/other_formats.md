# Data in other file formats

We understand that Excel may not be everyone's preferred option but:

* it is extremely commonly used,
* it is widely available,
* the format is well supported in a range of programming languages, and
* it is well-suited to storing tabular data frames, as long as the data is well
  structured and validated.

However, you can include files in any format in a dataset published through the
`safedata` system. These could include zipped archives containing raw data that has been
used to construct the analysis data, GIS files, or really anything that is included in
the dataset.

Even if all of your actual data is in other file formats, you will still need to fill in
an Excel template containing the key metadata for your dataset:

* Fill in a standard Excel template, using the Summary, Taxa (GBIFTaxa and/or NCBITaxa)
  and Locations worksheets to provide basic metadata and spatial and taxonomic indexing.
* In the Summary sheet, use the **External file** and **External file description** rows
  to provide a name and description for **every** external file you want to upload.
* Use the standard data worksheets to provide field descriptions for any tabular data
  held in external files as described below.

Some key points:

* File names **must not** contain spaces and you should use the **exact** file name of
  the data file in the metadata description.
* You **must** also complete the Taxa and Locations worksheets if applicable. We will
  not be able to validate them against the external data but we will be able to add the
  taxa and locations to our data index.
* You can mix standard Excel data worksheets and references to external files. In this
  way you can provide bulk raw data and processed data tables.
* If you have have many small files, then zipping them up and providing a broad
  description is preferable.

## Tabular data in other formats

We would prefer data tables to be included as Excel worksheets but there are good
reasons not to use Excel for tabular data:

* **Very large data tables**. Excel struggles with very large numbers of rows and
  reading data from large Excel files can be very slow. We're happy to accept very large
  tables submitted in other commonly accessible formats as simple text files, compressed
  text files, R data files and the like.

* **Structured databases**. If your dataset is a SQL database or something similar, with
  table definitions, constraints, keys and the like, then forcing it into an Excel
  workbook is a poor idea. We are happy to accept database dump files, but these must
  be text dumps so that they are portable.

If you do provide tabular data in another format, then you will still need to provide an
Excel file to provide metadata and this should include completing a data worksheet
describing the table metadata. The data worksheet does not need to contain any of the
actual data but it will allow us to index the available fields. Briefly:

* Create a worksheet for the table and [describe the fields](data.md).
* In the [data worksheet block](summary.md#the-data-worksheet-block) on the Summary
  sheet, use the **Worksheet external file** row to give the name of the file containing
  the tabular data. This must appear in the list of external file names and description.
* If an external file contains more than one data table, then it is fine to refer to the
  same file in more than one data worksheet.

## Very large files

Zenodo are pretty generous with their filespace but there are limits. An obvious example
of very large supporting datasets are bulk media files like audio and camera trap
images. These may need long-term bulk storage outside of Zenodo and we are currently
exploring options to host such data.

## Submitting data in other file formats

The data submission process through the `safedata` system is only set up to validate and
publish correctly formatted Excel files and we use the contents of these files to index
the datasets and provide dataset descriptions. As files in other formats are described
within the formatted Excel files, they can be uploaded directly to Zenodo - this is
likely to be much faster and much kinder to the metadata server.

However, as we discuss [here](../availability.md#data-administration), all data
submitted to your project's Zenodo community must be published using its central
administrative account. The process to include other file formats in a dataset is
therefore:

1. The data provider submits their Excel file and provides access to the other data files.
2. The project's data administrator then creates a new Zenodo deposit using the
   administrative account, uploads the other data files and then publishes the dataset
   through the metadata server to link the information in the submitted Excel file to the
   Zenodo deposit.

Recipes for administrators to create new datasets and to update existing data can be
found [here](../../data_managers/using_safedata/overview.md).

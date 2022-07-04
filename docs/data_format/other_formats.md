# Data in other file formats

We understand that Excel may not be everyone's preferred option but:

* it is extremely commonly used,
* it is widely available,
* the format is well supported in a range of programming languages, and
* it is well-suited to storing tabular data frames, as long as the data is [well
  structured and validated](../usage/usage.md).

However, you can include files in any format in a dataset published through SAFE. These
could include zipped archives containing raw data that has been used to construct the
analysis data, GIS files, or really anything that is included in the dataset.

Even if all of your actual data is in other file formats, you will still need to fill in
an Excel template containing the key metadata for your dataset:

* Fill in a standard Excel template, using the Summary, Taxa and Locations worksheets to
  provide basic metadata and spatial and taxonomic indexing.
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
  workbook is a poor idea. We are happy to accept database dump files, but we these must
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
images. These will need long-term bulk storage outside of Zenodo and we are currently
exploring options to host such data.

## Submitting data in other file formats

The data submission process through the SAFE project is only set up to validate and
publish correctly formatted Excel files and we use the contents of these files to index
the datasets and provide dataset descriptions. As files in other formats are described
within the formatted Excel files, they can be uploaded directly to Zenodo - this is
likely to be much faster and much kinder to our web server.

However, as we discuss [here](../availability.md#data-administration), all data
submitted to the SAFE Zenodo community must be published using our central
administrative account. The process to include other file formats in a dataset is
therefore:

1. The data provider submits their Excel file and provides access to the other data files.
2. A SAFE dataset administrator then creates a new Zenodo deposit using the
   administrative account, uploads the other data files and then publishes the dataset
   through the SAFE website to link the information in the submitted Excel file to the
   Zenodo deposit.

The details for administrators are as follows:

### New datasets

If this is the first version of a dataset to be published, an administrator will need to
do the following:

* When the Excel metadata file passes checks, log in to Zenodo using our Zenodo curation
  account (`the_safe_project`, tied to info@safeproject.net).
* Navigate to the [Zenodo deposit page](https://zenodo.org/deposit) and click on 'New
  upload'.
* On the 'New upload' page, upload the external files and fill in the following fields:
  * Upload type: set to dataset.
  * Title, orders and description: these will be overwritten when the dataset is
    published, so these field just need to be non-blank so that the draft can be saved.
    It is fine to type `test` in all three!
* Click 'Save' (**NOT** 'Publish') and then note down the deposit number for the new
  draft that appears in the URL.
* Go to the [dataset administration
  page](https://safeproject.net/datasets/administer_datasets) at the SAFE project
  website, click the ‘Adopt’ button and then paste the deposit ID number in.
* The website will then check the filenames match up, fills in the Zenodo description
  and publishes the dataset.

### Updates

If this is an update to an existing dataset, then the process is broadly similar:

* When the Excel metadata file passes checks, log in to Zenodo using our Zenodo curation
  account (`the_safe_project`, tied to info@safeproject.net).
* Do not create a new upload. Instead find the published dataset and click on the 'New
  version' button.
* Upload any new or changed files and delete any outdated files. You **must** delete the
  existing Excel file, which will be replaced automatically by the new one.
* You do not need to set any of the other fields - they will have been filled in when
  the first version was published and will be refreshed from the new metadata version.
* Again, click 'Save' and note down the new version deposit ID.
* Adopt the deposit from the SAFE website as for a new dataset.

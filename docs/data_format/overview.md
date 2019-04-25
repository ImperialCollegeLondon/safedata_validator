# Data formatting

In most cases, we expect the data format to be provided as **tabular datasets stored in Excel
spreadsheets**. This accounts for the vast majority of the data files used by researchers. We do
support other data files but those files will need to be documented using the same Excel format.

!!! Warning
    
    The formatting details described in this documentation will be used to automatically publish your
    data to Zenodo. You should choose titles, descriptions and keywords that you would be happy to be
    **permanently** associated with your dataset!


# Excel format overview

The basic format for a SAFE dataset submission is an Excel Workbook, which will typically contain at least four worksheets. The first  three worksheets must use the standard names listed below. 

* [**Summary**](summary.md): This contains some simple information about the authors of the dataset, access rights and the individual data tables in the dataset. 
* [**Taxa**](taxa.md): This describes all the taxa used in the dataset. 
* [**Locations**](locations.md): This describes all the sampling locations used in the dataset.
* [**Data worksheets**](data.md): After these worksheets come your data tables. You should label these sheets with a sensible name (not 'Sheet1'!) and each data table must be described in the Summary worksheet. You can include as many data tables as you like in a single dataset: we don't want you to spend time rearranging your data and are happy just to take the data in the natural tables you already use.


## File naming

Use a simple short name for your spreadsheet - there will be a lot of information giving more
detail inside. Please **do not use spaces** in your file name - you can use underscores to
separate words.

## Spreadsheet Template and Examples

Click on this link to download the {{ :working_at_safe:template.xlsx |spreadsheet template }}
containing the required worksheets, labels and headers.

You can also look at existing published datasets to see how the format is used:

* [https://safeproject.net/datasets/view_datasets](https://safeproject.net/datasets/view_datasets)
* [https://zenodo.org/communities/safe/](https://zenodo.org/communities/safe/)

## Format checking

We've tried to make the description below as clear as possible but in order to help you prepare
your file:

* It is easier to follow an example than to follow a description, so please use the template and look at the examples. 
* We use a Python program to automatically check the formatting of datasets. When you submit a file, you will get a report back from this program that will highlight any problems with your dataset. If there are problems, fix them and replace the submitted file. Once the file passes through the checker without problem, we will double check the file and then publish your dataset.



## FAQ or 'What does this error mean?'

There is a separate page [[working_at_safe:data_format_faq|here]] for frequently asked
questions about issues with data formatting.

## Non-Excel data files 

If
you have data that is stored in other sorts of files, because of formatting or file size
issues, then you will still need to fill in an Excel template containing the key metadata for
your dataset. We will then upload your data files manually to Zenodo and link the metadata to
them.


Most datasets will be submitted as Excel workbooks - the vast majority of data
is submitted as Excel files (or is in some other spreadsheet format that could be saved as Excel) - but we also support the inclusion of other data files. Bulk media files are often very large and will need long bulk storage outside of Zenodo.



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

## Tabular data in external files 

Some possible external data files contain tabular data: for example, a zipped file containing a
large CSV dataset, a SQL data dump or a shapefile with an attribute table. If so, then please
do fill in a data worksheet describing the table metadata as described below. The data
worksheet does not need to contain any of the actual data but it will allow us to index the
available fields.



 * Create a worksheet for the table and describe the fields as described below. 
 * In the worksheet summary fields, use the ''Worksheet external file'' row to give the name of the file containing the tabular data. This must appear in the list of external file names and description. 
 * If an external file contains more than one data table, then it is fine to refer to the same file in more than one data worksheet.



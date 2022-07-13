# Data formatting

In most cases, we expect the data format to be provided as **tabular datasets stored in
Excel spreadsheets**. This accounts for the vast majority of the data files used by
researchers. We do support other data files but those files will need to be documented
using the same Excel format.

!!! Warning
    The formatting details described in this documentation will be used to automatically
    publish your data to Zenodo. You should choose titles, descriptions and keywords
    that you would be happy to be **permanently** associated with your dataset!

## Excel format overview

The basic format for a SAFE dataset submission is an Excel Workbook, which will
typically contain at least four worksheets. The first four worksheets must use the
standard names listed below. While both Taxa and NCBITaxa can be included, it is only
required to include one of them.

- [**Summary**](summary.md): This contains some simple information about the authors of
  the dataset, access rights and the individual data tables in the dataset.
- [**GBIFTaxa**](gbif_taxa.md): This describes taxa used in the Data worksheets, using
  the GBIF taxonomy backbone as a reference. This is best for observational data.
- [**NCBITaxa**](ncbi_taxa.md): This describes taxa used in the Data worksheets, using
  the NCBI taxonomy database as a reference. This is best for sequencing data.
- [**Locations**](locations.md): This describes all the sampling locations used in the
  dataset.
- [**Data worksheets**](data.md): After these worksheets come your data tables. You
  should label these sheets with a sensible name (not 'Sheet1'!) and each data table
  must be described in the Summary worksheet. You can include as many data tables as you
  like in a single dataset: we don't want you to spend time rearranging your data and
  are happy just to take the data in the natural tables you already use.

### File naming

Use a simple short name for your spreadsheet - there will be a lot of information giving
more detail inside. Please **do not use spaces** in your file name - you can use
underscores to separate words.

### Spreadsheet Template and Examples

Click on this link to download the [spreadsheet
template](https://github.com/ImperialCollegeLondon/safedata_validator/blob/master/test_files/Template.xlsx?raw=true)
containing the required worksheets, labels and headers.

You can also look at existing published datasets to see how the format is used:

- [https://safeproject.net/datasets/view_datasets](https://safeproject.net/datasets/view_datasets)
- [https://zenodo.org/communities/safe/](https://zenodo.org/communities/safe/)

# Overview

Part of the agreement for research projects working at the SAFE Project is that all project
data are submitted to the SAFE Project data repository, so that it is available to future
researchers. In order to make it easy for data to be found and used in the future, we need
researchers to provide some (relatively!) simple metadata information in their datafiles.

The process in overview is:

* Prepare your data following the formatting guide.
* Go to the SAFE Project website and [Submit your dataset](https://www.safeproject.net/datasets/submit_dataset).
* We will automatically validate the data formatting. 
* If the validation succeeds then we will publish it on Zenodo.
* If the validation fails then you will get an error report so you can fix the problems and resubmit.

!!! Note

    You **must not** publish your dataset directly to Zenodo. This skips the validation step and means
    that the datasets are not linked together under a common account.

This documentation describes the format required for datasets collected at the SAFE project and the
software that we have developed to validate formatted data.

## Authorship and funding

All datasets submitted to SAFE are published to the [Zenodo data
repository](https://zenodo.org/communities/safe). Zenodo issues DOIs for all data depositions, making your
datasets easily citable.

You will need to provide a list of authors for submitted datasets, which will form part of the permanent
citation for the published dataset. Authorship on published datasets should be treated in the same way as
you would consider authorship on papers: you should include not only the people responsible for physically
collecting the data but also other researchers who facilitated the work, such as project supervisors and
local collaborators. Similarly, you must acknowledge any funding that supported your research.

## Format validation

The `safe_dataset_checker` repository contains a Python module to validate submitted files and report on
any problems. Most datasets will be submitted as Excel workbooks - the vast majority of data
is submitted as Excel files (or is in some other spreadsheet format that could be saved as Excel) - but we also support the inclusion of other data files. Bulk media files are often very large and will need long
term bulk storage outside of Zenodo.

The code validates:

  1. The data submission formatting of the file.
  1. All taxonomic names against the GBIF taxonomy database.
  1. All location names against the SAFE Gazetteer.

Datasets can be submitted by registered researchers at the [SAFE
website](https://safeproject.net/datasets/submit_dataset) which will automatically use this code to check
that the file is formatted correctly. However, you can also download and run it yourself!

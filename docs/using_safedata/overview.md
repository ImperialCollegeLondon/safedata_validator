# Using the `safedata` system

The `safedata_validator` package is one component of the wider `safedata` system for
data management and discovery. This system comprises:

1. The `safedata_validator` package itself, which is used to validate submitted datasets
   and ensure that the data and metadata for those datasets are consistent and meet the
   minimum requirements. When a dataset is successfully validated, the package is used
   to both publish the dataset to the Zenodo community for your datasets and to upload
   the metadata to a seperate metadata server for the project.

2. The metadata server is a web server running the [`safedata_server` server web
   application](https://github.com/ImperialCollegeLondon/safedata_server). This provides
   an index of the published datasets along with a range of APIs to search the metadata,
   including text, taxonomy and spatial searches of the published datasets.

3. The project Zenodo community: this is a project specific grouping of Zenodo records
   which provides DOIs and download access for the actual data files.

4. The [`safedata` R
   package](https://cran.r-project.org/web/packages/safedata/index.html), which is an R
   package that makes it easy for users to discover and download datasets of interest
   from your community.

## Use recipes

The sections below provide examples of using the `safedata_validator` package to
administer datasets. Typical use will be from the command line using a Unix-like shell
or [Windows subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install),
but the examples also show how to use the programmatic API for `safedata_validator` from
within Python.

<!--- 
The include extension does not seem to work when including content within tabbed blocks
in the use case recipes below so for the moment keeping the code in sections not tabs
and including the code from source files that can be validated by pre-commit
-->

### Validating and publishing a dataset

These examples show the typical workflow for publishing a dataset and accompanying
metadata using `safedata_validator`.

#### Shell

```sh
{!docs/using_safedata/publish_script.sh!}
```

#### Python

```sh
{!docs/using_safedata/publish_script.py!}
```

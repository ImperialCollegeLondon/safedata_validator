# Managing datasets with `safedata_validator`

The `safedata_validator` package is one component of the wider `safedata` system for
data management and discovery. This system comprises:

1. The [`safedata_validator` Python
   package](https://pypi.org/project/safedata-validator) itself, which is used to
   validate submitted datasets and ensure that the data and metadata for those datasets
   are consistent and meet the minimum requirements. When a dataset is successfully
   validated, the package also provides tools to both publish the dataset to the Zenodo
   community for your datasets and to upload the metadata to a seperate metadata server
   for the project.

2. The metadata server is a web server running the [`safedata_server` web
   application](https://github.com/ImperialCollegeLondon/safedata_server). This provides
   an index of the published datasets along with a range of APIs to search the metadata,
   including text, taxonomy and spatial searches of the published datasets.

3. The project Zenodo community: this is a project specific grouping of Zenodo records
   which provides DOIs and download access for the actual data files. Each project using
   the `safedata` system will have it's own **separate** Zenodo community.

4. The [`safedata` R
   package](https://cran.r-project.org/web/packages/safedata/index.html), which is an R
   package that makes it easy for users to discover and download datasets of interest
   from your community.

## Installing and using `safedata_validator`

Installing and configuring the `safedata_validator` package has multiple steps. The
basic overview is:

- [x] Ensure that you have a recent version of Python installed on your computer and
   install the `safedata_validator` package from
   [PyPi](https://pypi.org/project/safedata-validator/) using the `pip` package
   installer tool.

- [x] Create a `safedata_validator` configuration file, which is used by the package and
   command line tools to locate required resources and settings.

- [x] Use the `safedata_build_local_gbif` and `safedata_build_local_ncbi` command line
   tools to create local taxonomic validation databases and add the locations of these
   files to the configuration.

- [x] Create a GeoJSON gazetteer for your project, defining named locations to be used
   across datasets, and add the location of the gazetteer file to the configuration.

At this point, you should be able to use the `safedata_validate` tool to validate
datasets. However, there are extra steps to allow datasets to be
published.

- [x] Create a Zenodo account, community and access token that will be used to publish
  validated datasets and add these details to your configuration.

- [x] Set up a `safedata_server` metadata server to provide a searchable API for the
  detailed dataset metadata and again add the details to your configuration.

Once you have installed and configured these tools, then you can use the provided
[command line tools](command_line_tools/overview.md) to validate and publish datasets.
The [usage recipes](using_safedata/overview.md) show how the tools are used.

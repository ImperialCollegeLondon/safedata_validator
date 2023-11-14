# Using the `safedata` system

<!-- markdownlint-disable MD046 -->

The `safedata_validator` package is one component of the wider `safedata` system for
data management and discovery. This system comprises:

1. The `safedata_validator` package itself, which is used to validate submitted datasets
   and ensure that the data and metadata for those datasets are consistent and meet the
   minimum requirements. When a dataset is successfully validated, the package is used
   to both publish the dataset to the Zenodo community for your datasets and to upload
   the metadata to a seperate metadata server for the project.

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

## Use case recipes

The sections below provide examples of using the `safedata_validator` package to
administer datasets. Typical use will be from the command line using a Unix-like shell
or [Windows subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install),
but the examples also show how to use the programmatic API for `safedata_validator` from
within Python. The examples assume that a user has provided a SAFE formatted dataset and
an additional ZIP file containing additional files:

* `SAFE_dataset.xlsx`
* `Supplementary_files.zip`

### Validating and publishing as a new dataset

These examples show the typical workflow for publishing these data and accompanying
metadata as a completely new dataset using `safedata_validator`.

=== "bash"

    ```sh
    {%
    include "data_managers/using_safedata/publish_script.sh"
    %}
    ```

=== "python"

    ```python
    {%
    include "data_managers/using_safedata/publish_script.py"
    %}
    ```

### Validating and publishing as an update

A new version of an existing dataset can be created by requesting a new Zenodo record
using the **Concept ID** of an existing dataset. This is a Zenodo record ID that
identifies a collection of versions of a dataset. In the previous example, the specific
record ID for the dataset was `1143714` - the concept record ID might be `1143713`.

The workflows below look almost identical, except for the initial step of creating the
deposit using an existing concept ID.

=== "bash"

    ```sh
    {%
    include "data_managers/using_safedata/update_script.sh"
    %}
    ```

=== "python"

    ```python
    {%
    include "data_managers/using_safedata/update_script.py"
    %}
    ```

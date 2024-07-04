# Publishing a dataset

<!-- markdownlint-disable MD046 MD033 -->

This page provide examples of using the `safedata_validator` package to publish
validated datasets. The examples assume that a user has provided a SAFE formatted
dataset and an additional ZIP file containing additional files:

* `Example.xlsx`
* `Supplementary_files.zip`

Both of the use cases below include the creation and includsion of a GEMINI compliant
XML metadata file in a published dataset. We recommend this as good practice, but it is
optional.

## Validating and publishing as a new dataset

The [`safedata_zenodo publish_dataset`
command](../command_line_tools/safedata_zenodo.md#the-publish_dataset-subcommand) is
the main function for publishing a dataset. The example below shows it being used to
publish a dataset and additional external files:

```sh
safedata_zenodo publish_dataset Example.xlsx -e  Supplementary_files.zip
```

That command packages up a set of operations needed to publish a dataset. The example
code in the tabs below show the underlying workflow, either using the command line
interface or working from within Python,  for publishing these data and accompanying
metadata as a completely new dataset using `safedata_validator`. You would typically not
need to use these individual commands: this information is here to show what is going on
under the hood.

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

## Validating and publishing as an update

Zenodo can hold multiple versions of a dataset, allowing you to publish updates and
corrections. When you create a new dataset, that dataset is assigned a record id but
is also automatically also given a **concept record id**. These are often sequential:
in the previous example, the specific record ID for the dataset was `1143714` - the
concept record ID might be `1143713`. The concept record id is used to group versions of
a dataset and the associated DOI is used to redirect to the most recent version of the
dataset.

Releasing a new version generates a new record on Zenodo that is simply copied from the
most recent version of a dataset, updates the files and metadata and creates a new
record grouped under the same concept id. In order to do this, you can simply provide
the record id of the most recent version to the data publication commands. Note that the
required record id not the concept ID, but the most recent version of the dataset being
updated.

Again, the most straightforward approach is to use the `publish_dataset` subcommand and
add the `--new-version` (or `-n`) argument.

```sh
safedata_zenodo publish_dataset Example.xlsx -e  Supplementary_files.zip \
    --new-version 1143714
```

As above, the tabs below show what is going on within that process.

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

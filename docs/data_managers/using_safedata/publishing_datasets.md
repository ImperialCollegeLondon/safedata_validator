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
command](../command_line_tools/safedata_zenodo.md#the-publish_dataset-subcommand) is the
main function for publishing a dataset. It is important to note that the dataset
metadata file must be provided first and the dataset file provided second, otherwise the
publication process will fail. The example below shows it being used to publish a
dataset and additional external files:

```sh
safedata_zenodo publish_dataset Example.json Example.xlsx \
    --external-file Supplementary_files.zip
```

The expected output from that command is shown below:

```sh
- Configuring Resources
    - Configuring resources from user config: configs/config.cfg
    - Validating gazetteer: spatial_resources/gazetteer.geojson
    - Validating location aliases: spatial_resources/location_aliases.csv
    - Validating GBIF database: gbif_databases/gbif_backbone_2021-11-26.sqlite
    - Validating NCBI database: ncbi_databases/ncbi_taxonomy_2023-11-01.sqlite
    - Validating project database: project_databases/safe_projects.csv
Deposit created: 1143714
XML created: 1143714_GEMINI.xml
Uploading files:
Uploading Example.xlsx
100%|███████████████████████████████████████| 160k/160k [00:00<00:00, 494kB/s]
Uploading Supplementary_files.zip
100%|████████████████████████████████████| 1.00k/1.00k [00:00<00:00, 3.47kB/s]
Uploading 1143714_GEMINI.xml
100%|████████████████████████████████████| 27.1k/27.1k [00:00<00:00, 83.4kB/s]
Uploading deposit metadata
Dataset published: https://zenodo.org/records/1143714
```

The `publish_dataset` subcommand packages up a set of operations needed to publish a
dataset. The example code in the tabs below show the underlying workflow, either using
the command line interface or working from within Python,  for publishing these data and
accompanying metadata as a completely new dataset using `safedata_validator`. You would
typically not need to use these individual commands: this information is here to show
what is going on under the hood.

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

## Validating and publishing a new version of a dataset

Zenodo can hold multiple versions of a dataset, allowing you to publish updates and
corrections. Each version of a dataset will have a different record ID and they are also
grouped together under a shared **concept record ID**. One of the versions is always
identified as the **latest version** - and actually the concept ID works as a DOI that
always redirects to the latest version.

When you create a new version of a dataset, the Zenodo system creates an exact copy of
the most recent version. Users can then update any files that need changing, remove
outdated files, and update the metadata for the new deposit before publishing it.

In order to do this, you can provide the record ID of the most recent version of a
dataset that you want to update. The most straightforward approach is to use the
`publish_dataset` subcommand and add the `--new-version` (or `-n`) argument.

```sh
safedata_zenodo publish_dataset Example.json Example.xlsx \
    --external-file  Supplementary_files.zip \
    --new-version 1143714
```

The output from that command would look like:

```sh
- Configuring Resources
    - Configuring resources from user config: configs/config.cfg
    - Validating gazetteer: spatial_resources/gazetteer.geojson
    - Validating location aliases: spatial_resources/location_aliases.csv
    - Validating GBIF database: gbif_databases/gbif_backbone_2021-11-26.sqlite
    - Validating NCBI database: ncbi_databases/ncbi_taxonomy_2023-11-01.sqlite
    - Validating project database: project_databases/safe_projects.csv
Preparing new version of deposit 1143714
 - Unchanged files: Supplementary_files.zip
 - Removing outdated files: 1143714_GEMINI.xml, Example.xlsx
 - Uploading new or updated files: Example.xlsx
Deposit created: 1143900
XML created: 1143900_GEMINI.xml
Removing outdated files: 1143714_GEMINI.xml, Example.xlsx
Uploading files:
Uploading Example.xlsx
100%|███████████████████████████████████████| 160k/160k [00:00<00:00, 479kB/s]
Uploading 1143900_GEMINI.xml
100%|████████████████████████████████████| 27.1k/27.1k [00:00<00:00, 89.0kB/s]
Uploading deposit metadata
Dataset published: https://sandbox.zenodo.org/records/1143900
```

The `publish_dataset` subcommand does more complex checking when creating a new version
of an existing dataset. Because the newly created deposit already contains copies of the
most recent files, the command needs to check for:

* completely new files to be uploaded,
* existing files where the content has changed and which should be updated,
* existing files that are not in the new publication request and which should be
  deleted, and
* existing files that have not changed and can be left as is.

The subcommand will fail under a few circumstances:

* The provided record ID is not the most recent version. The command automatically
  checks for most recent version of the provided ID and will stop if the provided
  version does not match the most recent version. It will print out what that most
  recent ID is, but it does not automatically assume this is what you meant!

* If the files that you have provided to publish are identical to the existing files on
  the most recent version, then it will stop to avoid creating duplicate identical
  deposits. This step checks the name and the MD5 hash of each file against the existing
  file. The hash provides a unique signature for the contents of a file that allows the
  code to test for identical files. Using a different name for the file is currently
  accepted as a change but we don't advise publishing identical files under different
  names!
  
  This check ignores any GEMINI XML file in most recent version. Since these are named
  using the record ID of the deposit and are generated from the data files, they will
  only differ in their file name.

As above, the tabs below show what is going on within that process. This is more
involved than creating a new dataset because the existing files need to be deleted.

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

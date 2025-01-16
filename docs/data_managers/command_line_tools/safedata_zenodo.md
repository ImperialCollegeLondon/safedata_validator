# Publishing data on Zenodo

The process of publishing a dataset to the [Zenodo](https://zenodo.org) data repository
requires a validated dataset and the subcommands of the `safedata_zenodo` tool. These
commands use two _different_ JSON metadata files.

## Dataset metadata

The [`safedata_validate` command](./safedata_validate.md) generates a JSON file
containing a standard JSON description of the metadata in the dataset and of the data
tables it contains. Some of this metadata is used to populate the Zenodo description of
the published dataset files.

The dataset metadata is _also_ used to populate the database of a **metadata server**.
This is a separate website that provides the API for searching available data and forms
the main data discovery backend for the `safedata` R package.

### Zenodo deposit metadata

The Zenodo API uses JSON metadata to return key details on a Zenodo
deposit that is being prepared or published. It contains key API links that are
used to provide file details. This file will be generated when a new deposit is
generated and is then used to carry out the other publication steps.

## The `safedata_zenodo` tool

!!! info
    The subcommands of the `safedata_zenodo` tools require that the `zenodo`
    section of the [resources configuration](../install/configuration.md#publication-configuration)
    be completed. This is not required for simply validating datasets.

The `safedata_zenodo` command line tool provides the following subcommands which
are used to publish data, post metadata and help maintain and document published
datasets.

The top level command line help for the tool, showing the available subcommands
is shown below:

```bash
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_top.txt"
%}
```

### The `safedata_zenodo` subcommands

The command line help for each of the various subcommands is shown below:

#### The `publish_dataset` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_publish_dataset.txt"
%}
```

#### The `create_deposit` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_create_deposit.txt"
%}
```

#### The `discard_deposit` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_discard_deposit.txt"
%}
```

#### The `get_deposit` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_get_deposit.txt"
%}
```

#### The `publish_deposit` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_publish_deposit.txt"
%}
```

#### The `upload_files` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_upload_files.txt"
%}
```

#### The `delete_files` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_delete_files.txt"
%}
```

#### The `upload_metadata` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_upload_metadata.txt"
%}
```

#### The `amend_metadata` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_amend_metadata.txt"
%}
```

#### The `sync_local_dir` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_sync_local_dir.txt"
%}
```

#### The `maintain_ris` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_maintain_ris.txt"
%}
```

#### The `generate_html` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_generate_html.txt"
%}
```

#### The `generate_xml` subcommand

See also [here](../install/configuration.md#xml-configuration).

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_generate_xml.txt"
%}
```

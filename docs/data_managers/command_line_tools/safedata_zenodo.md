# Zenodo and metadata tools overview

The `safedata_validator` package also contains tools to publish validated
datasets to the [Zenodo](https://zenodo.org) data repository and to update
the metadata server.

## Data publication process

The process of publishing a dataset involves both the `safedata_validate` tool
and the various subcommands of the `safedata_zenodo` tool. These commands make
use of two _different_ JSON metadata files.

### Dataset metadata

The `safedata_validate` command generates a JSON file containing a standard JSON
description of the metadata in the dataset and of the data tables it contains.
Some of this metadata is used to populate the Zenodo description of the published
dataset files.

The dataset metadata is _also_ used to populate the database of a **metadata
server**. This is a separate website that provides the API for searching
available data and forms the main data discovery backend for the `safedata` R
package.

### Zenodo deposit metadata

The Zenodo API returns JSON metadata that provides key details on a Zenodo
deposit that is being prepared or published. It contains key API links that are
used to provide file details.

## The `safedata_zenodo` tool

!!! info
    The subcommands of the `safedata_zenodo` tools require that the `zenodo`
    section of the [resources configuration](../install/configuration.md#zenodo) be
    completed. This is not required for simply validating datasets.

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

### Simple publication process

As an initial example, the process for publishing a simple dataset
(without any external data files) would be:

```sh
# Validate the file, creating the Test_format_good.json metadata file
safedata_validate Test_format_good.xlsx

# Create a new deposit, creating a JSON file of metadata for the Zenodo
# deposit as - for example - zenodo_1059375.json
safedata_zenodo create_deposit

# Upload the file
safedata_zenodo upload_file zenodo_1059375.json Test_format_good.xlsx

# Populate the Zenodo deposit metadata from the dataset metadata
safedata_zenodo upload_metadata zenodo_1059375.json Test_format_good.json

# Publish the deposit
safedata_zenodo publish_deposit zenodo_1059375.json

# Post the metadata to the metadata server
safedata_zenodo post_metadata zenodo_1059375.json Test_format_good.json

```

### The `safedata_zenodo` subcommands

The command line help for each of the various subcommands is shown below:

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

#### The `upload_file` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_upload_file.txt"
%}
```

#### The `delete_file` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_delete_file.txt"
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

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_zenodo_generate_xml.txt"
%}
```

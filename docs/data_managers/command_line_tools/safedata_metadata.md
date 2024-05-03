
# The `safedata_metadata` tool

!!! info
    The subcommands of the `safedata_metadata` tools require that the `metadata`
    section of the [resources configuration](../install/configuration.md#zenodo) be
    completed. This is not required for simply validating datasets.

The `safedata_metadata` command line tool provides the following subcommands which
are used to update a server running the safedata server API.

The top level command line help for the tool, showing the available subcommands
is shown below:

```bash
{%
include "data_managers/command_line_tools/command_line_usage/safedata_metadata_top.txt"
%}
```

## The `safedata_metadata` subcommands

The command line help for each of the various subcommands is shown below:

### The `post_metadata` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_metadata_post_metadata.txt"
%}
```

### The `update_resources` subcommand

```sh
{%
include "data_managers/command_line_tools/command_line_usage/safedata_metadata_update_resources.txt"
%}
```

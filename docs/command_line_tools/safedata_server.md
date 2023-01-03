
# The `safedata_server` tool

!!! info The subcommands of the `safedata_server` tools require that the `metadata`
    section of the [resources configuration](../install/configuration.md#zenodo) be
    completed. This is not required for simply validating datasets.

The `safedata_server` command line tool provides the following subcommands which
are used to update a server running the safedata server API.

The top level command line help for the tool, showing the available subcommands
is shown below:

```bash
{!docs/command_line_tools/command_line_usage/safedata_server_top.txt!}
```

## The `safedata_server` subcommands

The command line help for each of the various subcommands is shown below:

### The `post_metadata` subcommand

```sh
{!docs/command_line_tools/command_line_usage/safedata_server_post_metadata.txt!}
```

### The `update_gazetteer` subcommand

```sh
{!docs/command_line_tools/command_line_usage/safedata_server_update_gazetteer.txt!}
```

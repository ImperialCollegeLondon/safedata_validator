# Command line tools overview

The `safedata_validator` package provides three main command line tools, but see also
the installation documentation for tools to build local copies of the
[GBIF](../install/build_local_gbif.md) and [NCBI](../install/build_local_ncbi.md)
taxonomic databases.

## `safedata_validate`

This tool is simply used to validate a dataset that is in [the `safedata`
format](../../data_providers/data_format/overview.md). The tool needs some simple
configuration, but is intended to be fairly easy to set up and use.

For more details, see [here](safedata_validate.md).

## `safedata_zenodo`

These tools are used to manage the process of publishing a validated dataset to Zenodo.
The configuration is a little more complex, as it requires the setup of a Zenodo
community account and access tokens.

For more details, see [here](safedata_zenodo.md).

## `safedata_server`

These tools are used to manage the updating an instance of the safedata metadata server,
including adding the details of published datasets and configuring core data such as the
gazetteer. The configuration requires a URL for a server running the safedata server API
and an access token for posting data to the server.

For more details, see [here](safedata_server.md).

# Command line tools overview

The `safedata_validator` package provides three main command line tools, but see also
the installation documentation for tools to build local copies of the
[GBIF](../install/build_local_gbif.md) and [NCBI](../install/build_local_ncbi.md)
taxonomic databases.

## Resources

All of the `safedata_validator` command line tools require access to a shared
configuration file providing validation resources and configuration information. These
files are described in the [configuration](../install/configuration.md) page and
normally all of the commands on this page will automatically load your configuration
data from the locations given on that page.

However, all three of the commands below also accept a manual path to a resource
configuration file using `--resources /path/to/file.cfg` or `-r /path/to/file.cfg`. This
can be useful, for example, if you want to have a separate configuration set up to use
the Zenodo sandbox for checking your publication process. If you want to check the
contents of a resources file, then all three commands also accept the `--show-resources`
or `-s` option, which will validate the resources and then print a summary to screen.

## `safedata_validate`

This tool is simply used to validate a dataset that is in [the `safedata`
format](../../data_providers/data_format/overview.md). The tool needs some simple
configuration, but is intended to be fairly easy to set up and use.

For more details, see [here](safedata_validate.md).

## `safedata_zenodo`

These tools are used to manage the process of publishing a validated dataset to Zenodo.
The configuration is a little more complex, as it requires the setup of a Zenodo
community account and access tokens and you may also want to edit the default
[HTML description](../install/configuration.md#html-description-template) used on Zenodo
and provide additional information to include XML summary metadata.

For more details, see [here](safedata_zenodo.md).

## `safedata_metadata`

These tools are used to manage the updating an instance of the safedata metadata server,
including adding the details of published datasets and configuring core data such as the
gazetteer. The configuration requires a URL for a server running the safedata server API
and an access token for posting data to the server.

For more details, see [here](safedata_metadata.md).

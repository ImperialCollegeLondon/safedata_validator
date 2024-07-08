# Command line tools overview

The `safedata_validator` package provides three main command line tools, but see also
the installation documentation for tools to build local copies of the
[GBIF](../install/build_local_gbif.md) and [NCBI](../install/build_local_ncbi.md)
taxonomic databases.

!!! tip

    These pages provides documentation on the details of the various command line tools
    but for guidance on _using the tools_, see these examples of:

    * [using the validator](../using_safedata/validating_datasets.md),
    * [publishing validated datasets](../using_safedata/publishing_datasets.md), and
    * [updating the metadata server](../using_safedata/posting_metadata.md).

## Resources

All of the `safedata_validator` command line tools require access to a shared
configuration file providing validation resources and configuration information. These
files are described in the [configuration](../install/configuration.md) page and
normally all of the commands on this page will automatically load your configuration
data from the locations given on that page.

However, it can sometimes be useful to use a different configuration. Example use cases
might be if you want to have a separate configuration set up to use the Zenodo sandbox
for checking your publication process, or if you want to add datasets to more than one
project. To support this, all three of the commands below also accept a manual path to a
resource configuration file using `--resources /path/to/file.cfg` or `-r
/path/to/file.cfg`. Alternatively, the tools all check if the environment variable
`SAFEDATA_VALIDATOR_CONFIG` has been set - this avoids having to use the `-r` flag with
every command.

If you want to verify the contents of a resources file, then you can obviously just open
and read the file, but all three commands also accept the `--show-resources` or `-s`
option, which will validate the resources and then print a summary to screen.

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

# Command line tools overview

The `safedata_validator` package provides two command line tools:

## `safedata_validate`

This tool is simply used to validate a dataset that is in the
[`safedata_validator` format](../data_providers/data_format/overview.md). The tool needs
some simple configuration, but is intended to be fairly easy to set up and use.

For more details, see [here](safedata_validate.md).

## `safedata_zenodo`

This tools is used to manage the process of publishing a validated dataset to Zenodo and
posting the metadata to a metadata server. The configuration is a little more complex,
as it requires the setup of a Zenodo community account and access tokens.

For more details, see [here](safedata_zenodo.md).

# Updating the metadata server

## Posting data for new published datasets

If you are maintaining a `safedata` metadata server then, after a dataset has been
published, you should post the metadata to the server to make it discoverable using the
`safedata` package. You will need to have configured `safedata_validator` to include the
required [metadata server
information](../install/configuration.md#metadata-configuration) but then the
[`safedata_metadata` command](../command_line_tools/safedata_metadata.md)
line usage is:

```sh
# Post the zenodo and dataset metadata to the metadata server
safedata_server post_metadata zenodo_1143714.json Example.json
```

## Updating the server resources

The metadata server is also used to distribute some of the validation resources used
both `safedata_validator` and `safedata`, including the
[gazetteer](../install/gazetteer.md) and the [Project ID
database](../install/configuration.md#validation-configuration). When you update these
resources for validation locally, those new data also need to be posted to the metadata
server. There are no command line arguments for this command - it simply updates the
server to the existing configuration file.

```sh
safedata_server update_resources
```

# Using the `safedata` system

<!-- markdownlint-disable MD046 -->

There are three stages for data managers in using the `safedata` system, all of which
use functionality provided by the `safedata_validator` package.

1. **Validation of datasets** using the `safedata_validate` tool. It is perfectly
   possible for data providers themselves to install and configure `safedata_validator`
   and this may be useful if you have a lot of data coming from few sources. However, it
   is more usual for the data manager to receive datasets and then go through a few
   cycles of validation and checking before a dataset is ready.

   See [here](./validating_datasets.md) for details.

2. **Publication of datasets** using the `safedata_zenodo` tool. Once a dataset has been
   validated, then it is ready for publication to the project Zenodo community. This
   involves creating a new Zenodo record, uploading the data files themselves and then
   filling in the required Zenodo metadata.

   See [here](./publishing_datasets.md) for details.

3. **Uploading dataset metadata** using the `safedata_server` tool. If you are
   maintaining a metadata server to support data discovery and use of the `safedata` R
   package, then the metadata for the published dataset needs to be uploaded to the
   server.

    See [here](./posting_metadata.md) for details.

The pages linked above provide examples of using the `safedata_validator` package.
Typically, data managers will use the command line interface using a Unix-like shell
or [Windows subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install),
but the pages also show how to use the programmatic API for `safedata_validator` from
within Python.

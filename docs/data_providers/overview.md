# Introduction

The `safedata_validator` package is built around making sure that datasets can
be:

* validated to a common standard,
* archived in an easily accessible location (currently [Zenodo](https://zenodo.org)),
  and
* have an easily queried metadata database.

Together, these aims help ensure that the datasets are safely stored and can be easily
found and reused by future researchers. However that does mean that the researchers
collecting those datasets provide some (relatively!) simple metadata.

If you are planning to submit data which needs to pass `safedata_validator`, the
process in overview is:

* Prepare your data following the [formatting information](data_format/overview.md).
* Submit the dataset to the data manager for the organisation you are submitting
  to, who will run the `safedata_validator` package on your dataset and provide
  automated feedback to help you address any issues:

    * If the validation succeeds then the dataset will be published on Zenodo. The
    dataset will be associated with the Zenodo community of the organisation you are
    submitting to.
    * If the validation fails then you will get an error report so you can fix the
    problems and resubmit.

* Once the dataset is published, the metadata for your dataset will be added to the
  organisation metadata server, so that it can be searched to help data discovery.

!!! Warning
    You **must not** publish your dataset directly to Zenodo. This skips validation
    and means that the datasets are not linked together under a single [curation
    account](availability.md#data-curation).

!!! Note "`safedata` at the SAFE Project"
    You can see published datasets for the [SAFE Project Zenodo
    community](https://zenodo.org/communities/safe/).

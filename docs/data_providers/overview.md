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
* Submit the dataset for validation. For the SAFE project, this is via the [SAFE Project
  website](https://safeproject.net/): you will need to be registered on the site and be
  a member of the project team that you are submitting data for.
  * If the validation succeeds then the dataset will be published on Zenodo
  * If the validation fails then you will get an error report so you can fix the
    problems and resubmit.

!!! Warning
    You **must not** publish your dataset directly to Zenodo. This skips validation
    and means that the datasets are not linked together under a single [curation account](availability.md#data-curation).

You can see published datasets at the [SAFE Project Zenodo
community](https://zenodo.org/communities/safe/).

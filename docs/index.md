# Overview

The `safedata_validator` package is a dataset validation and publishing tool for use
with large collections of related but heterogenous datasets. It was originally built to
handle the field data collections from the [SAFE Project](https://safeproject.net),
which has generated a [large number of related
datasets](https://zenodo.org/communities/safe).

It forms one element of the wider `safedata` system, which provide a general framework
for use in providing data capture, validation, publication and discovery within long
term ecological projects.

## User starting points

[Data providers](data_providers/overview.md)
: These pages provide an overview of the data preparation and submission process.

[Data managers](data_managers/command_line_tools/overview.md)
: These pages provide an introduction to the command line tools used to validate and
publish datasets.

[Developers](developers/api/overview.md)
: These pages provide the API for the package classes and methods.

## Package components

The `safedata_validator` package uses the following elements:

- A [defined format](data_providers/data_format/overview.md) for annotating data
  tables with metadata, validating sampling locations and taxonomic references and
  providing a standard set of summary metadata.
- The `safedata_validator` package itself, written in Python, which is used to check
  that a dataset meets the standard format. The package provides:
    - A programmatic API: `safedata_validator` can be imported and then used to
      implement validation within another program.
    - An implementation of the format checking process for data stored in Excel files,
      although implementations for other formats are possible.
    - Command line tools for both the validation (`safedata_validate`) and publication
      (`safedata_zenodo`) of datasets
- A data store to archive validated datasets. The `safedata_validator` packages
  currently uses [Zenodo](https://zenodo.org) to store datasets.

The wider `safedata` system then also uses:

- The [`safedata_server`](https://github.com/ImperialCollegeLondon/safedata_server) web
  application, which implements a metadata server for published datasets, used to
  maintain a searchable interface to the full metadata extracted from datasets.
- The [`safedata`](https://imperialcollegelondon.github.io/safedata/) package, which
  provides data discovery and download for published datasets into the [R statistical
  computing framework](https://www.r-project.org/).

## Code availability

The `safedata_validator` package is open source Python and  is maintained on
[GitHub](https://github.com/ImperialCollegeLondon/safedata_validator). It can
be installed using  [PyPI](https://pypi.org/project/safedata-validator).
See the [installation notes](data_managers/install/install.md) for setup instructions.

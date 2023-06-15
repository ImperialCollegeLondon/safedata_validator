# Overview

The `safedata_validator` package is part of a system for validating and publishing
heterogeneous datasets (the `safedata` system). It was originally built to handle the
field data collections from the [SAFE Project](https://safeproject.net), which generates
large numbers of related datasets. However, it is a general tool suitable for use within
long term ecological projects.

## User starting points

[Data providers](data_providers/overview.md)
: These pages provide an overview of the data preparation and submission process.

[Data managers](command_line_tools/overview.md)
: These pages provide an introduction to the command line tools used to validate and
publish datasets.

[Developers](api/overview.md)
: These pages provide the API for the package classes and methods.

## Package components

The `safedata_validator` package uses the following elements:

- A defined format for annotating data tables with metadata, validating sampling
  locations and taxonomic references and a standard set of summary metadata.
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
- A separate metadata server, which is used to maintain a searchable interface to the
  full metadata extracted from datasets.

## Code availability

The `safedata_validator` package is open source Python and  is maintained on
[GitHub](https://github.com/ImperialCollegeLondon/safedata_validator). It can
be installed using  [PyPI](https://pypi.org/project/safedata-validator).
See the [installation notes](install/install.md) for setup instructions.

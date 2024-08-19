# Something strange has gone wrong

This document exists to collect possible strange failures in the package testing or
documentation. These are issues we were aware of when implementing the relevant features
but couldn't think of a good way to fix. They are collected here for future reference.

## `test_example_dataset` suddenly fails

If an embargo date is provided, validation can only pass if this date is not in the past
and is not more than two years into the future (at least for the default SAFE Project
setup). `test_example_dataset` validates the `Example.xlsx` file which contains an
embargo date. This means that once this date is in the past the test will suddenly fail
without any modifications being made to any file.

## Docs build fails with message "Docs build process has resulted in file changes."

As part of the documentation build process the [`update_command_line_usage.sh`
script](../data_managers/command_line_tools/command_line_usage/update_command_line_usage.sh)
is run. This generates updated versions of the command line usage `.txt` files. The docs
build process then runs a check that the committed versions of the `.txt` files match
with the updated versions. If files have been altered by the docs build process the job
fails in order to indicate that the updated versions should be committed.

## Tables in data format summary documentation suddenly incorrect

The tables in the [data format documentation
pages](../data_providers/data_format/overview.md) are all created based on the contents
of `Example.xlsx` file. Most of the tables are generated from a specific sheet from the
file. However, the small tables in the [Summary documentation
page](../data_providers/data_format/summary.md) are generated based on particular lines
in the Summary sheet of the `Example.xlsx` file. The lines used to generate these tables
are hardcoded, meaning that if the order of fields in the `Example.xlsx` is changed the
table contents will become incorrect.

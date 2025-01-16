"""This file contains mkdocs hooks.

Details of how these hooks work can be found
[here](https://www.mkdocs.org/user-guide/configuration/#hooks). At present we use them
to automatically generate the files used to document the command line usage of
`safedata_validator`, and to generate a `.csv` file containing the information from the
Summary sheet of the Example.xlsx but with correctly formatted dates.

The `on_post_build` hook is then used to delete the `.csv` file once it has been used in
the docs build process.
"""

import os
import subprocess
from contextlib import contextmanager

import pandas as pd


@contextmanager
def work_in_dir(path):
    """Context manager for changing directories."""
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)


def on_pre_build(config, **kwargs) -> None:
    """Generate the files needed to fully render the documentation.

    This function calls the `update_command_line_usage.sh` so that the command line
    usage text files are regenerated every time the documentation is built.

    After that it loads in the summary sheet of the Example Excel file, reformats the
    dates to `DD-MM-YYYY` format, and outputs the reformatted data as a separate `.csv`
    file.
    """

    # Run the script to generate the command line usage text files
    with work_in_dir("docs/data_managers/command_line_tools/command_line_usage/"):
        subprocess.run(
            ["bash", "update_command_line_usage.sh"], stdout=subprocess.DEVNULL
        )

    # Generate tables for docs from the example Excel file.
    with work_in_dir("docs/data_providers/data_format/"):
        # Read in summary sheet as a pandas dataframe
        df = pd.read_excel(
            "Example.xlsx", sheet_name="Summary", header=None, keep_default_na=False
        )

        # Find rows that define dates and are not empty. Empty condition needed because
        # embargo date doesn't have to be provided.
        date_rows = df[
            (df[0].str.contains("date", case=False))
            & (~df[1].str.fullmatch("", na=False))
        ].index.tolist()

        # Parse dates as explicit dates using pandas and convert to DD-MM-YYYY format
        dates = pd.to_datetime(df.loc[date_rows, 1]).dt.strftime("%d-%m-%Y")

        # Replace the "dates" in the original dataframe with properly formatted dates
        df.loc[date_rows, 1] = dates

        # Output summary table as a csv
        df.to_csv("Summary.csv", index=False, header=False)


def on_post_build(config, **kwargs) -> None:
    """Delete temporary csv file that pre_build hook output."""

    # Delete the temporary csv file
    with work_in_dir("docs/data_providers/data_format/"):
        os.remove("Summary.csv")

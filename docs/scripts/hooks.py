"""This file contains mkdocs hooks.

Details of how these hooks work can be found
[here](https://www.mkdocs.org/user-guide/configuration/#hooks). At present we use them
to automatically generate the files used to document the command line usage of
`safedata_validator`, and to TODO - INSERT USE CASE 2.
"""

import os
import subprocess


def on_pre_build(config, **kwargs) -> None:
    """Generate the files needed to fully render the documentation.

    This function calls the `update_command_line_usage.sh` so that the command line
    usage text files are regenerated every time the documentation is built. It then
    TODO - Explain second thing that it does
    """

    # Change directory to command line usage folder
    os.chdir("docs/data_managers/command_line_tools/command_line_usage/")
    # Then run script to generate the command line usage text files
    subprocess.run(["sh", "update_command_line_usage.sh"], stdout=subprocess.DEVNULL)

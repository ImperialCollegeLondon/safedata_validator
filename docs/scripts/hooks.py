"""This file contains mkdocs hooks.

Details of how these hooks work can be found
[here](https://www.mkdocs.org/user-guide/configuration/#hooks). At present we use them
to INSERT USE CASE 1, INSERT USE CASE 2.
"""


def on_pre_build(config, **kwargs) -> None:
    """Generate the files needed to fully render the documentation."""
    # THINK I SHOULD LOOK INTO
    print("Hello world.")

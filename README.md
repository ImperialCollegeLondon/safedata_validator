# The safe_dataset_checker Python module

We're introducing a new format for tabular datasets collected at the SAFE project, which
involves including some fairly simple metadata in submitted files. This repository contains
Python code used to double check submitted files and report on any problems.

We currently only support datasets submitted as Excel workbooks: the vast majority of data
is submitted as Excel files (or is in some other spreadsheet format that could be saved
as Excel). We will work on other kinds of data - typically media files - but these often
are very large and will need long term bulk storage.

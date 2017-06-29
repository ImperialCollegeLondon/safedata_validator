# The safe_dataset_checker Python module

We're introducing a new format for tabular datasets collected at the SAFE project, which
involves including some fairly simple metadata in submitted files. This repository contains
Python code used to double check submitted files and report on any problems.

This code only supports datasets submitted as Excel workbooks: the vast majority of data
is submitted as Excel files (or is in some other spreadsheet format that could be saved
as Excel). We will work on other kinds of data - typically media files - but these often
are very large and will need long term bulk storage.

## Submitting datasets

Datasets can be submitted by registered researchers at the SAFE website, (https://safeproject.net/datasets/submit_dataset)
which will automatically use this code to check that the file is formatted correctly.

## Local use

The following steps should allow you to run the checks yourself before submitting a dataset:

1. Download the `safe_file_checker.py` file: 

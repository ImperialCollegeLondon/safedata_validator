# The safe_dataset_checker Python module

We're introducing a new format for tabular datasets collected at the SAFE project, which
involves including some fairly simple metadata in submitted files. This repository contains
Python code used to double check submitted files and report on any problems.

This code only supports datasets submitted as Excel workbooks: the vast majority of data
is submitted as Excel files (or is in some other spreadsheet format that could be saved
as Excel). We will work on other kinds of data - typically media files - but these often
are very large and will need long term bulk storage.

## Submitting datasets

Datasets can be submitted by registered researchers at the [SAFE website](https://safeproject.net/datasets/submit_dataset)
which will automatically use this code to check that the file is formatted correctly.

## Local use

The following steps should allow you to run the checks yourself before submitting a dataset:

1. Download the `safe_file_checker.py` file from [here](https://raw.githubusercontent.com/ImperialCollegeLondon/safe_dataset_checker/master/safe_dataset_checker.py) into the same folder as the Excel file containing your dataset.
2. Open a command line terminal and move to the directory where you saved `safe_file_checker.py`.
3. You should now be able to run the command below - a report on the file formatting will be outputted to the screen.

        python safe_file_checker.py MyDataset.xlsx
4. If this doesn't work, the most likely cause is that you need to install the `openpyxl` Python package, probably by using the `pip` package manager:

        pip install openpyxl

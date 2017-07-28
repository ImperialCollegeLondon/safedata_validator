# The `safe_dataset_checker` Python module

We're introducing a new format for tabular datasets collected at the
SAFE project, which involves including some fairly simple metadata in
submitted files. A description of the format can be found here:

https://www.safeproject.net/dokuwiki/working_at_safe/data_submission_format

This repository contains a Python module to validate submitted files and
report on any problems. We currently only support datasets submitted as
Excel workbooks: the vast majority of data is submitted as Excel files
(or is in some other spreadsheet format that could be saved as Excel).
We will work on other kinds of data - typically media files - but these
often are very large and will need long term bulk storage.

The code validates:

  1. The data submission formatting of the file.
  1. All taxonomic names against the NCBI taxonomy database.
  1. All location names against the SAFE Gazetteer.

Datasets can be submitted by registered researchers at the 
[SAFE website](https://safeproject.net/datasets/submit_dataset)
which will automatically use this code to check that the file is formatted correctly.
However, you may also want to run it yourself!

## Installation

The following steps should allow you to run the checks yourself before submitting a dataset:

1. Download the `safe_file_checker.py` file from [here](https://raw.githubusercontent.com/ImperialCollegeLondon/safe_dataset_checker/master/safe_dataset_checker.py) into the same folder as the Excel file containing your dataset.
2. Open a command line terminal and move to the directory where you saved `safe_file_checker.py`.
3. For guidance and testing, run the following:

       ./safe_dataset_checker.py -h

## Usage

The line above should show the message below:

    usage: safe_dataset_checker.py [-h] [-l LOCATIONS_JSON] [--use_entrez]
                                   [--check_all_ranks] fname

    This program validates an Excel file formatted as a SAFE dataset. As it runs,
    it outputs a report that highlights any problems with the formatting. The
    program validates taxonomic names against the NCBI taxonomy database. It will
    use the ete2 python package to use a local version if possible, but will
    otherwise attempt to use the Entrez web service to validate names. The program
    also validate sampling location names: by default, this is loaded
    automatically from the SAFE website so requires an internet connection, but a
    local copy can be provided for offline use.

    positional arguments:
      fname                 Path to the Excel file to be validated.

    optional arguments:
      -h, --help            show this help message and exit
      -l LOCATIONS_JSON, --locations_json LOCATIONS_JSON
                            Path to a locally stored json file of valid location
                            names
      --use_entrez          Use entrez queries for taxon validation, even if ete2
                            is available.
      --check_all_ranks     Check the validity of all taxonomic ranks included,
                            not just the standard required ranks.



In most cases, it is used simply like this:

    ./safe_dataset_checker.py path/to/My_Excel_File.xlsx

Note that the program __requires a web connection__ both to get a list 
of valid location names for SAFE and to validate taxon names via the
Entrez interface to the NCBI Taxonomy database:

https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

You can set the program up for offline use.

1. You will need to have installed and set up `ete2` (see below). Once 
installed, the program will use `ete2` unless you explicitly tell it not
to by adding the `--use_entrez` argument. You might want to do this
if you think your local copy of the NCBI database is out of date.

2. The current list of valid location names is automatically downloaded 
directly from the SAFE Gazetteer. To work offline, get a copy of
this list from the following link:

    https://www.safeproject.net/call/json/get_locations

    Save the output as a file (e.g. `SAFE_locations.json`). You will then
be able to run the program using the following:

        ./safe_dataset_checker.py My_Excel_File.xlsx --location_json SAFE_locations.json

## Requirements

You will need to install a couple of non-standard python packages, which ought
to be simple using the `pip` package manager.

  * `openpyxl`, which is needed to read the Excel files.

        pip install openpyxl

  * Optionally, `ete2`. Without this package, taxonomic names are
  validated using the NCBI Entrez web service. As an alternative, `ete2`
  provides an interface to a locally installed copy of the NCBI
  database. You'll need a good internet connection to install and set up
  `ete2` as it downloads about 40 MB from NCBI. This then needs
  processing and about 300 MB of disk space to build the local database
  and this can take some time (many minutes): only go down this route
  if you really need offline validation!

        pip install ete2
        # setup the local NCBI database
        python -c 'from ete2 import NCBITaxa;ncbi=NCBITaxa()'

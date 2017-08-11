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

    usage: safe_dataset_checker.py [-h] [-l LOCATIONS_JSON]
                                   [--ete3_database ETE3_DATABASE]
                                   [--check_all_ranks]
                                   fname

    This program validates an Excel file formatted as a SAFE dataset. As it runs,
    it outputs a report that highlights any problems with the formatting. The
    program validates taxonomic names against the NCBI taxonomy database. By
    default, it uses the Entrez web service to validate names, but if the ete3
    package is installed and the path to a built ete3 database is provided, then
    this will be used instead: this will work offline and is much faster but
    requires installation and setup. The program also validate sampling location
    names: by default, this is loaded automatically from the SAFE website so
    requires an internet connection, but a local copy can be provided for offline
    use.

    positional arguments:
      fname                 Path to the Excel file to be validated.

    optional arguments:
      -h, --help            show this help message and exit
      -l LOCATIONS_JSON, --locations_json LOCATIONS_JSON
                            Path to a locally stored json file of valid location
                            names
      --ete3_database ETE3_DATABASE
                            The path to a local NCBI Taxonomy database built for
                            the ete3 package.
      --check_all_ranks     Check the validity of all taxonomic ranks included,
                            not just the standard required ranks.




In most cases, it is used simply like this:

    ./safe_dataset_checker.py path/to/My_Excel_File.xlsx

Note that the program __requires a web connection__ both to get a list 
of valid location names for SAFE and to validate taxon names via the
Entrez interface to the NCBI Taxonomy database:

https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

You can set the program up for offline use.

1. You will need to have installed and set up `ete3` (see below) and
 then provide the path to the `ete3` local taxonomy database.

2. The current list of valid location names and bounding boxes is
automatically downloaded directly from the SAFE Gazetteer. To work
offline, get a copy of this data from the following link:

    https://www.safeproject.net/call/json/get_locations_bbox

    Save the output as a file (e.g. `SAFE_locations.json`). You will then
be able to run the program using the following:

        ./safe_dataset_checker.py My_Excel_File.xlsx --location_json SAFE_locations.json

## Requirements

You will need to install `openpyxl` to allow Python to read Excel files.
This ought to be simple using the `pip` package manager.

        pip install openpyxl

Optionally, you can use the `ete3` package to provide offline taxonomy
validation. This is not a simple option, so only go down this route
if you really need offline validation!

`ete3` is a great package but the way that the taxon checking code is
setup means that it is difficult to check the validity of the database
without triggering a complete rebuild. This is a bad thing to happen
if you are not expecting it since it downloads ~40MB of data from NCBI
and then spends several minutes building that into a local database
of around 300MB. As a result, `safe_dataset_checker` is very cautious
before it attempts to use it!

You will need a very recent version of `ete3` that provides a method
of checking this database without risking a complete rebuild. As of
August 2017, that means using the current `git` master:

    curl -O https://github.com/etetoolkit/ete/archive/master.zip
    unzip ete-master.zip
    cd ete-master
    python setup.py install

You will then need to build the database yourself. This isn't hard, but
does require a good internet connection and some time for the
processing. Start `python` and then run the following:

```python
import ete3
db = ete3.NCBITaxa(dbfile='path/to/loc/for/taxa.sqlite')
quit()
```

If you run `db = ete3.NCBITaxa()`, then `ete3` will store the database
by default in a hidden directory in your home directory
(`~/.etetoolkit/taxa.sqlite`). You can use this or specify your own
path as above. Either way, you will need to explicitly provide
`safe_dataset_checker` with  the location of this file before it will
use `ete3` for taxon validation.



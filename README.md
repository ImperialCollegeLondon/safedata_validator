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

## Usage

Most users will use the package directly as a command line program. For
guidance, run the following:

    ./safe_dataset_checker.py -h

In most cases, it is used simply like this:

    ./safe_dataset_checker.py My_Excel_File.xlsx

This requires a web connection to get the list of valid location names
and to validate taxon names via the Entrez interface to the NCBI
Taxonomy database:

https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

If you have installed and set up `ete2` (see below) then the same
command will use the local database unless you explicitly tell it not
to by adding the `--force_entrez` argument. You might want to do this
if you think your local copy of the NCBI database is out of date.

The list of valid location names is automatically downloaded directly
from the SAFE website. If you want to work offline, then get a copy of
the list from this link:

https://www.safeproject.net/call/json/get_locations

Save the output as a file (e.g. `SAFE_locations.json`). You will then
be able to run the program using the following:

    ./safe_dataset_checker.py My_Excel_File.xlsx --location_json SAFE_locations.json


## Requirements

If you want to use this code yourself, then will need to install a
couple of non-standard python packages, which ought to be simple using
the `pip` package manager.

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
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
  1. All taxonomic names against the GBIF taxonomy database.
  1. All location names against the SAFE Gazetteer.

Datasets can be submitted by registered researchers at the 
[SAFE website](https://safeproject.net/datasets/submit_dataset)
which will automatically use this code to check that the file is formatted
correctly. However, you may also want to run it yourself!

## Installation

The following steps should allow you to run the checks yourself before
submitting a dataset:

1. If you're using a Mac or Linux, then it is almost certain that you
already have Python installed. You may also do if you are using Windows.
To find out, open a command line window (run Terminal on a Mac, or `cmd`
on Windows) and then type `python` (on Mac or Unix) or `python.exe`
(on Windows).

1. If some text and a prompt (`>>>`) appears then you have Python.
First check the version number in the first line: if it doesn't start
Python 2.7 then you currently need to install Python 2.7 to run the
checker. You can have multiple versions of Python installed, but it is
going to be more complicated than is covered here.

  If you've got Python 2.7 then type `quit()` and skip to step 4.

  If you get a line that says the command is not found then you need to
install Python. Download a copy from here:

  https://www.python.org/downloads/

  The code is currently written to use Python 2.7, so make sure you
download an installer for Python 2.7.14 and not the more recent
Python 3.6 versions. For Windows, choose one of the MSI installer
options.

1. Repeat the command line check from the first step: if this still
doesn't work then you probably just need to tell the computer where to
find Python: search online for instructions to add python to the
`PATH` environment variable. On Windows, you will want to add (using
the typical install location) `C:\Python27` and `C:\Python27\scripts`.

1. The `safe_dataset_checker` program mostly uses commands from the
Python Standard Library - a set of code packages that are installed
with Python - but does use three extra packages that can be installed
using the `pip` package installer. At the command line, type:

        pip install openpyxl requests simplejson

  Those packages allow Python to: read Excel files, get validation data
over the internet and handle JSON formatted data.

1. Now create a folder to keep your data checking code in and download
the `safe_file_checker.py` file into it:

       https://raw.githubusercontent.com/ImperialCollegeLondon/safe_dataset_checker/master/safe_dataset_checker.py

1. Now open a command line terminal, change to your data checking
directory and run the following:

        python ./safe_dataset_checker.py -h

  In Windows, you will need to change `python` to `python.exe` here and
  in commands below that start `python`. You should see the usage
  instructions shown below.

## Usage

The program usage instructions are:

    usage: safe_dataset_checker.py [-h] [-l LOCATIONS_JSON]
                                   [--gbif_database GBIF_DATABASE]
                                   [--validate_doi]
                                   fname

    This program validates an Excel file formatted as a SAFE dataset. As it runs,
    it outputs a report that highlights any problems with the formatting. Much of
    the validation is to check that the data meets our metadata standards and is
    internally consistent. However, it uses external sources to perform validation
    in three areas.

    1. Taxon validation. The program validates taxonomic names against the GBIF
    taxonomy backbone. By default, it uses the GBIF web API to validate names,
    but can also use a local copy of the backbone provided in a sqlite database:
    this will work offline and is much faster but requires some simple setup.

    2. Location names. The program also validate sampling location names against
    the SAFE gazeteer. By default, this is loaded automatically from the SAFE
    website so requires an internet connection, but a local copy can be provided
    for offline use.

    3. DOI checking. Optionally, the program will validate any DOIs provided as
    having used the database. This requires a web connection and cannot be
    performed offline.

    positional arguments:
      fname                 Path to the Excel file to be validated.

    optional arguments:
      -h, --help            show this help message and exit
      -l LOCATIONS_JSON, --locations_json LOCATIONS_JSON
                            Path to a locally stored json file of valid location
                            names
      -g GBIF_DATABASE, --gbif_database GBIF_DATABASE
                            The path to a local sqlite database containing the
                            GBIF taxonomy backbone.
      --validate_doi        Check the validity of any publication DOIs, provided
                            by the user. Requires a web connection.


In most cases, it is used simply like this:

    python safe_dataset_checker.py path/to/My_Excel_File.xlsx

Note that the program __uses a web connection__  to get a list
of valid location names for SAFE, to validate taxon names via the
API to the [GBIF backbone taxonomy](https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c)
and to check any DOIs provided in the metadata. However, you can also
set the program up for offline use: this is also faster than using the
online services to validate the dataset. You cannot validate DOIs when
using the program offline.

## Offline usage

1. You will need to download a copy of the backbone taxonomy and build
a SQLite3 database from it. This isn't particularly hard, but the
resulting database is around 2GB, so you'll need file space! The steps
are:

    * You may need to download and install SQLite3.        You want the 'bundle' installer for your operating system.
      https://www.sqlite.org/download.html

    * Download the current zipped backbone from here:
      http://rs.gbif.org/datasets/backbone/backbone-current.zip

    * Extract the contents of the archive, which includes a tab
    delimited file `Taxon.tsv` that needs to be loaded into an SQLite3
    database. Move the `Taxon.tsv` file to your data checker folder and
    then, in the command line, type the following to create this
    new database and enter SQLite3:

           sqlite3 backbone-current.sqlite

    * Once you're in the new SQLite database, enter the following commands
    to import the data and then build two indices to speed up searches:
    on the canonical names of taxa and taxon ranks, and on the taxon id.
    This is a big file, so these might take several minutes to complete.

           .mode tab
           .import backbone-current/Taxon.tsv backbone
           # create covering indices
           create index backbone_name_rank on backbone (canonicalName, taxonRank);
           create index backbone_id on backbone (taxonID);
           .quit

    * You can now delete the `Taxon.tsv` file. You should now be able
    to use local taxonomy validation by telling `safe_dataset_checker.py`
    about the SQLite3 database:

           python safe_dataset_checker.py -g backbone-current.sqlite My_Excel_File.xlsx


2. The current list of valid location names and bounding boxes is
automatically downloaded directly from the SAFE Gazetteer. To work
offline, get a copy of this data from the following link:
https://www.safeproject.net/call/json/get_locations_bbox

    Save the output as a file in your data checker folder (e.g.
    `SAFE_locations.json`). You will then be able to run the program
    using the following:

          python safe_dataset_checker.py -l SAFE_locations.json My_Excel_File.xlsx

3. If you've done both the above steps then the following example would
validate a file using both local datasets, and won't need the internet
at all.

          python safe_dataset_checker.py -g backbone-current.sqlite -l SAFE_locations.json My_Excel_File.xlsx

4. You might want to create a shortcut to this command. On Windows,
create a file called `local_dataset_check.cmd` with the following
contents:

          python.exe safe_dataset_checker.py -g backbone-current.sqlite -l SAFE_locations.json %1

  You can now run the program using the local checking from a terminal by
typing:

          local_data_check.cmd My_Excel_File.xlsx

## GBIF validation

The basic idea is to use the taxon name and rank to find if there is a
match and the taxonomic status of the match. GBIF has the following
status codes:

    doubtful, accepted, homotypic synonym, synonym,
    heterotypic synonym, proparte synonym, misapplied

Misapplications and synonyms always have a suggested accepted usage,
doubtful taxa never do (at least when I checked a local backbone DB!).

### Online use
For online use, the program uses the following API to look for matches
against the GBIF backbone taxonomy. More details are provided at the
[GBIF developer site](https://www.gbif.org/developer/species).

http://api.gbif.org/v1/species/match?name=XX&strict=true&rank=YY

Validation against the GBIF species/match API returns a JSON dictionary
containing (along with the other results):
 - "matchType":"EXACT": an exact single match of the name XX at the
   taxonomic rank YY.
 - "matchType":"NONE": No exact single match found of that name at
    that rank.

Note that the matches can be made to any GBIF name field but with strict
matching this is almost certainly going to be to the canonical name
field (scientific name minus authorship). This can bring in multiple
matches, where a single canonical name is used in multiple synonyms: for
example _Zenicomus photuroides_ has one accepted use as _Zenicomus
photuroides_ Thomson, 1868 and about 10 synonyms with different
authors.

The GBIF API automatically looks for a single 'best' match - it will
return the single 'accepted' taxon with that canonical name over a set
of synonyms. The `verbose=true` argument can used to also return those
alternative matches but we don't use that information here as we don't
expect authorship details.

When GBIF returns a synonym, the program then looks up the accepted
usage. This is a little involved: the JSON data contains a usage key,
which refers to the name searched for and also a set of taxon keys,
which point to the accepted backbone taxon at each level.

You _can_ search for the accepted usage information using the usage key

http://api.gbif.org/v1/species/usagekey

The JSON contains an entry 'accepted' but this is the full scientific
name with authorship, and parsing those is a pain. It also contains
'acceptedKey', which you could feed back into the same API to get the
'canonicalName' for the accepted taxon. That is two API calls, so
instead the program uses the 'rank' entry from the match API to find
which accepted key to look up.

### Offline use

Validation against the local GBIF database works slightly differently
as it is a straight match against the canonical name with no preference
filtering. The equivalent database query is:

```{SQL}
select * from backbone
    where canonicalName = 'XX'
    and taxonRank = 'YY';
```

This will return all exact matches, including the less favoured
alternatives that GBIF omits, so that preference is reimplemented by the
program. The exact algorithm used by the API is unknown, so the program
just  looks for a single accepted use and then falls back to other
statuses. Note that currently the SQL query is case sensitive where the
API query is not.

Synonym handling is easier here because the database contains an
'acceptedNameUsageID' field that can be easily used to get the accepted
name. However, the higher taxonomy is less easy because, although higher
taxon names are provided, their taxon IDs aren't.

### Problems

Rare edge cases include taxon names with two equally approved usages:
for example, the genus _Morus_ is an accepted usage for both mulberries
and gannets. In these rare cases, an accepted usage would require a
parent taxon to discriminate between them.


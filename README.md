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


We want three things from validation:

* The status of the name provided by the user.

* The accepted usage (which might be the same).

* The backbone taxonomic hierarchy for the taxa, so we can index at
  higher taxonomic levels. However, this is complicated as there is
  no guarantee that taxa will always hook into the backbone at the next
  level above their own.

  For example, the accepted taxon _Wanosuchus atresus_ is only hooked
  in at order level: it is accepted as a species, but its parent is
  Crocodylia as the genus is doubtful. Similarly, _Goniopholis
  tenuidens_ is a synonym at species level but again has Crocodylia as
  a parent (and is considered a synonym for the family Goniopholidae).

  In some cases, taxa can hook in at levels below their own: the
  species _Brittonastrum greenei_ is a synonym of the **subspecies**
  _Agastache pallidiflora pallidiflora_.

  The `parent_taxon_level_analysis.R` file in this repository contains
  some code to check this:
  
  * All **accepted taxa** map to a more nested parent but 5%  map to a
  more nested parent more than one step up the hierarchy. The table
  below shows child taxon level as rows and parent taxon level as
  columns.

  |           | kingdom| phylum| class| order| family|   genus| species| subspecies| variety| form|
  |:----------|-------:|------:|-----:|-----:|------:|-------:|-------:|----------:|-------:|----:|
  |kingdom    |       0|      0|     0|     0|      0|       0|       0|          0|       0|    0|
  |phylum     |     100|      0|     0|     0|      0|       0|       0|          0|       0|    0|
  |class      |       5|    316|     0|     0|      0|       0|       0|          0|       0|    0|
  |order      |       7|     45|  1327|     0|      0|       0|       0|          0|       0|    0|
  |family     |    2191|   1339|  4267| 14423|      0|       0|       0|          0|       0|    0|
  |genus      |    3427|   4985|  5584|  6260| 220735|       0|       0|          0|       0|    0|
  |species    |    1567|    706|  1529|   696|   8944| 2449414|       0|          0|       0|    0|
  |subspecies |      41|      7|     3|     2|    832|     268|  200902|          0|       0|    0|
  |variety    |      53|     10|     0|    26|   2661|      50|   82914|         32|       0|    0|
  |form       |      12|      4|     0|     4|    815|      18|   19272|          0|      56|    0|

  * Only 77% of **unaccepted taxa** map to a parent at the next most
  nested taxonomic level and 4.5% map to a parent at the same or a less
  nested level, as in the example above.

  |           | kingdom| phylum| class| order| family|   genus| species| subspecies| variety| form|
  |:----------|-------:|------:|-----:|-----:|------:|-------:|-------:|----------:|-------:|----:|
  |kingdom    |       0|      0|     0|     0|      0|       0|       0|          0|       0|    0|
  |phylum     |      22|      8|     0|     0|      0|       0|       0|          0|       0|    0|
  |class      |       0|     14|     1|     0|      0|       0|       0|          0|       0|    0|
  |order      |       0|      5|    32|     0|      0|       0|       0|          0|       0|    0|
  |family     |      21|    157|   481|  3599|      0|       0|       0|          0|       0|    0|
  |genus      |    8555|  24242| 25055| 31010| 185911|       0|       0|          0|       0|    0|
  |species    |      64|     24|   173|   405|   2142| 1886329|  121225|         84|       5|    0|
  |subspecies |       3|      0|     1|     0|    151|   77512|   26266|         13|       0|    0|
  |variety    |       2|      1|     0|     2|    367|  212954|   50062|         47|       4|    0|
  |form       |       0|      0|     0|     0|    128|   48126|   10449|          3|       2|    0|

The two different APIs for offline and online usage work to try and get
the same answer.

### Online use

The GBIF API provides the `species/match` endpoint, which allows searching for
names against the GBIF backbone taxonomy. More details are provided at
the [GBIF developer site](https://www.gbif.org/developer/species), but
the endpoint is handy as it accepts variables to turn off fuzzy matching
(`strict=true`) and to check the name is at the stated rank (e.g.
`rank=species`). For example:

http://api.gbif.org/v1/species/match?name=Panthera%20leo&strict=true&rank=species

The JSON data returned looks like this:

```
{"usageKey":5219404,
 "scientificName":"Panthera leo (Linnaeus, 1758)",
 "canonicalName":"Panthera leo",
 "rank":"SPECIES",
 "status":"ACCEPTED",
 "confidence":99,
 "matchType":"EXACT",
 "kingdom":"Animalia",
 "phylum":"Chordata",
 "order":"Carnivora",
 "family":"Felidae",
 "genus":"Panthera",
 "species":"Panthera leo",
 "kingdomKey":1,
 "phylumKey":44,
 "classKey":359,
 "orderKey":732,
 "familyKey":9703,
 "genusKey":2435194,
 "speciesKey":5219404,
 "synonym":false,
 "class":"Mammalia"}
```

With strict searching, here we have a single exact single match
(`"matchType":"EXACT"`) of the name _Panthera leo_ at species level.
If no exact match is found at the rank then instead we would get
`"matchType":"NONE"` and possibly some `"notes"` with further details.
Note that the `match` endpoint checks names against multiple fields but
with strict matching this is almost certainly going to be to the canonical name
field (scientific name minus authorship).

A search can bring in multiple matches, where a single canonical name is
used in multiple synonyms: for example _Zenicomus photuroides_ has one
accepted use as _Zenicomus photuroides_ Thomson, 1868 and about 10
synonyms with different authors.

The GBIF `species/match` API automatically looks for a single 'best' match - it
will return the single taxon with that canonical name over a set of
synonyms: the JSON data will contain `"status": "ACCEPTED"`. The
`verbose=true` argument to the API can used to also return those
alternative matches but we don't use that information here as we don't
expect authorship details. If there are only non-accepted usages, then
the `status` will be something else: for example `SYNONYM` or
`MISAPPLIED`.

For accepted taxa, this information is enough to get the information
needed: we've got the accepted usage and a set of backbone taxa. The
only things to be wary of are potential missing levels and the fact
that taxon keys only go up to species level. In the example above the
`usageKey` matches the `speciesKey`, but the backbone includes
subspecies, variety and form levels. We can hook subspecies in, since
the `speciesKey` is the least nested parent possible but variety and
form do not have enough information to build the full chain without
further queries. These levels are not supported by the checker.

When GBIF returns a synonym, there is more to be done. The JSON data
from the `species/match` endpoint contains the `usageKey` for the search
name. That _probably_ means that the highest taxon key provided is the
accepted taxon and the next highest is the parent. However, the program
currently checks this using the `species/{usageKey}` endpoint, which
explictly contains an `acceptedKey` and a `parentKey` which ought to
match the previous information. For example:

http://api.gbif.org/v1/species/8548608


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

In many ways, this is easier to use because the database rows contain
both `acceptedNameUsageID` and `parentNameUsageID` fields that can be
easily used to get the accepted name and hook the taxon into the higher
taxonomy. Although canonical higher taxon names are provided, their
taxon IDs are not, so higher taxa are saved into a set containing unique
pairs of names and ranks and then added to the index when all taxa have
been processed.

### Problems

Rare edge cases include taxon names with two equally approved usages:
for example, the genus _Morus_ is an accepted usage for both mulberries
and gannets. This kind of problem is described in the JSON `"note"`
field and provided to the user. In these rare cases, an accepted
usage would require a parent taxon to discriminate between them.


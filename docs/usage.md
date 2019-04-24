# Using the SAFE Dataset Checker

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

Note that, by default, the program **uses a web connection** to:

1. get a list of valid location names for SAFE, 
2. validate taxon names via the API to the [GBIF backbone taxonomy](https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c), and
3. check any DOIs provided in the metadata. 

However, you can also set the program up for offline use: this is also much faster than using the online
services to validate the dataset but you cannot validate DOIs when using the program offline. To use the program offline, you need to provide paths to a local copy of the GBIF database and a local copy of the SAFE gazetteer.

    python safe_dataset_checker.py -g backbone-current.sqlite -l SAFE_locations.json My_Excel_File.xlsx

Details on obtaining those two files are provided [here]().

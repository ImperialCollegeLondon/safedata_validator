import argparse
import textwrap
from safedata_validator.version import __version__
from safedata_validator.field import Dataset

def _safedata_validator_cli():
    """
    This program validates an Excel file formatted as a SAFE dataset. As it runs, it outputs
    a report that highlights any problems with the formatting.
    Much of the validation is to check that the data meets our metadata standards and is
    internally consistent. However, it uses external sources to perform validation in three
    areas.
    1. Taxon validation. The program validates taxonomic names against the GBIF taxonomy
    backbone. By default, it uses the GBIF web API to validate names, but can also use a
    local copy of the backbone provided in a sqlite database: this will work offline and
    is much faster but requires some simple setup.
    2. Location names. The program also validate sampling location names against the SAFE gazeteer.
    By default, this is loaded automatically from the SAFE website so requires an internet
    connection, but a local copy can be provided for offline use.
    3. DOI checking. Optionally, the program will validate any DOIs provided as having used
    the database. This requires a web connection and cannot be performed offline.
    """

    desc = textwrap.dedent(_safedata_validator_cli.__doc__)
    fmt = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, formatter_class=fmt)

    parser.add_argument('filename', help="Path to the Excel file to be validated.")
    # parser.add_argument('-c', '--check', default='sltwf',
    #                     help='Which of the summary, locations, taxa, worksheets and '
    #                          'finalisation should be checked.')
    parser.add_argument('-p', '--project_id', default=None, type=int, action='append',
                        help='If provided, check that the project ID within the file '
                             'matches this integer. Multiple values can be provided '
                             'to generate a set of valid IDs.', dest='valid_pid')
    parser.add_argument('-l', '--locations', default=None,
                        help='A path to a locally stored locations data file or the URL '
                             'of a web service that provides the same data.')
    parser.add_argument('-g', '--gbif_database', default=None,
                        help=('The path to a local sqlite database containing the GBIF '
                              'taxonomy backbone.'))
    parser.add_argument('--validate_doi', action="store_true", default=False,
                        help=('Check the validity of any publication DOIs, '
                              'provided by the user. Requires a web connection.'))
    parser.add_argument('--chunk_size', default=1000, type=int, 
                        help=('Data are loaded from worksheets in chunks: the'
                              'number of rows in a chunk is set by this argument'))
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    
    args = parser.parse_args()
    
    
    ds = Dataset()
    ds.load_from_workbook(filename=args.filename, valid_pid=args.valid_pid,
                          validate_doi=args.validate_doi, chunk_size=args.chunk_size)

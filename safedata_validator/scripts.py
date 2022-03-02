import argparse
import textwrap
from safedata_validator.version import __version__
from safedata_validator.field import Dataset

def _safedata_validator_cli():
    """
    This program validates an Excel file formatted as a SAFE dataset. As it
    runs, it outputs a report that highlights any problems with the formatting.
    Much of the validation is to check that the data meets our metadata
    standards and is internally consistent. 
    
    However, the package uses external resources to perform validation of taxa
    and sampling locations and to provide other information. For this reason,
    using this program requires you to provide a configuration file for these
    resources or to have installed a configuration file in a standard location
    (see the package website API for details.)
    """

    desc = textwrap.dedent(_safedata_validator_cli.__doc__)
    fmt = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, formatter_class=fmt)

    parser.add_argument('filename', help="Path to the Excel file to be validated.")
    parser.add_argument('-p', '--project_id', default=None, type=int, action='append',
                        help='If provided, check that the project ID within the file '
                             'matches this integer. Multiple values can be provided '
                             'to generate a set of valid IDs.', dest='valid_pid')
    parser.add_argument('--validate_doi', action="store_true", default=False,
                        help=('Check the validity of any publication DOIs, '
                              'provided by the user. Requires a web connection.'))
    parser.add_argument('--chunk_size', default=1000, type=int, 
                        help=('Data are loaded from worksheets in chunks: the '
                              'number of rows in a chunk is set by this argument'))
    
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    
    args = parser.parse_args()
    
    ds = Dataset()
    ds.load_from_workbook(filename=args.filename, valid_pid=args.valid_pid,
                          validate_doi=args.validate_doi, chunk_size=args.chunk_size)

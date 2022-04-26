import os
import sys
import argparse
import textwrap
from safedata_validator.version import __version__
from safedata_validator.field import Dataset
from safedata_validator.zenodo import download_ris_data
from safedata_validator.logger import LOGGER

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

    If validation is successful, then a JSON format file containing key metadata
    will be saved to the same location as the valdiated file. The JSON metadata
    is used in the dataset publication process.
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

    if ds.passed:
        json_file =   os.path.splitext(args.filename)[0] + '.json'

        with open(json_file, 'w') as json_out:
            json_out.write(ds.to_json())
        
        sys.stdout.write('------------------------\n')
        sys.stdout.write(f'JSON metadata written to {json_file}\n')
        sys.stdout.write('------------------------\n')





def _safedata_download_ris_cli():
    """
    This program maintains a RIS format bibliography file of the datasets
    uploaded to a Zenodo community. It can update an existing RIS format file
    to add new records or it can create the file from scratch.

    The program uses both the Zenodo API (to find the records in the community)
    and the Datacite API to access machine readable biblioigraphic records.
    """

    desc = textwrap.dedent(_safedata_download_ris_cli.__doc__)
    fmt = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, formatter_class=fmt)
    
    # positional argument inputs
    parser.add_argument('ris_file', type=str,
                        help="The file path to populate with RIS records. If this file "
                             "already exists, it is assumed to be RIS file to update "
                             "with any new records not already included in the file.")

    args = parser.parse_args()

    # Run the download RIS data function
    data = download_ris_data(ris_file=args.ris_file)
    if data is None:
        LOGGER.info(f'No new records found')
    elif os.path.exists(args.ris_file):
        LOGGER.info(f'Appending RIS data for {len(data)} new records')
        write_mode = 'a'
    else:
        LOGGER.info(f'Downloading RIS data for {len(data)} records')
        write_mode = 'w'

    with open(args.ris_file, write_mode) as ris_file:
        for this_entry in data:
            ris_file.write(this_entry)
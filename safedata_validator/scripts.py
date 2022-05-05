import os
import sys
import argparse
import textwrap
import simplejson

from safedata_validator.version import __version__
from safedata_validator.field import Dataset
from safedata_validator.zenodo import (create_deposit, delete_file, discard_deposit, get_deposit, download_ris_data, 
                                       sync_local_dir, upload_file)
from safedata_validator.resources import Resources
from safedata_validator.logger import FORMATTER, LOGGER, CONSOLE_HANDLER

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


def _safedata_zenodo_cli():
    """
    This is a the command line interface for publishing safedata 
    validated datasets to Zenodo, downloading information and maintaining a
    local copy of the datasets in the file structure required by the R safedata
    package.
    
    See the individual help for each of the subcommands
    """

    # create the top-level parser and configure to take subparsers
    desc = textwrap.dedent(_safedata_zenodo_cli.__doc__)
    fmt = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, formatter_class=fmt)
    parser.add_argument('-d', '--debug', action='store_true', 
                        help="Show debugging for resource loading")
    parser.add_argument('-r', '--resources', type=str, default=None,
                        help="Path to a safedata_validator resource configuration file")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Surpress normal information messages. ")

    subparsers = parser.add_subparsers(dest='subcommand')

    # Setup the parser for the "create" command
    create_desc ="""Create a new deposit draft, possibly as a new version of an
    existing _published_ record. When successful, the function also downloads
    and saves a JSON file containing the Zenodo deposit metadata. This file is
    used as an input to other subcommands that work with an existing deposit.
    """
    parser_create = subparsers.add_parser('create', 
                                          description = textwrap.dedent(create_desc),
                                          help="Create a new deposit")
    
    parser_create.add_argument('-d', '--deposit_id', type=int, default=None, 
                         help='The draft created will be a new version of the '
                         'existing published Zenodo record with this ID')

    # Setup the parser for the "upload" command
    upload_desc = """Uploads the contents of a specified file to an
    _unpublished_ Zenodo deposit, optionally using an alternative filename. 
    """

    parser_upload = subparsers.add_parser('upload', 
                                          description = textwrap.dedent(upload_desc),
                                          help='Upload a file')
        
    parser_upload.add_argument('deposit_json', type=str, default=None, 
                               help='Path to a Zenodo metadata JSON for the deposit')
    parser_upload.add_argument('filepath', type=str, default=None, 
                                help='The path to the file to be uploaded')
    parser_upload.add_argument('--zenodo_filename', type=str, default=None, 
                                help='An optional alternative file name to be used on Zenodo')

    # INFO subcommand
    info_desc="""Print out summary information on a deposit and download the
    Zenodo metadata for that deposit."""

    parser_info = subparsers.add_parser('info', 
                                        description = textwrap.dedent(info_desc),
                                        help="Display and download deposit metadata")
    parser_info.add_argument('deposit_id', type=int, default=None, 
                             help='An ID for an existing Zenodo deposit')

    # Discard subcommand
    discard_desc="""Discard an unpublished deposit. The deposit and all uploaded
    files will be removed from Zenodo."""

    parser_discard = subparsers.add_parser('discard', 
                                        description = textwrap.dedent(discard_desc),
                                        help="Discard an unpublished deposit")
    
    parser_discard.add_argument('deposit_json', type=str, default=None, 
                               help='Path to a Zenodo metadata file for the deposit to discard')

    # Delete_file subcommand
    discard_desc="""Delete an uploaded file from an unpublished deposit. The
    deposit metadata will be re-downloaded to ensure an up to date list of files
    in the deposit.
    """

    parser_discard = subparsers.add_parser('delete_file', 
                                        description = textwrap.dedent(discard_desc),
                                        help="Delete a file from an unpublished deposit")
    
    parser_discard.add_argument('deposit_json', type=str, default=None, 
                               help='Path to a Zenodo metadata file for the deposit to discard')
    parser_discard.add_argument('filename', type=str, default=None, 
                               help='The name of the file to delete')

    # sync_local_dir subcommand
    
    sync_local_desc = """
    The safedata R package defines a directory structure used to store metadata
    and files downloaded from a safedata community on Zenodo and from a safedata
    metadata server. This tool allows a safedata developer or community
    maintainer to create or update such a directory with _all_ of the resources
    in the Zenodo community, regardless of their public access status. This
    forms a backup (although Zenodo is heavily backed up) but also provides
    local copies of the files for testing and development of the code packages.

    If this is a new data directory, you also need to provide the API url for
    the a server providing safedata metadata.
    """

    parser_sync = subparsers.add_parser('sync_local_dir', 
                                    description = textwrap.dedent(sync_local_desc),
                                    formatter_class=fmt,
                                    help="Create or update a local safedata directory")
    
    parser_sync.add_argument('datadir', type=str, help='The path to a local directory containing '
                        'an existing safedata directory or an empty folder in which to create one')
    parser_sync.add_argument('-a', '--api', type=str, default=None,
                        help='An API from which JSON dataset metadata can be downloaded. If '
                        'datadir is an existing safedata directory, then the API will be read '
                        'from `url.json`.')
    parser_sync.add_argument('--xlsx_only', type=bool, default=True,
                        help='Should the download ignore large non-xlsx files, defaulting '
                        'to True.')

    # ris subcommand

    ris_desc = """
    This command maintains a RIS format bibliography file of the datasets
    uploaded to a Zenodo community. It can update an existing RIS format file
    to add new records or it can create the file from scratch.

    The program uses both the Zenodo API (to find the records in the community)
    and the Datacite API to access machine readable biblioigraphic records.
    """

    parser_ris = subparsers.add_parser('ris', 
                                description = textwrap.dedent(ris_desc),
                                formatter_class=fmt,
                                help="Maintain a RIS bibliography file for datasets")
    
    # positional argument inputs
    parser_ris.add_argument('ris_file', type=str,
                        help="The file path to populate with RIS records. If this file "
                             "already exists, it is assumed to be RIS file to update "
                             "with any new records not already included in the file.")

    # Parse the arguments and set the verbosity
    args = parser.parse_args()

    if args.debug and args.quiet:
        args.quiet = False 
        
    # Load the package config 
    if not args.debug:
        CONSOLE_HANDLER.setLevel('CRITICAL')
    else:
        CONSOLE_HANDLER.setLevel('DEBUG')

    resources = Resources(args.resources)

    if args.quiet:
        # Don't suppress error messages
        CONSOLE_HANDLER.setLevel('ERROR')
    else:
        CONSOLE_HANDLER.setLevel('DEBUG')

    # Handle the different subcommands
    if args.subcommand == 'create':

        response, error = create_deposit(deposit_id=args.deposit_id,
                                         resources=resources)

        if error is not None:
            LOGGER.error(f'Failed to create deposit: {error}')
            return
        
        rec_id = response['record_id']
        LOGGER.info(f'Created deposit: {recid}')
        outfile = os.path.join(os.getcwd(), f"zenodo_{recid}.json")
        with open(outfile, 'w') as outf:
            simplejson.dump(links, outf)
            LOGGER.info(f'Metadata downloaded to: {outfile}')
        
        return
    
    if args.subcommand == 'upload':

        # Load the Zenodo deposit JSON, which contains API links
        with open(args.deposit_json) as dep_json:
            metadata = simplejson.load(dep_json)
        
        # Report on proposed upload
        if args.zenodo_filename is None:
            LOGGER.info(f'Uploading: {os.path.basename(args.filepath)}')
        else:
            LOGGER.info(f'Uploading: {os.path.basename(args.filepath)} as {args.zenodo_filename}')
        
        # Run the command
        response, error = upload_file(metadata=metadata,
                                      filepath=args.filepath, 
                                      zenodo_filename=args.zenodo_filename,
                                      resources=resources,
                                      progress_bar = not args.quiet)
        
        # Report on the ouctome.
        if error is not None:
            LOGGER.error(f'Failed to upload file: {error}')
        else:
            LOGGER.info('File uploaded')

        return

    if args.subcommand == 'info':

        # Run the command
        response, error = get_deposit(deposit_id=args.deposit_id,
                                       resources=resources)
    
        if error is not None:
            LOGGER.error(f'Failed to get info: {error}')
            return

        # Dump the response to a JSON file
        outfile = os.path.join(os.getcwd(), f"zenodo_{args.deposit_id}.json")
        with open(outfile, 'w') as outf:
            simplejson.dump(response, outf)

        # Print a short summary
        LOGGER.info(f"Record ID: {response['record_id']}")
        LOGGER.info(f"Concept ID: {response['conceptrecid']}")
        LOGGER.info(f"Status: {'published' if response['submitted'] else 'Not published'}")
        LOGGER.info(f"Title: {response['title'] if response['title'] else 'Not set'}")

        if len(response['files']) == 0:
            LOGGER.info(f"Files: none uploaded")
        else:        
            LOGGER.info(f"Files:")
            FORMATTER.push()
            for this_file in response['files']:
                LOGGER.info(f"{this_file['filename']} ({this_file['filesize']} bytes, "
                            f"md5: {this_file['checksum']})")
            FORMATTER.pop()

        LOGGER.info(f'Metadata downloaded to: {outfile}')

    if args.subcommand == 'discard':

        # Load the Zenodo deposit JSON, which contains API links
        with open(args.deposit_json) as dep_json:
            metadata = simplejson.load(dep_json)
        
        # Run the command
        response, error = discard_deposit(metadata=metadata,
                                          resources=resources)
        
        # Report on the ouctome.
        if error is not None:
            LOGGER.error(f'Failed to discard deposit: {error}')
        else:
            LOGGER.info('Deposit discarded')

        return

    if args.subcommand == 'delete_file':

        # Load the Zenodo deposit JSON, which contains API links
        with open(args.deposit_json) as dep_json:
            metadata = simplejson.load(dep_json)
        
        # Run the command
        response, error = delete_file(metadata=metadata,
                                      filename=args.filename,
                                      resources=resources)
        
        # Report on the ouctome.
        if error is not None:
            LOGGER.error(f'Failed to delete file: {error}')
        else:
            LOGGER.info('File deleted')

        return

    if args.subcommand == 'sync_local_dir':

        sync_local_dir(datadir=args.datadir,
                       api=args.api,
                       xlsx_only=args.xlsx_only,
                       resources=resources)

    if args.subcommand == 'ris':

        # Run the download RIS data function
        download_ris_data(ris_file=args.ris_file,
                          resources=resources)
        

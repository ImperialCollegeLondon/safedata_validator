"""Provide command line scripts.

This module provides command line scripts to expose validation and publishing
functionality.

* _safedata_validate_cli, exposed as `safedata_validate`
* _safedata_zenodo_cli, exposed as `safedata_zenodo`
"""

import argparse
import os
import sys
import tempfile
import textwrap
from pathlib import Path

import simplejson

from safedata_validator import __version__
from safedata_validator.field import Dataset
from safedata_validator.logger import (
    FORMATTER,
    LOGGER,
    get_handler,
    use_file_logging,
    use_stream_logging,
)
from safedata_validator.resources import Resources
from safedata_validator.taxondb import (
    build_local_gbif,
    build_local_ncbi,
    download_gbif_backbone,
    download_ncbi_taxonomy,
)
from safedata_validator.zenodo import (
    create_deposit,
    dataset_description,
    delete_file,
    discard_deposit,
    download_ris_data,
    generate_inspire_xml,
    get_deposit,
    post_metadata,
    publish_deposit,
    sync_local_dir,
    update_gazetteer,
    update_published_metadata,
    upload_file,
    upload_metadata,
)


def _safedata_validator_cli():
    """Validate a dataset using a command line interface.

    This program validates an Excel file formatted as a `safedata` dataset.
    As it runs, it outputs a report that highlights any problems with the
    formatting. Much of the validation is to check that the data meets our
    metadata standards and is internally consistent.

    However, the package uses external resources to perform validation of
    taxa and sampling locations and to provide other information. For
    this reason, using this program requires you to provide a configuration
    file for these resources or to have installed a configuration file in a
    standard location. If you run `safedata_validate` without a
    configuration file, the output will report the standard locations for
    your operating system.

    If validation is successful, then a JSON format file containing key
    metadata will be saved to the same location as the validated file.
    The JSON metadata is used in the dataset publication process.
    """

    desc = textwrap.dedent(_safedata_validator_cli.__doc__)
    fmt = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, formatter_class=fmt)

    parser.add_argument("filename", help="Path to the Excel file to be validated.")

    parser.add_argument(
        "-r",
        "--resources",
        type=str,
        default=None,
        help=("A path to a resources configuration file"),
    )
    parser.add_argument(
        "--validate_doi",
        action="store_true",
        default=False,
        help=(
            "Check the validity of any publication DOIs, "
            "provided by the user. Requires a web connection."
        ),
    )
    parser.add_argument(
        "--chunk_size",
        default=1000,
        type=int,
        help=(
            "Data are loaded from worksheets in chunks: the "
            "number of rows in a chunk is set by this argument"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        type=Path,
        help=("Save the validation report to a file, not print to the console."),
    )

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )

    args = parser.parse_args()

    # Configure the logging location
    if args.output is None:
        use_stream_logging()
    else:
        use_file_logging(args.output)

    # Create the dataset and load from workbook
    ds = Dataset(resources=Resources(args.resources))
    ds.load_from_workbook(
        filename=args.filename,
        validate_doi=args.validate_doi,
        chunk_size=args.chunk_size,
    )

    if ds.passed:
        json_file = os.path.splitext(args.filename)[0] + ".json"

        with open(json_file, "w") as json_out:
            json_out.write(ds.to_json())

        sys.stdout.write("------------------------\n")
        sys.stdout.write("File validation passed\n")
        sys.stdout.write(f"JSON metadata written to {json_file}\n")
        sys.stdout.write("------------------------\n")


def _desc_formatter(prog):
    """Bespoke argparse description formatting."""
    return argparse.RawDescriptionHelpFormatter(prog, max_help_position=16)


def _safedata_zenodo_cli():
    """Publish validated datasets to Zenodo using a command line interface.

    This is a the command line interface for publishing safedata validated
    datasets to Zenodo, downloading information and maintaining a local
    copy of the datasets in the file structure required by the R safedata
    package.

    The safedata_zenodo command is used by providing subcommands for the
    different actions required to publish a validated dataset. The list of
    subcommands (with aliases) is shown below and individual help is
    available for each of the subcommands:

        safedata_zenodo subcommand -h

    The subcommands for this tool use two different JSON format metadata
    files:

    * A dataset metadata file (`dataset_json`). This is the output from
        using the `safedata_validate` tool. Some of the information in this
        file is used to create the Zenodo dataset description, and all of
        the data is used to describe a dataset on the separate metadata
        server.

    * A Zenodo metadata file (`zenodo_json`), that describes the metadata
        associated with a Zenodo deposit or published record.

    Note that most of these actions are also available via the Zenodo website.
    """

    # create the top-level parser and configure to take subparsers
    desc = textwrap.dedent(_safedata_zenodo_cli.__doc__)
    fmt = _desc_formatter
    parser = argparse.ArgumentParser(
        prog="safedata_zenodo",
        description=desc,
        formatter_class=fmt,
        usage="safedata_zenodo [-h] [-r RESOURCES] [-q] SUBCOMMAND ...",
    )

    parser.add_argument(
        "-r",
        "--resources",
        type=str,
        default=None,
        help="Path to a safedata_validator resource configuration file",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress normal information messages. ",
    )

    subparsers = parser.add_subparsers(dest="subcommand", metavar="")

    # CREATE DEPOSIT subcommand
    create_deposit_desc = """
    Create a new deposit draft. The concept_id option uses a provided Zenodo
    concept ID to creates a draft as a new version of an existing data set.

    When successful, the function downloads and saves a JSON file containing the
    resulting Zenodo deposit metadata. This file is used as an input to other
    subcommands that work with an existing deposit.
    """
    create_deposit_parser = subparsers.add_parser(
        "create_deposit",
        description=textwrap.dedent(create_deposit_desc),
        help="Create a new Zenodo draft deposit",
        formatter_class=fmt,
    )

    create_deposit_parser.add_argument(
        "-c",
        "--concept_id",
        type=int,
        default=None,
        help="A Zenodo concept ID",
    )

    # DISCARD DEPOSIT subcommand
    discard_deposit_desc = """
    Discard an unpublished deposit. The deposit and all uploaded files will be
    removed from Zenodo.
    """

    discard_deposit_parser = subparsers.add_parser(
        "discard_deposit",
        description=textwrap.dedent(discard_deposit_desc),
        help="Discard an unpublished deposit",
        formatter_class=fmt,
    )

    discard_deposit_parser.add_argument(
        "zenodo_json",
        type=str,
        default=None,
        help="Path to a Zenodo metadata file for the deposit to discard",
    )

    # GET DEPOSIT subcommand
    get_deposit_desc = """
    Download the Zenodo metadata for a deposit and print out summary information.
    """

    get_deposit_parser = subparsers.add_parser(
        "get_deposit",
        description=textwrap.dedent(get_deposit_desc),
        help="Download and display deposit metadata",
    )
    get_deposit_parser.add_argument(
        "zenodo_id",
        type=int,
        default=None,
        help="An ID for an existing Zenodo deposit",
    )

    # PUBLISH DEPOSIT subcommand
    publish_deposit_desc = """
    Publishes a Zenodo deposit. This is the final step in publishing a dataset and
    is not reversible. Once a dataset is published, the DOI associated with the
    record is published to Datacite.

    It may be worth reviewing the deposit webpage (https://zenodo.org/deposit/###)
    before finally publishing.
    """

    publish_deposit_parser = subparsers.add_parser(
        "publish_deposit",
        description=textwrap.dedent(publish_deposit_desc),
        help="Publish a draft deposit",
        formatter_class=fmt,
    )

    publish_deposit_parser.add_argument(
        "zenodo_json", type=str, help="Path to a Zenodo metadata file"
    )

    # UPLOAD FILE subcommand
    upload_file_desc = """
    Uploads the contents of a specified file to an _unpublished_ Zenodo deposit,
    optionally using an alternative filename. If you upload a new file to the same
    filename, it will replace the existing uploaded file.
    """

    upload_file_parser = subparsers.add_parser(
        "upload_file",
        description=textwrap.dedent(upload_file_desc),
        help="Upload a file to an unpublished deposit",
        formatter_class=fmt,
    )

    upload_file_parser.add_argument(
        "zenodo_json",
        type=str,
        default=None,
        help="Path to a Zenodo metadata JSON for the deposit",
    )
    upload_file_parser.add_argument(
        "filepath", type=str, default=None, help="The path to the file to be uploaded"
    )
    upload_file_parser.add_argument(
        "--zenodo_filename",
        type=str,
        default=None,
        help="An optional alternative file name to be used on Zenodo",
    )

    # DELETE FILE subcommand
    delete_file_desc = """
    Delete an uploaded file from an unpublished deposit. The deposit metadata will
    be re-downloaded to ensure an up to date list of files in the deposit.
    """

    delete_file_parser = subparsers.add_parser(
        "delete_file",
        description=textwrap.dedent(delete_file_desc),
        help="Delete a file from an unpublished deposit",
        formatter_class=fmt,
    )

    delete_file_parser.add_argument(
        "zenodo_json",
        type=str,
        default=None,
        help="Path to a Zenodo metadata file for the deposit to discard",
    )
    delete_file_parser.add_argument(
        "filename", type=str, default=None, help="The name of the file to delete"
    )

    # UPLOAD METADATA subcommand
    upload_metadata_desc = """
    Uses the dataset metadata created using `safedata_validate` to populate the
    required Zenodo metadata for an unpublished deposit.
    """

    upload_metadata_parser = subparsers.add_parser(
        "upload_metadata",
        description=textwrap.dedent(upload_metadata_desc),
        help="Populate the Zenodo metadata",
        formatter_class=fmt,
    )

    upload_metadata_parser.add_argument(
        "zenodo_json", type=str, help="Path to a Zenodo metadata file"
    )
    upload_metadata_parser.add_argument(
        "dataset_json", type=str, help="Path to a dataset metadata file"
    )

    # AMEND METADATA subcommand
    amend_metadata_desc = """
    Updates the Zenodo metadata for an published deposit. To use this, make sure
    you have the most recent Zenodo metadata for the deposit and then edit the
    JSON file to the new values. You can only edit the contents of the
    `metadata` section.

    Caution: this command should only be used to make urgent changes - such as
    access restrictions. It is also easy to submit invalid metadata!
    """

    amend_metadata_parser = subparsers.add_parser(
        "amend_metadata",
        description=textwrap.dedent(amend_metadata_desc),
        help="Update published Zenodo metadata",
        formatter_class=fmt,
    )

    amend_metadata_parser.add_argument(
        "zenodo_json_update", type=str, help="Path to an updated Zenodo metadata file"
    )

    # SYNC LOCAL DIR subcommand

    sync_local_dir_desc = """
    Synchronize a local data directory

    This subcommand allows a safedata developer or community maintainer to
    create or update such a directory with _all_ of the resources in the Zenodo
    community, regardless of their public access status. This forms a backup
    (although Zenodo is heavily backed up) but also provides local copies of the
    files for testing and development of the code packages.

    The file structure of the directory follows that used by the safedata R
    package, used to store metadata and files downloaded from a safedata
    community on Zenodo and from a safedata metadata server. The
    `safedata_validator` configuration file will need to include the metadata
    API.

    By default, only the XLSX files containing metadata and data tables are
    downloaded, ignoring any additional files, which are often large.
    """

    sync_local_dir_parser = subparsers.add_parser(
        "sync_local_dir",
        description=textwrap.dedent(sync_local_dir_desc),
        help="Create or update a local safedata directory",
        formatter_class=fmt,
    )

    sync_local_dir_parser.add_argument(
        "datadir",
        type=str,
        help="The path to a local directory containing "
        "an existing safedata directory or an empty folder in which to create one",
    )

    sync_local_dir_parser.add_argument(
        "--not-just-xlsx",
        action="store_true",
        default=False,
        help="Should large non-xlsx files also be downloaded.",
    )

    sync_local_dir_parser.add_argument(
        "--replace-modified",
        action="store_true",
        default=False,
        help="Should locally modified files be overwritten with the archive version",
    )
    # MAINTAIN RIS subcommand

    maintain_ris_desc = """
    This command maintains a RIS format bibliography file of the datasets
    uploaded to a Zenodo community. It can update an existing RIS format file
    to add new records or it can create the file from scratch.

    The program uses both the Zenodo API (to find the records in the community)
    and the Datacite API to access machine readable bibliographic records.
    """

    maintain_ris_parser = subparsers.add_parser(
        "maintain_ris",
        description=textwrap.dedent(maintain_ris_desc),
        help="Maintain a RIS bibliography file for datasets",
        formatter_class=fmt,
    )

    # positional argument inputs
    maintain_ris_parser.add_argument(
        "ris_file",
        type=str,
        help="The file path to populate with RIS records. If this file "
        "already exists, it is assumed to be RIS file to update "
        "with any new records not already included in the file.",
    )

    # GENERATE HTML subcommand
    generate_html_desc = """
    Generates an html file containing a standard description of a dataset from the JSON
    metadata. Usually this will be generated and uploaded as part of the dataset
    publication process, but this subcommand can be used for local checking of the
    resulting HTML.
    """

    generate_html_parser = subparsers.add_parser(
        "generate_html",
        description=textwrap.dedent(generate_html_desc),
        help="Generate an HTML dataset description",
    )

    generate_html_parser.add_argument(
        "zenodo_json", type=str, help="Path to a Zenodo metadata file"
    )
    generate_html_parser.add_argument(
        "dataset_json", type=str, help="Path to a dataset metadata file"
    )

    # GENERATE XML subcommand
    generate_xml_desc = """
    Creates an INSPIRE compliant XML metadata file for a published dataset,
    optionally including a user provided lineage statement (such as project
    details).
    """

    generate_xml_parser = subparsers.add_parser(
        "generate_xml",
        description=textwrap.dedent(generate_xml_desc),
        help="Create INSPIRE compliant metadata XML",
        formatter_class=fmt,
    )

    generate_xml_parser.add_argument(
        "zenodo_json", type=str, help="Path to a Zenodo metadata file"
    )
    generate_xml_parser.add_argument(
        "dataset_json", type=str, help="Path to a dataset metadata file"
    )
    generate_xml_parser.add_argument(
        "-l",
        "--lineage-statement",
        type=str,
        help="Path to a text file containing a lineage statement",
    )

    # SHOW CONFIG subcommand
    show_config_desc = """
    Loads the safedata_validator configuration file, displays the config setup
    being used for safedata_zenodo and then exits.
    """

    # This has no arguments, so subparser object not needed
    _ = subparsers.add_parser(
        "show_config",
        description=textwrap.dedent(show_config_desc),
        help="Report the config being used and exit",
        formatter_class=fmt,
    )

    # ------------------------------------------------------
    # Parser definition complete - now handle the inputs
    # ------------------------------------------------------

    # Parse the arguments and set the verbosity
    args = parser.parse_args()

    if args.subcommand is None:
        parser.print_usage()
        return

    # Show the package config and exit if requested
    if args.subcommand == "show_config":
        resources = Resources(args.resources)
        print("\nZenodo configuration:")
        for key, val in resources.zenodo.items():
            print(f" - {key}: {val}")
        print("\nMetadata server configuration:")
        for key, val in resources.metadata.items():
            print(f" - {key}: {val}")
        return

    handler = get_handler()
    if args.quiet:
        # Don't suppress error messages
        handler.setLevel("ERROR")
    else:
        handler.setLevel("DEBUG")

    resources = Resources(args.resources)

    # Handle the remaining subcommands
    if args.subcommand in ["create_deposit", "cdep"]:
        # Run the command
        response, error = create_deposit(
            concept_id=args.concept_id, resources=resources
        )
        # Trap errors
        if error is not None:
            LOGGER.error(f"Failed to create deposit: {error}")
            return

        # Output the response as a deposit  JSON file
        rec_id = response["record_id"]
        LOGGER.info(f"Created deposit: {rec_id}")
        outfile = os.path.join(os.getcwd(), f"zenodo_{rec_id}.json")

        with open(outfile, "w") as outf:
            simplejson.dump(response, outf)
            LOGGER.info(f"Zenodo deposit metadata downloaded to: {outfile}")

    elif args.subcommand in ["discard_deposit", "ddep"]:
        # Load the Zenodo deposit JSON, which contains API links
        with open(args.zenodo_json) as dep_json:
            metadata = simplejson.load(dep_json)

        # Run the command
        response, error = discard_deposit(metadata=metadata, resources=resources)

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to discard deposit: {error}")
        else:
            LOGGER.info("Deposit discarded")

    elif args.subcommand in ["get_deposit", "gdep"]:
        # Run the command
        response, error = get_deposit(deposit_id=args.zenodo_id, resources=resources)

        if error is not None:
            LOGGER.error(f"Failed to get info: {error}")
            return

        # Dump the response to a JSON file
        outfile = os.path.join(os.getcwd(), f"zenodo_{args.zenodo_id}.json")
        with open(outfile, "w") as outf:
            simplejson.dump(response, outf)

        # Print a short summary
        LOGGER.info(f"Record ID: {response['record_id']}")
        LOGGER.info(f"Concept ID: {response['conceptrecid']}")
        LOGGER.info(
            f"Status: {'published' if response['submitted'] else 'Not published'}"
        )
        LOGGER.info(f"Title: {response['title'] if response['title'] else 'Not set'}")

        if len(response["files"]) == 0:
            LOGGER.info("Files: none uploaded")
        else:
            LOGGER.info("Files:")
            FORMATTER.push()
            for this_file in response["files"]:
                LOGGER.info(
                    f"{this_file['filename']} ({this_file['filesize']} bytes, "
                    f"md5: {this_file['checksum']})"
                )
            FORMATTER.pop()

        LOGGER.info(f"Metadata downloaded to: {outfile}")

    elif args.subcommand in ["publish_deposit", "pdep"]:
        with open(args.zenodo_json) as zn_json:
            zenodo_json_data = simplejson.load(zn_json)

        # Run the function
        response, error = publish_deposit(zenodo=zenodo_json_data, resources=resources)

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to publish deposit: {error}")
            return

        # Update the Zenodo JSON file with publication details
        LOGGER.info(f"Published to: {response['links']['record']}")
        with open(args.zenodo_json, "w") as zn_json:
            simplejson.dump(response, zn_json)
            LOGGER.info("Zenodo metadata updated")

    elif args.subcommand in ["upload_file", "ufile"]:
        # Load the Zenodo deposit JSON, which contains API links
        with open(args.zenodo_json) as dep_json:
            metadata = simplejson.load(dep_json)

        # Report on proposed upload
        if args.zenodo_filename is None:
            LOGGER.info(f"Uploading: {os.path.basename(args.filepath)}")
        else:
            LOGGER.info(
                f"Uploading: {os.path.basename(args.filepath)} as "
                f"{args.zenodo_filename}"
            )

        # Run the command
        response, error = upload_file(
            metadata=metadata,
            filepath=args.filepath,
            zenodo_filename=args.zenodo_filename,
            resources=resources,
            progress_bar=not args.quiet,
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to upload file: {error}")
        else:
            LOGGER.info("File uploaded")

    elif args.subcommand in ["delete_file", "dfile"]:
        # Load the Zenodo deposit JSON, which contains API links
        with open(args.zenodo_json) as dep_json:
            metadata = simplejson.load(dep_json)

        # Run the command
        response, error = delete_file(
            metadata=metadata, filename=args.filename, resources=resources
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to delete file: {error}")
        else:
            LOGGER.info("File deleted")

    elif args.subcommand in ["upload_metadata", "umeta"]:
        # Open the two JSON files
        with open(args.dataset_json) as ds_json:
            dataset_json = simplejson.load(ds_json)

        with open(args.zenodo_json) as zn_json:
            zenodo_json = simplejson.load(zn_json)

        # Run the function
        response, error = upload_metadata(
            metadata=dataset_json, zenodo=zenodo_json, resources=resources
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to add metadata file: {error}")
        else:
            LOGGER.info("Metadata uploaded")

    elif args.subcommand in ["amend_metadata", "ameta"]:
        with open(args.deposit_json_update) as zn_json_update:
            zenodo_json_update = simplejson.load(zn_json_update)

        # Run the function
        response, error = update_published_metadata(
            zenodo_md=zenodo_json_update,
            resources=resources,
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to update published metadata: {error}")
        else:
            LOGGER.info("Metadata updated")

    elif args.subcommand in ["sync_local_dir", "sync"]:
        sync_local_dir(
            datadir=args.datadir,
            xlsx_only=not args.not_just_xlsx,
            resources=resources,
            replace_modified=args.replace_modified,
        )

    elif args.subcommand in ["maintain_ris", "ris"]:
        # Run the download RIS data function
        download_ris_data(ris_file=args.ris_file, resources=resources)

    elif args.subcommand in ["generate_html", "html"]:
        # Run the download RIS data function
        with open(args.dataset_json) as ds_json:
            dataset_json = simplejson.load(ds_json)

        desc = dataset_description(metadata=dataset_json)

        with open(args.html_out, "w") as outf:
            outf.write(desc)

    elif args.subcommand in ["generate_xml", "xml"]:
        # Open the two JSON files
        with open(args.dataset_json) as ds_json:
            dataset_json = simplejson.load(ds_json)

        with open(args.zenodo_json) as zn_json:
            zenodo_json = simplejson.load(zn_json)

        if args.lineage_statement is not None:
            with open(args.lineage_statement) as lin_file:
                lineage_statement = lin_file.readlines()
        else:
            lineage_statement = None

        # Run the function
        response, error = generate_inspire_xml(
            metadata=dataset_json,
            zenodo=zenodo_json,
            resources=resources,
            lineage_statement=lineage_statement,
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to generate INSPIRE xml: {error}")
        else:
            LOGGER.info("Inspire XML generated")

    return


def _safedata_server_cli():
    """Post updated information to a safedata server instance.

    This command line tool provides functions to update a web server running
    the safedata server API. This includes updating the gazetteer data used
    for the instance and updating the database of published datasets.

    To use these tools, the safedata_validator Resources configuration must
    contain the URL for the server and an access token.
    """

    # create the top-level parser and configure to take subparsers
    desc = textwrap.dedent(_safedata_server_cli.__doc__)
    fmt = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, formatter_class=fmt)

    parser.add_argument(
        "-r",
        "--resources",
        type=str,
        default=None,
        help="Path to a safedata_validator resource configuration file",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress normal information messages. ",
    )

    subparsers = parser.add_subparsers(dest="subcommand", metavar="")

    # POST METADATA subcommand

    post_metadata_desc = """
    This function posts dataset and Zenodo metadata for a validated and published
    dataset to an server providing the API for accessing and searching safedata
    metadata. The metadata server then handles populating a database with the
    metadata and exposing that data via the server API.
    """

    post_metadata_parser = subparsers.add_parser(
        "post_metadata",
        description=textwrap.dedent(post_metadata_desc),
        help="Post dataset metadata to a safedata server",
        formatter_class=fmt,
    )

    # positional argument inputs
    post_metadata_parser.add_argument(
        "zenodo_json", type=str, help="Path to a Zenodo metadata file"
    )
    post_metadata_parser.add_argument(
        "dataset_json", type=str, help="Path to a dataset metadata file"
    )

    # UPDATE GAZETTEER subcommand

    update_gazetteer_desc = """
    This function updates the gazetteer resources being used by a safedata server
    instance. It uploads both the gazetteer geojson file and location aliases
    file and the server then validates the file contents.
    """

    update_gazetteer_parser = subparsers.add_parser(
        "update_gazetteer",
        description=textwrap.dedent(update_gazetteer_desc),
        help="Update the gazetteer data on safedata server",
        formatter_class=fmt,
    )

    # positional argument inputs
    update_gazetteer_parser.add_argument(
        "gazetteer_json", type=str, help="Path to a gazetteer geojson file"
    )
    update_gazetteer_parser.add_argument(
        "location_aliases", type=str, help="Path to a CSV of location aliases"
    )

    # ------------------------------------------------------
    # Parser definition complete - now handle the inputs
    # ------------------------------------------------------

    # Parse the arguments and set the verbosity
    args = parser.parse_args()

    if args.subcommand is None:
        parser.print_usage()
        return

    # Show the package config and exit if requested
    if args.subcommand == "show_config":
        resources = Resources(args.resources)
        print("\nZenodo configuration:")
        for key, val in resources.zenodo.items():
            print(f" - {key}: {val}")
        print("\nMetadata server configuration:")
        for key, val in resources.metadata.items():
            print(f" - {key}: {val}")
        return

    handler = get_handler()
    if args.quiet:
        # Don't suppress error messages
        handler.setLevel("ERROR")
    else:
        handler.setLevel("DEBUG")

    resources = Resources(args.resources)

    # Handle the remaining subcommands
    if args.subcommand in ["post_metadata", "cdep"]:
        # Open the two JSON files
        with open(args.dataset_json) as ds_json:
            dataset_json = simplejson.load(ds_json)

        with open(args.zenodo_json) as zn_json:
            zenodo_json = simplejson.load(zn_json)

        # Run the function
        response, error = post_metadata(
            metadata=dataset_json, zenodo=zenodo_json, resources=resources
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to post metadata: {error}")
        else:
            LOGGER.info("Metadata posted")

        return

    if args.subcommand in ["update_gazetteer"]:
        # Open the two JSON files
        with open(args.gazetteer_json) as gz_json:
            gazetteer_json = simplejson.load(gz_json)

        with open(args.location_aliases) as aliases_csv:
            location_aliases = aliases_csv.read()

        # Run the function
        response, error = update_gazetteer(
            gazetteer=gazetteer_json,
            location_aliases=location_aliases,
            resources=resources,
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to update gazetteer: {error}")
        else:
            LOGGER.info("Gazetteer updated")

        return


# Local Database building


def _build_local_gbif_cli():
    """Build a local GBIF database.

    This tool builds an SQLite database of the GBIF backbone taxonomy to use
    in validation by safedata_validate. There are multiple versions of the
    dataset, and the available versions can be seen here:

        https://hosted-datasets.gbif.org/datasets/backbone/

    The tool will optionally take a timestamp - using the format '2021-11-26'
    - to build a particular version, but defaults to the most recent version.
    """

    desc = textwrap.dedent(_build_local_gbif_cli.__doc__)
    fmt = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, formatter_class=fmt)

    parser.add_argument("outdir", help="Location to create database file.")
    parser.add_argument(
        "-t",
        "--timestamp",
        default=None,
        type=str,
        help="The time stamp of a database archive version to use.",
    )

    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as download_loc:
        file_data = download_gbif_backbone(
            outdir=download_loc, timestamp=args.timestamp
        )
        build_local_gbif(outdir=args.outdir, **file_data)

    return


def _build_local_ncbi_cli():
    """Build a local NCBI database.

    This tool builds an SQLite database of the NCBI  taxonomy to use in
    validation by safedata_validate. There are multiple archived versions
    of the dataset, and the available versions can be seen here:

        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/

    The tool will optionally take a timestamp - using the format '2021-11-26'
    - to build a particular version, but defaults to the most recent version.
    """

    desc = textwrap.dedent(_build_local_ncbi_cli.__doc__)
    fmt = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, formatter_class=fmt)

    parser.add_argument("outdir", help="Location to create database file.")
    parser.add_argument(
        "-t",
        "--timestamp",
        default=None,
        type=str,
        help="The time stamp of a database archive version to use.",
    )

    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as download_loc:
        file_data = download_ncbi_taxonomy(
            outdir=download_loc, timestamp=args.timestamp
        )
        build_local_ncbi(outdir=args.outdir, **file_data)

    return

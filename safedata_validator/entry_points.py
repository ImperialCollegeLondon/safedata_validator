"""Provide command line scripts.

This module provides functions exposed as command line entry points:

* ``_safedata_validate_cli``, exposed as `safedata_validate`
* ``_safedata_zenodo_cli``, exposed as `safedata_zenodo`
* ``_safedata_metadata_cli``, exposed as `safedata_metadata`
* ``_build_local_gbif_cli``, exposed as `safedata_build_local_gbif`
* ``_build_local_ncbi_cli``, exposed as `safedata_build_local_ncbi`
"""

import argparse
import datetime
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
from safedata_validator.server import post_metadata, update_resources
from safedata_validator.taxondb import (
    build_local_gbif,
    build_local_ncbi,
    download_gbif_backbone,
    download_ncbi_taxonomy,
    get_gbif_version,
    get_ncbi_version,
)
from safedata_validator.zenodo import (
    create_deposit,
    dataset_description,
    delete_files,
    discard_deposit,
    download_ris_data,
    generate_inspire_xml,
    get_deposit,
    publish_dataset,
    publish_deposit,
    sync_local_dir,
    update_published_metadata,
    upload_files,
    upload_metadata,
)


def _desc_formatter(prog):
    """Bespoke argparse description formatting."""
    return argparse.RawDescriptionHelpFormatter(prog, max_help_position=16)


def _show_resources(resources: Resources) -> None:
    """Simple function to print out the resources to the display."""

    print("\nValidation extents:")
    for key, val in resources.extents.items():
        # handle date formatting
        if val is not None and isinstance(val[0], datetime.date):
            val = [v.isoformat() for v in val]
        print(f" - {key}: {val}")
    print("\nZenodo configuration:")
    for key, val in resources.zenodo.items():
        print(f" - {key}: {val}")
    print("\nMetadata server configuration:")
    for key, val in resources.metadata.items():
        print(f" - {key}: {val}")


def _safedata_validate_cli(args_list: list[str] | None = None) -> int:
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
    metadata will be saved. This is used in the dataset publication process.
    By default, the JSON file is saved to the same directory as the input
    file, using the same filename but with the `.json` extension. This can
    be saved elsewhere using the `--json` option.

    The command also outputs a log of the validation process, which
    identifies validation issues. This defaults to being written to
    stderr but can be redirected to a file using the `--log` option.

    Args:
        args_list: This is a developer option used to simulate command line usage by
            providing a list of command line argument strings to the entry point
            function. For example, ``safedata_validate --show-resources`` can be
            replicated by calling ``_safedata_validate_cli(['--show-resources'])``.

    Returns:
        An integer code showing success (0) or failure (1).
    """

    # If no arguments list is provided
    if args_list is None:
        args_list = sys.argv[1:]

    # Check function docstring exists to safeguard against -OO mode, and strip off the
    # description of the function args_list, which should not be included in the command
    # line docs
    if _safedata_validate_cli.__doc__ is not None:
        desc = textwrap.dedent(
            "\n".join(_safedata_validate_cli.__doc__.splitlines()[:-10])
        )
    else:
        desc = "Python in -OO mode: no docs"

    parser = argparse.ArgumentParser(
        prog="safedata_validate",
        description=desc,
        formatter_class=_desc_formatter,
    )

    # Filename is 'optional' in that show resources can run without a filename
    parser.add_argument(
        "filename", help="Path to the Excel file to be validated", nargs="?"
    )

    parser.add_argument(
        "-r",
        "--resources",
        type=str,
        default=None,
        help="A path to a resources configuration file",
    )

    parser.add_argument(
        "-s",
        "--show-resources",
        action="store_true",
        default=False,
        help="Validate and display the selected resources and exit",
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
        "-l",
        "--log",
        default=None,
        type=Path,
        help="Save the validation log to a file, not print to the console.",
    )

    parser.add_argument(
        "-j",
        "--json",
        default=None,
        type=Path,
        help="An optional output path for the validated dataset JSON.",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    args = parser.parse_args(args=args_list)

    # Load resources
    resources = Resources(args.resources)

    # Show the resources and exit if requested
    if args.show_resources:
        _show_resources(resources)
        return 0

    # Catch the possibility that filename is empty, which is only allowed if
    # --show-resources is specified
    if args.filename is None:
        parser.print_usage()
        print(
            "safedata_validate: error: the following arguments are required: filename"
        )
        return 1

    # Configure the logging location
    if args.log is None:
        use_stream_logging()
    else:
        use_file_logging(args.log)

    # Create the dataset and load from workbook
    ds = Dataset(resources=resources)
    ds.load_from_workbook(
        filename=args.filename,
        validate_doi=args.validate_doi,
        chunk_size=args.chunk_size,
    )

    if not ds.passed:
        sys.stdout.write("File validation failed\n")
        return 1

    # Output JSON file
    json_file = args.json or os.path.splitext(args.filename)[0] + ".json"
    with open(json_file, "w") as json_out:
        json_out.write(ds.to_json())

    sys.stdout.write("File validation passed\n")
    sys.stdout.write(f"JSON metadata written to {json_file}\n")
    return 0


def _safedata_zenodo_cli(args_list: list[str] | None = None) -> int:
    """Publish validated datasets to Zenodo using a command line interface.

    This is a the command line interface for publishing safedata validated
    datasets to Zenodo, downloading information and maintaining a local
    copy of the datasets in the file structure required by the R safedata
    package.

    The safedata_zenodo command is used by providing subcommands for the
    different actions required to publish a validated dataset. The list of
    subcommands is shown below and individual help is available for each
    of the subcommands:

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

    The subcommands that send and receive data from Zenodo also accept the
    `--live` and `--sandbox` options which can be used to override the
    `use_sandbox` setting in the configuration file. If the configuration is set
    to `true` then `--live` will use the live site and if the configuration is set
    to `false` then the `--sandbox` can be used to use the sandbox instead.

    Note that most of these actions are also available via the Zenodo website.

    Args:
        args_list: This is a developer option used to simulate command line usage by
            providing a list of command line argument strings to the entry point
            function. For example, ``safedata_zenodo create_deposit`` can be
            replicated by calling ``_safedata_zenodo_cli(['create_deposit'])``.

    Returns:
        An integer code showing success (0) or failure (1).
    """

    # If no arguments list is provided
    if args_list is None:
        args_list = sys.argv[1:]

    # Check function docstring exists to safeguard against -OO mode, and strip off the
    # description of the function args_list, which should not be included in the command
    # line docs
    if _safedata_zenodo_cli.__doc__ is not None:
        desc = textwrap.dedent(
            "\n".join(_safedata_zenodo_cli.__doc__.splitlines()[:-10])
        )
    else:
        desc = "Python in -OO mode: no docs"

    # create the top-level parser and configure to take subparsers
    parser = argparse.ArgumentParser(
        prog="safedata_zenodo",
        description=desc,
        formatter_class=_desc_formatter,
    )

    parser.add_argument(
        "-r",
        "--resources",
        type=str,
        default=None,
        help="Path to a safedata_validator resource configuration file",
    )

    parser.add_argument(
        "-s",
        "--show-resources",
        action="store_true",
        default=False,
        help="Validate and display the selected resources and exit",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress normal information messages. ",
    )

    # Add a shared mutually exclusive switch between the sandbox and live Zenodo sites,
    # which will override the values set in the configuration
    sandbox_switches = argparse.ArgumentParser(add_help=False)
    switch_group = sandbox_switches.add_mutually_exclusive_group()
    switch_group.add_argument(
        "--live",
        dest="sandbox",
        action="store_false",
        default=None,
        help="Use the Zenodo live site, overriding the configuration file",
    )
    switch_group.add_argument(
        "--sandbox",
        dest="sandbox",
        action="store_true",
        default=None,
        help="Use the Zenodo sandbox site, overriding the configuration file",
    )

    # Create subparsers to add shared positional arguments across actions - these can be
    # added to individual actions via the parents argument
    parse_zenodo_metadata = argparse.ArgumentParser(add_help=False)
    parse_zenodo_metadata.add_argument(
        "zenodo_json",
        type=str,
        default=None,
        help="Path to a Zenodo JSON file for a deposit",
    )

    parse_dataset_metadata = argparse.ArgumentParser(add_help=False)
    parse_dataset_metadata.add_argument(
        "dataset_json",
        type=str,
        default=None,
        help="Path to a JSON metadata file for a dataset",
    )

    subparsers = parser.add_subparsers(dest="subcommand", metavar="")

    # CREATE DEPOSIT subcommand
    create_deposit_desc = """
    Create a new deposit draft. 
    
    The new version option takes the record ID of the most recent version of an
    existing dataset and creates a new deposit as a new version of that dataset. 
    Versions of datasets are grouped under a single concept ID, which always 
    redirects to the most recent version. Use the most recent version ID and 
    _not_ the concept ID here.

    When successful, the function downloads and saves a JSON file containing the
    resulting Zenodo deposit metadata. This file is used as an input to other
    subcommands that work with an existing deposit.

    The --id-to-stdout option can be provided to explicitly return the new
    deposit ID to stdout, where it can be captured for use in shell scripts.
    All other logging is written to stderr.
    """
    create_deposit_parser = subparsers.add_parser(
        "create_deposit",
        description=textwrap.dedent(create_deposit_desc),
        help="Create a new Zenodo draft deposit",
        formatter_class=_desc_formatter,
        parents=[sandbox_switches],
    )

    create_deposit_parser.add_argument(
        "-n",
        "--new-version",
        type=int,
        default=None,
        help="Create a new version of the dataset with the provided Zenodo ID.",
    )

    create_deposit_parser.add_argument(
        "-i",
        "--id-to-stdout",
        action="store_true",
        default=False,
        help="Write the deposit record ID to stdout.",
    )

    # DISCARD DEPOSIT subcommand
    discard_deposit_desc = """
    Discard an unpublished deposit. The deposit and all uploaded files will be
    removed from Zenodo.
    """

    discard_deposit_parser = subparsers.add_parser(  # noqa: F841
        "discard_deposit",
        description=textwrap.dedent(discard_deposit_desc),
        help="Discard an unpublished deposit",
        formatter_class=_desc_formatter,
        parents=[parse_zenodo_metadata, sandbox_switches],
    )

    # GET DEPOSIT subcommand
    get_deposit_desc = """
    Download the Zenodo metadata for a deposit and print out summary information.
    """

    get_deposit_parser = subparsers.add_parser(
        "get_deposit",
        description=textwrap.dedent(get_deposit_desc),
        help="Download and display deposit metadata",
        parents=[sandbox_switches],
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

    publish_deposit_parser = subparsers.add_parser(  # noqa: F841
        "publish_deposit",
        description=textwrap.dedent(publish_deposit_desc),
        help="Publish a draft deposit",
        formatter_class=_desc_formatter,
        parents=[parse_zenodo_metadata, sandbox_switches],
    )

    # UPLOAD FILES subcommand
    upload_files_desc = """
    Uploads a set of files to an _unpublished_ Zenodo deposit. If you upload a new file
    to the same filename, it will replace the existing uploaded file.
    """

    upload_files_parser = subparsers.add_parser(
        "upload_files",
        description=textwrap.dedent(upload_files_desc),
        help="Upload a file to an unpublished deposit",
        formatter_class=_desc_formatter,
        parents=[parse_zenodo_metadata, sandbox_switches],
    )

    upload_files_parser.add_argument(
        "filepaths",
        type=str,
        default=None,
        help="The paths to the file to be uploaded",
        nargs="*",
    )

    # DELETE FILE subcommand
    delete_files_desc = """
    Delete an list of uploaded file from an unpublished deposit. The deposit metadata
    will be re-downloaded to ensure an up to date list of files in the deposit.
    """

    delete_files_parser = subparsers.add_parser(
        "delete_files",
        description=textwrap.dedent(delete_files_desc),
        help="Delete files from an unpublished deposit",
        formatter_class=_desc_formatter,
        parents=[parse_zenodo_metadata, sandbox_switches],
    )

    delete_files_parser.add_argument(
        "filenames",
        type=str,
        default=None,
        help="The names of files to delete",
        nargs="*",
    )

    # UPLOAD METADATA subcommand
    upload_metadata_desc = """
    Uses the dataset metadata created using `safedata_validate` to populate the
    required Zenodo metadata for an unpublished deposit.
    """

    upload_metadata_parser = subparsers.add_parser(  # noqa: F841
        "upload_metadata",
        description=textwrap.dedent(upload_metadata_desc),
        help="Populate the Zenodo metadata",
        formatter_class=_desc_formatter,
        parents=[parse_zenodo_metadata, parse_dataset_metadata, sandbox_switches],
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
        formatter_class=_desc_formatter,
        parents=[sandbox_switches],
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

    # This does not use the sandbox switches - doesn't make a great deal of sense to
    # sync a sandbox collection.
    sync_local_dir_parser = subparsers.add_parser(
        "sync_local_dir",
        description=textwrap.dedent(sync_local_dir_desc),
        help="Create or update a local safedata directory",
        formatter_class=_desc_formatter,
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

    # This does not use the sandbox switches - doesn't make a great deal of sense to
    # sync a sandbox collection.
    maintain_ris_parser = subparsers.add_parser(
        "maintain_ris",
        description=textwrap.dedent(maintain_ris_desc),
        help="Maintain a RIS bibliography file for datasets",
        formatter_class=_desc_formatter,
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
    resulting HTML and developing custom templates.
    """

    # This does not use the sandbox switches - not applicable
    generate_html_parser = subparsers.add_parser(
        "generate_html",
        description=textwrap.dedent(generate_html_desc),
        help="Generate an HTML dataset description",
        parents=[parse_dataset_metadata],
    )

    generate_html_parser.add_argument(
        "html_out",
        type=str,
        default=None,
        help="Output path for the HTML file",
    )

    # GENERATE XML subcommand
    generate_xml_desc = """
    Creates an INSPIRE compliant XML metadata file for a published dataset,
    optionally including a user provided lineage statement (such as project
    details).
    """

    # This does not use the sandbox switches - not applicable
    generate_xml_parser = subparsers.add_parser(
        "generate_xml",
        description=textwrap.dedent(generate_xml_desc),
        help="Create INSPIRE compliant metadata XML",
        formatter_class=_desc_formatter,
        parents=[parse_zenodo_metadata, parse_dataset_metadata],
    )

    generate_xml_parser.add_argument(
        "-l",
        "--lineage-statement",
        type=str,
        help="Path to a text file containing a lineage statement",
    )

    generate_xml_parser.add_argument(
        "xml_out",
        type=str,
        default=None,
        help="Output path for the XML file",
    )

    # PUBLISH DATASET subcommand
    publish_dataset_desc = """
    This subcommand runs through the complete publication process for a validated
    dataset and any external files. It does not provide all of the options of the
    subcommands for individual steps but covers the main common usage of the
    safedata_zenodo command. If the publication process fails, the resulting partial
    deposit is discarded. All of the files in the deposit must be replaced.
    """
    publish_dataset_parser = subparsers.add_parser(
        "publish_dataset",
        description=textwrap.dedent(publish_dataset_desc),
        help="Publish a validated dataset",
        formatter_class=_desc_formatter,
        parents=[parse_dataset_metadata, sandbox_switches],
    )

    publish_dataset_parser.add_argument(
        "dataset",
        type=Path,
        help="Path to the Excel file to be published",
    )

    publish_dataset_parser.add_argument(
        "-n",
        "--new-version",
        type=int,
        default=None,
        help="Create a new version of the dataset with the provided Zenodo ID.",
    )

    publish_dataset_parser.add_argument(
        "--no-xml",
        action="store_true",
        default=False,
        help="Do not include metadata XML in the published record",
    )

    publish_dataset_parser.add_argument(
        "-e",
        "--external-files",
        default=[],
        type=Path,
        nargs="*",
        help="A set of external files documented in the dataset to be included.",
    )

    # ------------------------------------------------------
    # Parser definition complete - now handle the inputs
    # ------------------------------------------------------

    # Parse the arguments
    args = parser.parse_args(args=args_list)

    # Load resources
    resources = Resources(args.resources)

    # Show the resources and exit if requested
    if args.show_resources:
        _show_resources(resources)
        return 0

    if args.subcommand is None:
        parser.print_usage()
        return 0

    # Handle the sandbox switches:
    # If the command parser provides the sandbox argument and one of the flags has been
    # set (the default is None) then update the loaded resources
    if "sandbox" in args and args.sandbox is not None:
        resources.zenodo.use_sandbox = args.sandbox

    # Set the verbosity
    handler = get_handler()
    if args.quiet:
        # Don't suppress error messages
        handler.setLevel("ERROR")
    else:
        handler.setLevel("DEBUG")

    # Handle the remaining subcommands
    if args.subcommand == "create_deposit":
        # Run the command
        response, error = create_deposit(
            new_version=args.new_version, resources=resources
        )
        # Trap errors
        if error is not None:
            LOGGER.error(f"Failed to create deposit: {error}")
            return 1

        # Output the response as a deposit  JSON file
        rec_id = response["record_id"]
        LOGGER.info(f"Created deposit: {rec_id}")
        outfile = os.path.join(os.getcwd(), f"zenodo_{rec_id}.json")

        with open(outfile, "w") as outf:
            simplejson.dump(response, outf)
            LOGGER.info(f"Zenodo deposit metadata downloaded to: {outfile}")

        # If requested, write the record ID to standard out for capture in scripts
        if args.id_to_stdout:
            sys.stdout.write(str(rec_id))

    elif args.subcommand == "discard_deposit":
        # Load the Zenodo deposit JSON, which contains API links
        with open(args.zenodo_json) as dep_json:
            metadata = simplejson.load(dep_json)

        # Run the command
        response, error = discard_deposit(zenodo_metadata=metadata, resources=resources)

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to discard deposit: {error}")
            return 1

        LOGGER.info("Deposit discarded")

    elif args.subcommand == "get_deposit":
        # Run the command
        response, error = get_deposit(deposit_id=args.zenodo_id, resources=resources)

        if error is not None:
            LOGGER.error(f"Failed to get info: {error}")
            return 1

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

    elif args.subcommand == "publish_deposit":
        with open(args.zenodo_json) as zn_json:
            zenodo_json_data = simplejson.load(zn_json)

        # Run the function
        response, error = publish_deposit(zenodo=zenodo_json_data, resources=resources)

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to publish deposit: {error}")
            return 1

        # Update the Zenodo JSON file with publication details
        LOGGER.info(f"Published to: {response['links']['record']}")
        with open(args.zenodo_json, "w") as zn_json:
            simplejson.dump(response, zn_json)
            LOGGER.info("Zenodo metadata updated")

    elif args.subcommand == "upload_files":
        # Load the Zenodo deposit JSON, which contains API links
        with open(args.zenodo_json) as dep_json:
            metadata = simplejson.load(dep_json)

        # Run the command - the response here has a different type so use a different
        # placeholder variable.
        upload_response, error = upload_files(
            metadata=metadata,
            filepaths=args.filepaths,
            resources=resources,
            progress_bar=not args.quiet,
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to upload files: {error}")
            return 1
        else:
            LOGGER.info("Files uploaded")

    elif args.subcommand == "delete_files":
        # Load the Zenodo deposit JSON, which contains API links
        with open(args.zenodo_json) as dep_json:
            metadata = simplejson.load(dep_json)

        # Run the command
        response, error = delete_files(
            metadata=metadata, filenames=args.filenames, resources=resources
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to delete files: {error}")
        else:
            LOGGER.info("Files deleted")

    elif args.subcommand == "upload_metadata":
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

    elif args.subcommand == "amend_metadata":
        with open(args.deposit_json_update) as zn_json_update:
            zenodo_json_update = simplejson.load(zn_json_update)

        # Run the function
        response, error = update_published_metadata(
            zenodo=zenodo_json_update,
            resources=resources,
        )

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to update published metadata: {error}")
        else:
            LOGGER.info("Metadata updated")

    elif args.subcommand == "sync_local_dir":
        sync_local_dir(
            datadir=args.datadir,
            xlsx_only=not args.not_just_xlsx,
            resources=resources,
            replace_modified=args.replace_modified,
        )

    elif args.subcommand == "maintain_ris":
        # Run the download RIS data function
        download_ris_data(ris_file=args.ris_file, resources=resources)

    elif args.subcommand == "generate_html":
        # Run the download RIS data function
        with open(args.dataset_json) as ds_json:
            dataset_json = simplejson.load(ds_json)

        generated_html = dataset_description(
            dataset_metadata=dataset_json, resources=resources
        )

        out_path = Path(args.html_out)
        if out_path.exists():
            LOGGER.error("HTML output file already exists")
            return 1

        with open(out_path, "w") as outf:
            outf.write(generated_html)

        LOGGER.info("HTML generated")

    elif args.subcommand == "generate_xml":
        # Open the two JSON files
        with open(args.dataset_json) as ds_json:
            dataset_json = simplejson.load(ds_json)

        with open(args.zenodo_json) as zn_json:
            zenodo_json = simplejson.load(zn_json)

        if args.lineage_statement is not None:
            with open(args.lineage_statement) as lin_file:
                lineage_statement = lin_file.read()
        else:
            lineage_statement = None

        # Run the function
        generated_xml = generate_inspire_xml(
            dataset_metadata=dataset_json,
            zenodo_metadata=zenodo_json,
            resources=resources,
            lineage_statement=lineage_statement,
        )

        out_path = Path(args.xml_out)
        if out_path.exists():
            LOGGER.error("XML output file already exists")
            return 1

        with open(out_path, "w") as outf:
            outf.write(generated_xml)

        LOGGER.info("Inspire XML generated")

    elif args.subcommand == "publish_dataset":
        with open(args.dataset_json) as ds_json:
            dataset_json = simplejson.load(ds_json)

        # Publish the dataset, trapping the possible exceptions to simply print a
        # message and return a failure exit code.
        try:
            publish_id, publish_url = publish_dataset(
                resources=resources,
                dataset=args.dataset,
                dataset_metadata=dataset_json,
                external_files=args.external_files,
                new_version=args.new_version,
            )
        except (FileNotFoundError, ValueError, RuntimeError) as excep:
            print(excep)
            return 1

    return 0


def _safedata_metadata_cli(args_list: list[str] | None = None) -> int:
    """Post updated information to a safedata server instance.

    This command line tool provides functions to update a web server running
    the safedata server API. This includes updating the gazetteer data used
    for the instance and updating the database of published datasets.

    To use these tools, the safedata_validator Resources configuration must
    contain the URL for the server and an access token.

    Args:
        args_list: This is a developer option used to simulate command line usage by
            providing a list of command line argument strings to the entry point
            function. For example, ``safedata_metadata update_resources`` can be
            replicated by calling ``_safedata_metadata_cli(['update_resources'])``.

    Returns:
        An integer code showing success (0) or failure (1).
    """

    # If no arguments list is provided
    if args_list is None:
        args_list = sys.argv[1:]

    # Check function docstring exists to safeguard against -OO mode, and strip off the
    # description of the function args_list, which should not be included in the command
    # line docs
    if _safedata_metadata_cli.__doc__ is not None:
        desc = textwrap.dedent(
            "\n".join(_safedata_metadata_cli.__doc__.splitlines()[:-10])
        )
    else:
        desc = "Python in -OO mode: no docs"

    # create the top-level parser and configure to take subparsers
    parser = argparse.ArgumentParser(
        prog="safedata_metadata",
        description=desc,
        formatter_class=_desc_formatter,
    )

    parser.add_argument(
        "-r",
        "--resources",
        type=str,
        default=None,
        help="Path to a safedata_validator resource configuration file",
    )

    parser.add_argument(
        "-s",
        "--show-resources",
        action="store_true",
        default=False,
        help="Validate and display the selected resources and exit",
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
        formatter_class=_desc_formatter,
    )

    # positional argument inputs
    post_metadata_parser.add_argument(
        "zenodo_json", type=str, help="Path to a Zenodo metadata file"
    )
    post_metadata_parser.add_argument(
        "dataset_json", type=str, help="Path to a dataset metadata file"
    )

    # UPDATE RESOURCES subcommand, which has no arguments

    update_resources_desc = """
    This function updates the resources being used by a safedata server instance. It
    uploads the gazetteer geojson file, location aliases CSV and any project database
    CSV being used. The server then validates the file contents.

    The resource file being used sets both the metadata server instance IP address and
    provides the paths to the files to be uploaded.
    """

    subparsers.add_parser(
        "update_resources",
        description=textwrap.dedent(update_resources_desc),
        help="Update the gazetteer data on safedata server",
        formatter_class=_desc_formatter,
    )

    # ------------------------------------------------------
    # Parser definition complete - now handle the inputs
    # ------------------------------------------------------

    # Parse the arguments
    args = parser.parse_args(args=args_list)

    # Load resources
    resources = Resources(args.resources)

    # Show the resources and exit if requested
    if args.show_resources:
        _show_resources(resources)
        return 0

    if args.subcommand is None:
        parser.print_usage()
        return 0

    # Set the verbosity
    handler = get_handler()
    if args.quiet:
        # Don't suppress error messages
        handler.setLevel("ERROR")
    else:
        handler.setLevel("DEBUG")

    resources = Resources(args.resources)

    # Handle the remaining subcommands
    if args.subcommand == "post_metadata":
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
            return 1

        LOGGER.info("Metadata posted")

    if args.subcommand == "update_resources":
        # Run the function
        response, error = update_resources(resources=resources)

        # Report on the outcome.
        if error is not None:
            LOGGER.error(f"Failed to update resources: {error}")
            return 1

        LOGGER.info("Resources updated")

    return 0


# Local Database building


def _build_local_gbif_cli(args_list: list[str] | None = None) -> int:
    """Build a local GBIF database.

    This tool builds an SQLite database of the GBIF backbone taxonomy to use
    in validation by safedata_validate. There are multiple versions of the
    dataset, and the available versions can be seen here:

        https://hosted-datasets.gbif.org/datasets/backbone/

    The tool will optionally take a timestamp - using the format '2021-11-26'
    - to build a particular version, but defaults to the most recent version.

    Args:
        args_list: This is a developer option used to simulate command line usage by
            providing a list of command line argument strings to the entry point
            function. For example, ``safedata_build_local_gbif -t 2012-01-01 file``
            can be replicated using ``_build_local_gbif_cli(['-t', '2012-01-01',
            'file'])``.

    Returns:
        An integer code showing success (0) or failure (1).
    """

    # If no arguments list is provided
    if args_list is None:
        args_list = sys.argv[1:]

    # Check function docstring exists to safeguard against -OO mode, and strip off the
    # description of the function args_list, which should not be included in the command
    # line docs
    if _build_local_gbif_cli.__doc__ is not None:
        desc = textwrap.dedent(
            "\n".join(_build_local_gbif_cli.__doc__.splitlines()[:-10])
        )
    else:
        desc = "Python in -OO mode: no docs"

    # create the parser
    parser = argparse.ArgumentParser(
        prog="safedata_build_local_gbif",
        description=desc,
        formatter_class=_desc_formatter,
    )

    parser.add_argument("outfile", help="Filename to use for database file.", type=Path)
    parser.add_argument(
        "-t",
        "--timestamp",
        default=None,
        type=str,
        help="The time stamp of a database archive version to use.",
    )

    args = parser.parse_args(args=args_list)

    # Validate output file
    # Look that file can be created without clobbering existing file
    outdir = args.outfile.absolute().parent

    if not (outdir.exists() and outdir.is_dir()):
        LOGGER.error("The output directory does not exist.")
        return 1

    if args.outfile.absolute().exists():
        LOGGER.error("The output file already exists.")
        return 1

    # Validate the timestamp
    try:
        timestamp, url = get_gbif_version(args.timestamp)
    except ValueError as excep:
        LOGGER.error(str(excep))
        return 1

    # Download and build from validated url
    with tempfile.TemporaryDirectory() as download_loc:
        file_data = download_gbif_backbone(
            outdir=download_loc, timestamp=timestamp, url=url
        )
        build_local_gbif(outfile=args.outfile, **file_data)

    return 0


def _build_local_ncbi_cli(args_list: list[str] | None = None) -> int:
    """Build a local NCBI database.

    This tool builds an SQLite database of the NCBI  taxonomy to use in
    validation by safedata_validate. There are multiple archived versions
    of the dataset, and the available versions can be seen here:

        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/

    The tool will optionally take a timestamp - using the format '2021-11-26'
    - to build a particular version, but defaults to the most recent version.

    Args:
        args_list: This is a developer option used to simulate command line usage by
            providing a list of command line argument strings to the entry point
            function. For example, ``safedata_build_local_ncbi -t 2012-01-01 file``
            can be replicated by calling ``_build_local_ncbi_cli(['-t', '2012-01-01',
            'file']))``.

    Returns:
        An integer code showing success (0) or failure (1).
    """

    # If no arguments list is provided
    if args_list is None:
        args_list = sys.argv[1:]

    # Check function docstring exists to safeguard against -OO mode, and strip off the
    # description of the function args_list, which should not be included in the command
    # line docs
    if _build_local_ncbi_cli.__doc__ is not None:
        desc = textwrap.dedent(
            "\n".join(_build_local_ncbi_cli.__doc__.splitlines()[:-10])
        )
    else:
        desc = "Python in -OO mode: no docs"

    # create the parser
    parser = argparse.ArgumentParser(
        prog="safedata_build_local_ncbi",
        description=desc,
        formatter_class=_desc_formatter,
    )

    parser.add_argument("outfile", help="Filename to use for database file.", type=Path)
    parser.add_argument(
        "-t",
        "--timestamp",
        default=None,
        type=str,
        help="The time stamp of a database archive version to use.",
    )

    args = parser.parse_args(args=args_list)

    # Validate output file
    # Look that file can be created without clobbering existing file
    outdir = args.outfile.absolute().parent

    if not (outdir.exists() and outdir.is_dir()):
        LOGGER.error("The output directory does not exist.")
        return 1

    if args.outfile.absolute().exists():
        LOGGER.error("The output file already exists.")
        return 1

    # Validate the timestamp
    try:
        filename, timestamp, filesize = get_ncbi_version(timestamp=args.timestamp)
    except ValueError as excep:
        LOGGER.error(str(excep))
        return 1

    # Download and build
    with tempfile.TemporaryDirectory() as download_loc:
        file_data = download_ncbi_taxonomy(
            outdir=download_loc,
            timestamp=timestamp,
            filename=filename,
            filesize=filesize,
        )
        build_local_ncbi(outfile=args.outfile, **file_data)

    return 0

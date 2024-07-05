"""This module provides functions to:

1. handle the publication of datasets after they have been validated
   using safedata_validate, including the generation of HTML descriptions
   of datasets.
2. maintain local copies of datasets in the folder structure expected
   by the safedata R package.
3. compile a RIS format bibliographic file for published datasets.
"""  # noqa D415

from __future__ import annotations

import decimal
import hashlib
import re
import shutil
from datetime import datetime as dt
from importlib import resources as il_resources  # avoid confusion with sdv.resources
from itertools import groupby
from pathlib import Path

import requests
import rispy
import simplejson
from dominate import tags
from jinja2 import Environment, FileSystemLoader, select_autoescape
from tqdm import tqdm
from tqdm.utils import CallbackIOWrapper

from safedata_validator.logger import FORMATTER, LOGGER
from safedata_validator.resources import Resources
from safedata_validator.taxa import taxon_index_to_text

# Constant definition of zenodo action function response type
# TODO - should all of these functions just Raise? Then we could just use try
# blocks rather than using the this approach.
ZenodoFunctionResponseType = tuple[dict, str | None]
"""Function return value

The functions interacting with Zenodo all return a common format of tuple of length 2:

* A dictionary containing the response content. For responses that do not generate a
  response content but just indicate success via HTTP status codes, an empty dictionary
  is returned. An empty dictionary is also returned when the function results in an
  error. For upload files, the response is a list of such dictionaries, one for each
  file uploaded.
* An error message on failure or None on success

So, for example:

```{python}
({'key': 'value'}, None)
({}, 'Something went wrong')
```

The expected use pattern is then:

```{python}
response, error = zenodo_function(args)
```
"""

ZenodoFunctionResponseListType = tuple[list[dict], str | None]
"""Function return value

This follows the ZenodoFunctionResponseType but returns a list of response dictionaries.
This is only used by the upload_files function, where each file upload has a response
containing useful information.
"""


# TODO: would make a lot more sense, particularly given the integrated
#       publish_dataset option to have functions take the output of
#       _resources_to_zenodo_api rather than continually recreating it inside each
#       function. Some functions need more resources, but most do not.


def _resources_to_zenodo_api(resources: Resources | None = None) -> dict:
    """Get a dictionary of the Zenodo and Metadata config from Resources.

    Args:
        resources: An instance of Resources or the default None to search for a local
            Resources configuration file.
    """

    # Get resource configuration
    if resources is None:
        resources = Resources()

    # Get the Zenodo API
    config_fail = False

    if resources.zenodo.use_sandbox is None:
        config_fail = True
    elif resources.zenodo.use_sandbox:
        zenodo_api = "https://sandbox.zenodo.org/api"
        token = resources.zenodo.zenodo_sandbox_token
    else:
        zenodo_api = "https://zenodo.org/api"
        token = resources.zenodo.zenodo_token

    if zenodo_api is None or token is None:
        config_fail = True

    # Get the contact details if used
    contact_name = resources.zenodo.contact_name
    contact_affiliation = None
    contact_orcid = None

    if contact_name is not None:
        contact_affiliation = resources.zenodo.contact_affiliation
        contact_orcid = resources.zenodo.contact_orcid

    # # Get the metadata api
    # metadata_api = resources.metadata.api
    # metadata_token = resources.metadata.token
    # if metadata_api is None or metadata_token is None:
    #     config_fail = True
    # metadata_ssl = resources.metadata.ssl_verify

    if config_fail:
        raise RuntimeError("safedata_validator not configured for Zenodo functions")

    return {
        "zapi": zenodo_api,
        "ztoken": {"access_token": token},
        "zcomm": resources.zenodo.community_name,
        "zcname": contact_name,
        "zcaffil": contact_affiliation,
        "zcorc": contact_orcid,
        # "mdapi": metadata_api,
        # "mdtoken": metadata_token,
        # "mdssl": metadata_ssl,
    }


def _compute_md5(fname: Path) -> str:
    """Calculate the md5 hash for a local file."""
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as fname_obj:
        for chunk in iter(lambda: fname_obj.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def _zenodo_error_message(response) -> str:
    """Format a Zenodo JSON error response as a string."""

    # Report the immediate reason and code along with the URL endpoint with the access
    # token redacted
    url = re.sub("(?<=access_token=).*$", "<redacted>", response.url)
    return_string = (
        f"\n\nZenodo error: {response.reason} ({response.status_code})\nURL: {url}\n"
    )

    # Attempt to get any JSON data
    try:
        response_json = response.json()
    except requests.exceptions.JSONDecodeError:
        return return_string

    return_string += f"Message: {response_json['message']}\n"

    # Report on error messages in response object
    errors = response_json.get("errors", [])
    if errors:
        return_string += "Errors:\n"
        for e in errors:
            messages = "\n    - ".join(e["messages"])
            return_string += f" * Messages for field {e['field']}:\n    - {messages}"
        return_string += "\n"

    return return_string


def _min_dp(val: float, min_digits: int = 2) -> str:
    """Display a number with a minimum number of decimal places.

    INSPIRE/GEMINI XML requires a bounding box with at least two decimal places
    provided, even if that is 4.50 rather than 4.5. This function tries to ensure that
    representation.
    """

    decimal_val = decimal.Decimal(str(val))
    round_precision = decimal_val.as_tuple().exponent
    format_precision = max(abs(round_precision), min_digits)  # type: ignore [arg-type]
    string_format = f"0.{format_precision}f"

    return format(val, string_format)


# Zenodo action functions


def get_deposit(
    deposit_id: int, resources: Resources | None = None
) -> ZenodoFunctionResponseType:
    """Download the metadata of a Zenodo deposit.

    Args:
        deposit_id: The Zenodo record id of an existing dataset.
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    zres = _resources_to_zenodo_api(resources)
    zenodo_api = zres["zapi"]
    params = zres["ztoken"]

    # request the deposit
    dep = requests.get(
        f"{zenodo_api}/deposit/depositions/{deposit_id}", params=params, json={}
    )

    # check for success and return the information.
    if dep.status_code == 200:
        return dep.json(), None
    else:
        return {}, _zenodo_error_message(dep)


def create_deposit(
    new_version: int | None = None,
    resources: Resources | None = None,
) -> ZenodoFunctionResponseType:
    """Create a new deposit.

    Creates a new deposit draft, possibly as a new version of an existing published
    record. Creating a new version requires the Zenodo ID of an existing dataset: this
    has to be the ID of the most recently published version of a dataset, not the
    concept ID used to group datasets or any of the older versions.

    Args:
        new_version: Optionally, create a new version of the dataset with the provided
            Zenodo ID.
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)
    zenodo_api = zres["zapi"]
    params = zres["ztoken"]

    # get the correct draft api
    if new_version is None:
        api = f"{zenodo_api}/deposit/depositions"
    else:
        api = f"{zenodo_api}/deposit/depositions/{new_version}/actions/newversion"

    # Create the draft
    new_draft = requests.post(api, params=params, json={})

    # trap errors in creating the new version (not 201: created)
    if new_draft.status_code != 201:
        return {}, _zenodo_error_message(new_draft)

    if new_version is None:
        return new_draft.json(), None

    # For new versions, the response is an update to the existing copy,
    # so need to separately retrieve the new draft
    api = new_draft.json()["links"]["latest_draft"]
    dep = requests.get(api, params=params, json={})

    # trap errors in creating the resource - successful creation of new version
    #  drafts returns 200
    if dep.status_code != 200:
        return {}, _zenodo_error_message(dep)
    else:
        return dep.json(), None


def upload_metadata(
    metadata: dict, zenodo: dict, resources: Resources | None = None
) -> ZenodoFunctionResponseType:
    """Upload dataset metadata.

    Takes a dictionary of dataset metadata, converts it to a JSON payload of Zenodo
    metadata and uploads it to a deposit.

    Args:
        metadata: The metadata dictionary for a dataset
        zenodo: The zenodo metadata dictionary for a deposit
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)

    # basic contents
    zen_md = {
        "metadata": {
            "upload_type": "dataset",
            # "publication_date": datetime.date.today().isoformat(),
            "title": metadata["title"],
            "keywords": metadata["keywords"],
            "license": "cc-by",
            "communities": [{"identifier": zres["zcomm"]}],
        }
    }

    # Add a contact name to contributors if provided in config
    if zres["zcname"] is not None:
        zen_md["metadata"]["contributors"] = [
            {
                "name": zres["zcname"],
                "type": "ContactPerson",
                "affiliation": zres["zcaffil"],
                "orcid": zres["zcorc"],
            }
        ]

    # set up the access rights
    dataset_access = metadata["access"].lower()
    if dataset_access == "embargo":
        zen_md["metadata"]["access_right"] = "embargoed"
        zen_md["metadata"]["embargo_date"] = metadata["embargo_date"]
    elif dataset_access == "open":
        zen_md["metadata"]["access_right"] = "open"
    elif dataset_access == "restricted":
        zen_md["metadata"]["access_right"] = "restricted"
        zen_md["metadata"]["access_conditions"] = metadata["access_conditions"]
    else:
        raise ValueError("Unknown access status")

    # set up the dataset creators - the format has already been checked and names
    # should be present and correct. Everything else is optional, so strip None
    # values and pass the rest to Zenodo
    zen_md["metadata"]["creators"] = [
        {ky: auth[ky] for ky in auth if auth[ky] is not None and ky != "email"}
        for auth in metadata["authors"]
    ]

    # Add the html description
    zen_md["metadata"]["description"] = dataset_description(
        dataset_metadata=metadata, resources=resources
    )

    # attach the metadata to the deposit resource
    mtd = requests.put(zenodo["links"]["self"], params=zres["ztoken"], json=zen_md)

    # trap errors in uploading metadata and tidy up
    if mtd.status_code != 200:
        return {}, mtd.reason
    else:
        return {}, None


def update_published_metadata(
    zenodo: dict,
    resources: Resources | None = None,
) -> ZenodoFunctionResponseType:
    """Update published deposit metadata.

    Updates the metadata on a published deposit, for example to modify the access status
    of deposit. In general, metadata should be updated by releasing a new version of the
    dataset, and this function should only be used where it is essential that the
    published version by altered.

    Args:
        zenodo: A Zenodo metadata dictionary, with an updated metadata section
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)

    links = zenodo["links"]

    # Unlock the published deposit for editing
    edt = requests.post(links["edit"], params=zres["ztoken"])

    if edt.status_code != 201:
        return {}, edt.json()

    # # Amend the metadata
    # for key, val in new_values.items():
    #     if val is not None:
    #         metadata[key] = val
    #     elif key in metadata:
    #         metadata.pop(key)

    # If any API calls from now fail, we need to tidy up the edit
    # status of the record, or it will block subsequent attempts

    upd = requests.put(
        links["self"],
        params=zres["ztoken"],
        headers={"Content-Type": "application/json"},
        data=simplejson.dumps({"metadata": zenodo["metadata"]}),
    )

    success_so_far = 0 if upd.status_code != 200 else 1
    ret = upd.json()

    # Republish to save the changes
    if success_so_far:
        pub = requests.post(links["publish"], params=zres["ztoken"])
        success_so_far = 0 if pub.status_code != 202 else 1
        ret = pub.json()

    # If all steps have been successful, return a 0 code, otherwise
    # try to discard the edits and return the most recent failure
    # notice

    if success_so_far:
        return ret, None
    else:
        dsc = requests.post(links["discard"], params=zres["ztoken"])
        success_so_far = 0 if dsc.status_code != 201 else 1
        if not success_so_far:
            ret = dsc.json()

        return {}, ret


def upload_files(
    metadata: dict,
    filepaths: list[Path],
    progress_bar: bool = True,
    resources: Resources | None = None,
) -> ZenodoFunctionResponseListType:
    """Upload a file to Zenodo.

    Uploads the a list of files to an unpublished Zenodo deposit. If any filenames
    already exists in the deposit, they will be replaced with the new content

    Args:
        metadata: The Zenodo metadata dictionary for a deposit
        filepaths: The path to the file to be uploaded
        progress_bar: Should the upload progress be displayed
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)
    params = zres["ztoken"]

    # Resolve filepaths and check they are all existing files
    filepaths = [f.resolve() for f in filepaths]
    bad_paths = [str(f) for f in filepaths if not (f.exists() and f.is_file())]

    if bad_paths:
        raise OSError(f"Filepaths unknown or not a file: {','.join(bad_paths)} ")

    # Collect response for each file
    response_content = []

    # Upload each file
    for fpath in filepaths:
        # upload the file
        # - https://gist.github.com/tyhoff/b757e6af83c1fd2b7b83057adf02c139
        file_size = fpath.stat().st_size
        api = f"{metadata['links']['bucket']}/{fpath.name}"

        with open(fpath, "rb") as file_io:
            print(f"Uploading {fpath.name}")
            if progress_bar:
                with tqdm(
                    total=file_size, unit="B", unit_scale=True, unit_divisor=1024
                ) as upload_monitor:
                    # Upload the wrapped file
                    wrapped_file = CallbackIOWrapper(
                        upload_monitor.update, file_io, "read"
                    )
                    fls = requests.put(api, data=wrapped_file, params=params)
            else:
                fls = requests.put(api, data=file_io, params=params)

        # trap errors in uploading file
        # - no success or mismatch in md5 checksums
        if fls.status_code != 201:
            return list({}), _zenodo_error_message(fls)

        # TODO - could this be inside the tqdm with call?
        #      - both are looping over the file contents
        # https://medium.com/codex/chunked-uploads-with-binary-files-in-python-f0c48e373a91
        local_hash = _compute_md5(fpath)

        fls_json: dict = fls.json()
        if fls_json["checksum"] != f"md5:{local_hash}":
            return list({}), "Mismatch in local and uploaded MD5 hashes"

        response_content.append(fls_json)

    return response_content, None


def discard_deposit(
    zenodo_metadata: dict, resources: Resources | None = None
) -> ZenodoFunctionResponseType:
    """Discard a deposit.

    Deposits can be discarded - the associated files and metadata will be deleted and
    the Zenodo ID no longer exists. Once deposits are published to records, they cannot
    be deleted via the API - contact the Zenodo team for help.

    Args:
        zenodo_metadata: The Zenodo metadata dictionary for a deposit
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)
    params = zres["ztoken"]

    delete = requests.delete(zenodo_metadata["links"]["self"], params=params)

    if delete.status_code == 204:
        return {"result": "success"}, None
    else:
        return {}, _zenodo_error_message(delete)


def publish_deposit(
    zenodo: dict, resources: Resources | None = None
) -> ZenodoFunctionResponseType:
    """Publish a created deposit.

    Args:
        zenodo: The dataset metadata dictionary for a deposit
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)
    params = zres["ztoken"]

    # publish
    pub = requests.post(zenodo["links"]["publish"], params=params)

    # trap errors in publishing, otherwise return the publication metadata
    if pub.status_code != 202:
        return {}, pub.json()
    else:
        return pub.json(), None


def delete_files(
    metadata: dict, filenames: list[str], resources: Resources | None = None
) -> ZenodoFunctionResponseType:
    """Delete an uploaded file from an unpublished Zenodo deposit.

    Args:
        metadata: The Zenodo metadata dictionary for a deposit
        filenames: A list of files to delete from the deposit
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)
    params = zres["ztoken"]

    # get an up to date list of existing files (metadata
    # might be outdated)
    files = requests.get(metadata["links"]["files"], params=params)

    # check the result of the files request
    if files.status_code != 200:
        # failed to get the files
        return {}, _zenodo_error_message(files)

    # Get a dictionary of the available file links
    files_dict = {f["filename"]: f["links"]["self"] for f in files.json()}

    # Get matching files
    unknown_files = []
    delete_links = []
    for file in filenames:
        if file in files_dict:
            delete_links.append(files_dict[file])
        else:
            unknown_files.append(file)

    if unknown_files:
        return {}, f"Files not found in the deposit: {','.join(unknown_files)}"

    return _delete_files_from_links(delete_links=delete_links, params=params)


def _delete_files_from_links(delete_links: list[str], params: dict):
    """Delete files from Zenodo deposit from a list of API links.

    This private method is used to delete files from a deposit using a list of the
    Zenodo file "self" links. This is following the API below, where the link text
    provides the specific API link to the file to be deleted

        DELETE /api/deposit/depositions/:id/files/:file_id

    Args:
        delete_links: A list of Zenodo deposit file 'self' links
        params: A dictionary of authentication parameters for the Zenodo API
    """

    for link in delete_links:
        file_del = requests.delete(link, params=params)

        if file_del.status_code != 204:
            return {}, _zenodo_error_message(file_del)

    return {"result": "success"}, None


"""
Dataset description generation (HTML and GEMINI XML)
"""


def dataset_description(
    dataset_metadata: dict,
    resources: Resources | None = None,
) -> tags.div | str:
    """Create an HTML dataset description.

    This function takes the dataset metadata exported by safedata_validate and uses it
    to populate an HTML template file. The resulting HTML can then be used to to provide
    a summary description of the dataset, either for local use or to upload as the
    description component of the Zenodo metadata,

    A default template is provided with the safedata_validator package, but users can
    provide bespoke templates via the configuration file.

    Args:
        dataset_metadata: The dataset metadata
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        Either a string of rendered HTML or a dominate.tags.div object.
    """

    # NOTE - this could conceivably just pass the complete dataset metadata dictionary
    #        straight to the Jinja template context. That would make the template more
    #        complex but would also expose all of the metadata.

    # Load resources if needed
    if resources is None:
        resources = Resources()

    # Get the template path elements
    if resources.zenodo.html_template is None:
        template_path = il_resources.files("safedata_validator.templates").joinpath(
            "description_template.html"
        )
    else:
        template_path = Path(resources.zenodo.html_template)
        if not template_path.exists():
            raise FileNotFoundError(
                f"Configured html template not found: {resources.zenodo.html_template}"
            )

    # Using autoescape=False is not generally recommended, but the title and taxa
    # context elements contain HTML tags
    # - mypy: importlib returns a Traversable, which is a protocol that Path complies
    #         with, but the attribute isn't being recognized
    env = Environment(
        loader=FileSystemLoader(template_path.parent),  # type: ignore [attr-defined]
        autoescape=False,
    )

    template = env.get_template(template_path.name)

    # PROJECT Title and authors are added by Zenodo from zenodo metadata
    # TODO - option to include here?

    context_dict = dict(
        description=dataset_metadata["description"].replace("\n", "</br>")
    )

    # proj_url = URL('projects', 'project_view', args=[metadata['project_id']],
    #               scheme=True, host=True)
    # desc += P(B('Project: '), 'This dataset was collected as part of the following '
    #                          'SAFE research project: ', A(B(title), _href=proj_url))
    ##

    # Funding information
    context_dict["funders"] = dataset_metadata["funders"]
    context_dict["permits"] = dataset_metadata["permits"]

    # Filenames associated with the dataset
    context_dict["dataset_filename"] = dataset_metadata["filename"]
    # TODO - the external file default in the metadata definition should be an
    #        empty list, not None
    context_dict["external_files"] = (
        []
        if dataset_metadata["external_files"] is None
        else dataset_metadata["external_files"]
    )

    context_dict["all_filenames"] = [context_dict["dataset_filename"]] + [
        f["file"] for f in context_dict["external_files"]
    ]

    # Group the sheets by their 'external' file - which is None for sheets in the
    # submitted workbook - and collect them into a dictionary by source file. Because
    # you can't sort a mix of strings and None elements, this substitutes in
    # '__internal__' to represent internal sheets.
    tables_by_source = [
        (sh["external"] or "__internal__", sh)
        for sh in dataset_metadata["dataworksheets"]
    ]

    # Now group into a dictionary keyed by __internal__ or external file names
    tables_by_source.sort(key=lambda sh: sh[0])
    tables_grouped_by_source = groupby(tables_by_source, key=lambda sh: sh[0])

    # Convert to a list of table information, keyed by file.
    tables_dict_by_source = {
        ky: [val[1] for val in tpl] for ky, tpl in tables_grouped_by_source
    }

    # We've now  a dictionary of table descriptions that might have an entry for each
    # provided file. Get the internal tables separately in the context
    if "__internal__" in tables_dict_by_source:
        context_dict["internal_tables"] = tables_dict_by_source.pop("__internal__")
    else:
        context_dict["internal_tables"] = []

    # Now need to pair any external table metadata with the external file descriptions.
    # TODO - the external file default in the metadata definition should be an
    #        empty list, not None
    if dataset_metadata["external_files"] is None:
        external_files = dict()
    else:
        # Repackage external metadata to be keyed by file name and provide description
        # and a default empty list of tables
        external_files = {
            vl["file"]: {"description": vl["description"], "tables": []}
            for vl in dataset_metadata["external_files"]
        }
        # Add the remaining table descriptions to the appropriate files.
        for extf_key, extf_tabs in tables_dict_by_source.items():
            external_files[extf_key]["tables"] = extf_tabs

    context_dict["external_file_data"] = external_files

    # Populate a list of filenames
    context_dict["all_filenames"] = [
        context_dict["dataset_filename"],
        *list(external_files.keys()),
    ]

    # Add extents if populated
    context_dict["temporal_extent"] = dataset_metadata["temporal_extent"]
    context_dict["latitudinal_extent"] = dataset_metadata["latitudinal_extent"]
    context_dict["longitudinal_extent"] = dataset_metadata["longitudinal_extent"]

    # Find taxa data from each database and convert to HTML representation. The metadata
    # will be an empty list if the dataset does not contain any taxa.
    context_dict["gbif_timestamp"] = dataset_metadata["gbif_timestamp"]
    context_dict["ncbi_timestamp"] = dataset_metadata["ncbi_timestamp"]

    gbif_taxon_index = dataset_metadata["gbif_taxa"]
    ncbi_taxon_index = dataset_metadata["ncbi_taxa"]

    context_dict["gbif_taxa"] = (
        taxon_index_to_text(taxa=gbif_taxon_index, html=True, auth="GBIF")
        if gbif_taxon_index
        else None
    )

    context_dict["ncbi_taxa"] = (
        taxon_index_to_text(taxa=ncbi_taxon_index, html=True, auth="NCBI")
        if ncbi_taxon_index
        else None
    )

    html = template.render(context_dict)

    return html


def generate_inspire_xml(
    dataset_metadata: dict,
    zenodo_metadata: dict,
    resources: Resources,
    lineage_statement: str | None = None,
) -> str:
    """Convert dataset and zenodo metadata into GEMINI XML.

    Produces an INSPIRE/GEMINI formatted XML record from dataset metadata,
    and Zenodo record metadata using a template XML file. The dataset URL
    defaults to the Zenodo record but can be replaced if a separate URL (such as
    a project specific website) is used. The Gemini XML standard requires a
    statement about the lineage of a dataset - this is automatically taken from the
    package configuration but can be overridden for individual datasets, for example to
    add dataset specific links, using the `lineage_statement` argument.

    Args:
        dataset_metadata: A dictionary of the dataset metadata
        zenodo_metadata: A dictionary of the Zenodo record metadata
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.
        lineage_statement: An optional alternative lineage statement about the data.

    Returns:
        A string containing GEMINI compliant XML.
    """

    template_path = il_resources.files("safedata_validator.templates").joinpath(
        "gemini_xml_template.xml"
    )

    # Get the Jinja environment and load the template
    # - mypy: importlib returns a Traversable, which is a protocol that Path complies
    #         with, but the attribute isn't being recognized
    env = Environment(
        loader=FileSystemLoader(template_path.parent),  # type: ignore [attr-defined]
        autoescape=select_autoescape(),
    )

    template = env.get_template(template_path.name)

    # Build some reused values from the metadata
    # URIs -  form the DOI URL from the prereserved DOI metadata
    doi_url = f"https://doi.org/{zenodo_metadata['metadata']['prereserve_doi']['doi']}"

    # A true "publication" date is not available until a record is published, so use the
    # creation date of the deposit as a reasonable replacement, with the caveat that you
    # should generate the XML and publish on the same day.
    pub_date = dt.fromisoformat(zenodo_metadata["created"]).date()

    # A citation string
    authors = [au["name"] for au in dataset_metadata["authors"]]
    author_string = ", ".join(authors)
    if len(authors) > 1:
        author_string = author_string.replace(", " + authors[-1], " & " + authors[-1])

    citation_string = (
        f"{author_string} ({pub_date.year}) "
        f"{dataset_metadata['title']} [Dataset] {doi_url}"
    )

    # Resource constraints text
    if dataset_metadata["access"] == "embargo":
        access_statement = (
            f"This data is under embargo until {dataset_metadata['embargo_date']}."
            "After that date there are no restrictions to public access."
        )
    elif dataset_metadata["access"] == "restricted":
        access_statement = (
            "This dataset is currently not publicly available, please contact the "
            "Zenodo community owner to request access."
        )
    else:
        access_statement = "There are no restrictions to public access."

    # Get a copy of the project wide XML configuration from the resources and update it
    # with the file specific elements from the zenodo and dataset metadata
    context_dict = resources.xml.copy()

    context_dict.update(
        # Values also used on the Zenodo information or duplicated in the xml
        contactName=resources.zenodo.contact_name,
        contactOrcID=resources.zenodo.contact_orcid,
        pointofcontactName=resources.zenodo.contact_name,
        pointofcontactCountry=resources.xml.contactCountry,
        pointofcontactEmail=resources.xml.contactEmail,
        pointofcontactOrcID=resources.zenodo.contact_name,
        # Dataset specific information
        citationRSIdentifier=doi_url,
        dateStamp=pub_date.isoformat(),
        publicationDate=pub_date.isoformat(),
        fileIdentifier=str(zenodo_metadata["id"]),
        title=dataset_metadata["title"],
        authors=dataset_metadata["authors"],
        abstract=dataset_metadata["description"],
        keywords=dataset_metadata["keywords"],
        citationString=citation_string,
        embargoValue=access_statement,
        startDate=dataset_metadata["temporal_extent"][0][:10],
        endDate=dataset_metadata["temporal_extent"][1][:10],
        westBoundLongitude=_min_dp(dataset_metadata["longitudinal_extent"][0], 2),
        eastBoundLongitude=_min_dp(dataset_metadata["longitudinal_extent"][1], 2),
        southBoundLatitude=_min_dp(dataset_metadata["latitudinal_extent"][0], 2),
        northBoundLatitude=_min_dp(dataset_metadata["latitudinal_extent"][1], 2),
        downloadLink=doi_url,
    )

    # Override global lineage statement
    if lineage_statement is not None:
        context_dict["lineageStatement"] = lineage_statement

    xml = template.render(context_dict)

    return xml


def publish_dataset(
    resources: Resources,
    dataset: Path,
    dataset_metadata: dict,
    external_files: list[Path],
    new_version: int | None,
    no_xml: bool = False,
) -> tuple[int, str]:
    """Publish a validated dataset.

    This function takes a dataset and its validated metadata, along with any additional
    files named in the dataset, and publishes them to a new Zenodo record. It merges
    several of the :mod:`~safedata_validator.zenodo` functions to provide a single
    interface to carry out the complete publication process. The function checks that
    the set of provided files (dataset and external files) matches the files documented
    in the dataset.

    When a new version of an existing database is created, the resulting Zenodo deposit
    contains copies of the most recent files. At present, the publish dataset command
    deletes all of these files and expects to upload the full set of replacement files.
    This will be inefficient if only some files need to be changed, and this function
    may be updated in the future to allow only new files to be updated.

    It returns the URL of the resulting published datset. If the publication process
    fails, the partly completed deposit is deleted to avoid cluttering the Zenodo
    deposit list.

    Args:
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.
        dataset: A path to the dataset file.
        external_files: A list of paths to external files named in the dataset.
        dataset_metadata: The dataset metadata.
        new_version: Optionally, create a new version of the dataset with the provided
            Zenodo ID. This must be the most recent version of the dataset concept.
        no_xml: A flag to suppress the automatic inclusion of Gemini XML metadata.

    Returns:
        A tuple containing the id number and URL of the new record.

    Raises:
        FileNotFoundError: dataset or external files not found.
        ValueError: provided files do not match documented files.
        RuntimeError: issues with Zenodo API.
    """

    # Check the files to upload exist.
    paths_to_upload = [dataset, *external_files]
    not_found = [str(f) for f in paths_to_upload if not f.exists()]

    if not_found:
        raise FileNotFoundError(f"Files not found: {', '.join(not_found)}")

    # Check for duplicated filenames
    filenames_only = [f.name for f in paths_to_upload]
    if len(filenames_only) != len(set(filenames_only)):
        raise ValueError("Files to be uploaded do not have unique names.")

    # Are all external files listed in the dataset provided
    if dataset_metadata["external_files"] is None:
        metadata_ext_files = set()
    else:
        metadata_ext_files = {val["file"] for val in dataset_metadata["external_files"]}
    provided_ext_files = {val.name for val in external_files}

    if metadata_ext_files != provided_ext_files:
        raise ValueError(
            "External file names in dataset do not match provided "
            f"external file names: {', '.join(metadata_ext_files)}"
        )

    # For new versions of an existing dataset, get the existing dataset metadata and
    # figure out which files are being changed.
    if new_version is not None:
        # Get the requested version data
        requested_version_metadata, error = get_deposit(
            deposit_id=new_version, resources=resources
        )
        if error:
            # Report on the failure to get the existing version data
            raise RuntimeError(error)

        # Check if that is the latest version, because that is what is _actually_ cloned
        # when a version is requested and we need to get the most recent file listing to
        # update the files sanely.
        latest_version_metadata = requests.get(
            requested_version_metadata["links"]["latest"]
        ).json()

        latest_id = latest_version_metadata["id"]

        if latest_id != new_version:
            raise RuntimeError(
                f"Provided version id ({new_version}) is not the "
                f"most recent version id({latest_id})"
            )

        # Collect the name and md5 sum of the existing files and the incoming files to
        # allow for checking of files with identical name and content
        incoming_files = {(p.name, _compute_md5(p)) for p in paths_to_upload}
        existing_files = {
            (p["key"], p["checksum"].removeprefix("md5:"))
            for p in latest_version_metadata["files"]
        }

        # Split files into files to be upload, files to be deleted from deposit and
        # files that are not being changed
        new_or_updated_files = [val[0] for val in (incoming_files - existing_files)]
        files_to_remove = [val[0] for val in (existing_files - incoming_files)]
        unchanged = [val[0] for val in (incoming_files & existing_files)]

        # Trap updates with no differences - no new files or updates and only a
        # GEMINI.xml file to remove
        no_changes = (not new_or_updated_files) or (
            len(files_to_remove) == 1 and files_to_remove[0].endswith("_GEMINI.xml")
        )
        if no_changes:
            raise RuntimeError(
                "No file changes: content identical to existing version."
            )

        # Reduce the upload paths to only the paths of new or updated files. Note the
        # earlier code checks that the filenames are unique.
        paths_to_upload = [p for p in paths_to_upload if p.name in new_or_updated_files]

        # Report on update plan.
        print(f"Preparing new version of deposit {new_version}")
        if unchanged:
            print(f" - Unchanged files: {', '.join(unchanged)}")
        if files_to_remove:
            print(f" - Removing outdated files: {', '.join(files_to_remove)}")
        if new_or_updated_files:
            print(
                f" - Uploading new or updated files: {', '.join(new_or_updated_files)}"
            )

    else:
        # No files will need to be removed from completely new datasets.
        files_to_remove = []

    # From here on, we need to trap failures and tidy up, so define an local function to
    # handle the bail process
    def _clean_up_and_bail(zenodo_metadata, resources) -> None:
        """Handle the tidy up if the publication process fails."""
        publish_response, error = discard_deposit(
            zenodo_metadata=zenodo_metadata, resources=resources
        )
        print("Issue with publication process - draft deposit discarded.")
        raise RuntimeError(error)

    # Create the new deposit to publish the dataset
    zenodo_metadata, error = create_deposit(
        resources=resources, new_version=new_version
    )

    if error:
        _clean_up_and_bail(zenodo_metadata=zenodo_metadata, resources=resources)

    # Report success
    zenodo_id = zenodo_metadata["id"]
    print(f"Deposit created: {zenodo_id}")

    # Generate XML if requested
    if not no_xml:
        xml_content = generate_inspire_xml(
            dataset_metadata=dataset_metadata,
            zenodo_metadata=zenodo_metadata,
            resources=resources,
        )

        xml_file = dataset.parent / f"{zenodo_id}_GEMINI.xml"
        with open(xml_file, "w") as xml_out:
            xml_out.write(xml_content)
        paths_to_upload.append(xml_file)
        print(f"XML created: {xml_file}")

    # Remove any outdated files from deposits created as new versions
    if files_to_remove:
        # Get the authentication parameters - see TODO at top about passing zres around
        zres = _resources_to_zenodo_api(resources)
        params = zres["ztoken"]

        # Get all of the file links returned by create_deposit and then delete them.
        print(f"Removing outdated files: {','.join(files_to_remove)}")
        removal_links = [
            f["links"]["self"]
            for f in zenodo_metadata["files"]
            if f["filename"] in files_to_remove
        ]

        delete_result, error = _delete_files_from_links(
            delete_links=removal_links, params=params
        )

        # Handle errors
        if error:
            _clean_up_and_bail(zenodo_metadata=zenodo_metadata, resources=resources)

    # Post the files
    print("Uploading files:")
    file_upload_response, error = upload_files(
        metadata=zenodo_metadata, filepaths=paths_to_upload, resources=resources
    )
    # Handle errors
    if error:
        _clean_up_and_bail(zenodo_metadata=zenodo_metadata, resources=resources)

    # Post the metadata
    print("Uploading deposit metadata")
    md_upload_response, error = upload_metadata(
        metadata=dataset_metadata, zenodo=zenodo_metadata, resources=resources
    )
    if error:
        _clean_up_and_bail(zenodo_metadata=zenodo_metadata, resources=resources)

    # Publish the deposit
    publish_response, error = publish_deposit(
        zenodo=zenodo_metadata, resources=resources
    )
    if error:
        _clean_up_and_bail(zenodo_metadata=zenodo_metadata, resources=resources)

    # Return the new publication ID and link
    zenodo_url = publish_response["links"]["html"]
    print(f"Dataset published: {zenodo_url}")

    return (zenodo_id, zenodo_url)


# Bibliographic and local database functions


def download_ris_data(
    resources: Resources | None = None, ris_file: Path | None = None
) -> None:
    """Downloads Zenodo records into a RIS format bibliography file.

    This function is used to maintain a bibliography file of the records
    uploaded to a safedata community on Zenodo. It accesses the Zenodo community
    specified in the resource configuration and downloads all records. It then
    optionally checks the list of downloaded DOIs against the content of an
    existing RIS file and then downloads citations for all new DOIs from
    datacite.org.

    Args:
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.
        ris_file: The path to an existing RIS format file containing previously
            downloaded records.

    Returns:
        A list of strings containing RIS formatted citation data.
    """

    if resources is None:
        resources = Resources()

    # Get a list of known DOI records from an existing RIS file if one is
    # provided
    known_recids = []
    new_doi = []

    if ris_file is not None and ris_file.exists():
        with open(ris_file) as bibliography_file:
            entries = rispy.load(bibliography_file)
            for entry in entries:
                record_id = int(entry["url"].split("/")[-1])
                known_recids.append(record_id)

    # Zenodo API call to return the records associated with the SAFE community
    zres = _resources_to_zenodo_api(resources)
    z_api = zres["zapi"]
    z_cname = zres["zcomm"]

    api = f"{z_api}/records/?q=communities:{z_cname}"

    # Provide feedback on DOI collection
    LOGGER.info(f"Fetching record DOIs from {api}:")
    FORMATTER.push()

    # The API is paged - it contains a set of records and a link that points
    # to the next page of records, so keep looping until there are no more next
    n_records = 0
    while True:
        # Get the data
        safe_data = requests.get(api)

        if safe_data.status_code != 200:
            raise OSError("Cannot access Zenodo API")
        else:
            # Retrieve the record data and store the DOI for each record
            safe_data_dict = safe_data.json()
            for hit in safe_data_dict["hits"]["hits"]:
                if hit["id"] not in known_recids:
                    new_doi.append(hit["doi"])

            # Reporting
            n_records += len(safe_data_dict["hits"]["hits"])
            LOGGER.info(f"{n_records}")

            # Update the link for the next page, unless there is no next page
            if "next" in safe_data_dict["links"]:
                api = safe_data_dict["links"]["next"]
            else:
                break

    # Use the datacite API to retrieve the citation data associated with the DOI
    # and save it out to a RIS format file
    if not new_doi:
        LOGGER.info("No new DOIs found")
        return

    # Get the DOI data
    data = []

    FORMATTER.pop()
    LOGGER.info(
        f"Retrieving citation data from Datacite for {len(new_doi)} new records"
    )
    FORMATTER.push()

    for doi in new_doi:
        ris_data = requests.get(
            f"https://data.datacite.org/application/x-research-info-systems/{doi}"
        )

        if ris_data.status_code != 200:
            LOGGER.warning(f"DOI {doi} not found in datacite.org")
        else:
            # Write the response content to the data list. It comes in as byte
            # data so needs to be decoded to a string variable
            LOGGER.info(f"Retrieved citation for DOI {doi}")
            data.append(ris_data.content.decode("utf-8") + "\r\n")

    FORMATTER.pop()

    # Writing only occurs if a ris file path has actually been provided
    if ris_file:
        if ris_file.exists():
            LOGGER.info(f"Appending RIS data for {len(data)} new records to {ris_file}")
            write_mode = "a"
        else:
            LOGGER.info(f"Writing RIS data for {len(data)} records to {ris_file}")
            write_mode = "w"

        with open(ris_file, write_mode) as ris_file_out:
            for this_entry in data:
                ris_file_out.write(this_entry)


def sync_local_dir(
    datadir: Path,
    xlsx_only: bool = True,
    replace_modified: bool = False,
    resources: Resources | None = None,
) -> None:
    """Synchronise a local data directory with a Zenodo community.

    The safedata R package defines a directory structure used to store metadata and
    files downloaded from a safedata community on Zenodo and from a safedata metadata
    server. This tool allows a safedata developer or community maintainer to create or
    update such a directory with _all_ of the resources in the Zenodo community,
    regardless of their public access status. This forms a backup (although Zenodo is
    heavily backed up) but also provides local copies of the files for testing and
    development of the code packages.

    This function requires that the resources are configured with access tokens for
    Zenodo and the details of the metadata server.

    Args:
        datadir: The path to a local directory containing an existing safedata
            directory or an empty folder in which to create one.
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.
        xlsx_only: Should the download ignore large non-xlsx files, defaulting
            to True.
        replace_modified: Should the synchronisation replace locally modified files with
            the archived version. By default, modified local files are left alone.
    """

    # Private helper functions
    def _get_file(url: str, outf: Path, params: dict | None = None) -> None:
        """Download a file from a URL."""
        resource = requests.get(url, params=params, stream=True)

        with open(outf, "wb") as outf_obj:
            shutil.copyfileobj(resource.raw, outf_obj)

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)
    zenodo_api = zres["zapi"]
    params = zres["ztoken"]

    # The dir argument should be an existing path
    if not (datadir.exists() and datadir.is_dir()):
        raise OSError(f"{datadir} is not an existing directory")

    # Get the configured metadata api
    api = zres["mdapi"]

    # Check for an existing API url file and check it is congruent with config
    url_file = datadir / "url.json"

    if url_file.exists():
        with open(url_file) as urlf:
            dir_api = simplejson.load(urlf)["url"][0]

        if api != dir_api:
            raise RuntimeError(
                "Configured api does not match existing api in directory"
            )
    else:
        with open(url_file, "w") as urlf:
            simplejson.dump({"url": [api]}, urlf)

    # Download index files - don't bother to check for updates, this isn't
    # a frequent thing to do
    LOGGER.info("Downloading index files")
    _get_file(f"{api}/api/index", datadir / "index.json")
    _get_file(f"{api}/api/gazetteer", datadir / "gazetteer.geojson")
    _get_file(f"{api}/api/location_aliases", datadir / "location_aliases.csv")

    # Get the deposits associated with the account, which includes a list of download
    # links
    params["page"] = 1
    deposits = []

    LOGGER.info("Scanning Zenodo deposits")
    while True:
        this_page = requests.get(
            zenodo_api + "/deposit/depositions",
            params=params,
            json={},
            headers={"Content-Type": "application/json"},
        )

        if not this_page.ok:
            raise RuntimeError("Could not connect to Zenodo API. Invalid token?")

        if this_page.json():
            deposits += this_page.json()
            LOGGER.info(f"Page {params['page']}")
            params["page"] += 1
        else:
            break

    LOGGER.info(f"Processing {len(deposits)} deposits")

    # Download the files
    for dep in deposits:
        con_rec_id = str(dep["conceptrecid"])
        rec_id = str(dep["record_id"])

        if not dep["submitted"]:
            LOGGER.info(f"Unsubmitted draft {con_rec_id}/{rec_id}")
            continue

        LOGGER.info(f"Processing deposit {con_rec_id}/{rec_id}")
        FORMATTER.push()

        # Create the directory structure if needed
        rec_dir = datadir / con_rec_id / rec_id
        if not rec_dir.exists():
            LOGGER.info("Creating directory")
            rec_dir.mkdir()
        else:
            LOGGER.info("Directory found")

        # loop over the files in the record
        for this_file in dep["files"]:
            if xlsx_only and not this_file["filename"].endswith(".xlsx"):
                LOGGER.info(f"Skipping non-excel file {this_file['filename']}")
                continue

            LOGGER.info(f"Processing {this_file['filename']}")
            FORMATTER.push()

            outf = rec_dir / this_file["filename"]
            local_copy = outf.exists()

            if not local_copy:
                LOGGER.info("Downloading")
                _get_file(this_file["links"]["download"], outf, params=params)
            elif local_copy and _compute_md5(outf) != this_file["checksum"]:
                if replace_modified:
                    LOGGER.info("Replacing locally modified file")
                    _get_file(this_file["links"]["download"], outf, params=params)
                else:
                    LOGGER.warning("Local copy modified")
            else:
                LOGGER.info("Already present")

            FORMATTER.pop()

        # Get the metadata json
        metadata = rec_dir / f"{rec_id}.json"
        if metadata.exists():
            LOGGER.info("JSON Metadata found")
        else:
            LOGGER.info("Downloading JSON metadata ")
            _get_file(f"{api}/api/record/{rec_id}", metadata)

        FORMATTER.pop()

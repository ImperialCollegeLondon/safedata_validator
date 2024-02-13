"""This module provides functions to:

1. handle the publication of datasets after they have been validated
   using safedata_validate, including the generation of HTML descriptions
   of datasets.
2. maintain local copies of datasets in the folder structure expected
   by the safedata R package.
3. compile a RIS format bibliographic file for published datasets.
"""  # noqa D415

import copy
import hashlib
import os
import shutil
from datetime import datetime as dt
from itertools import groupby
from pathlib import Path
from typing import Optional, Union

import requests  # type: ignore
import rispy
import simplejson
from dominate import tags
from dominate.util import raw
from lxml import etree
from tqdm import tqdm
from tqdm.utils import CallbackIOWrapper

from safedata_validator.logger import FORMATTER, LOGGER
from safedata_validator.resources import Resources
from safedata_validator.taxa import taxon_index_to_text

# Constant definition of zenodo action function response type
ZenodoFunctionResponseType = tuple[dict, Optional[str]]
"""Function return value

The functions interacting with Zenodo all return a common format of tuple of length 2:

* A dictionary containing the response content. For responses that do not generate a
  response content but just indicate success via HTTP status codes, an empty dictionary
  is returned. An empty dictionary is also returned when the function results in an
  error.
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


def _resources_to_zenodo_api(resources: Optional[Resources] = None) -> dict:
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
        zenodo_api = resources.zenodo.zenodo_sandbox_api
        token = resources.zenodo.zenodo_sandbox_token
    else:
        zenodo_api = resources.zenodo.zenodo_api
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


def _compute_md5(fname: str) -> str:
    """Calculate the md5 hash for a local file."""
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as fname_obj:
        for chunk in iter(lambda: fname_obj.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def _zenodo_error_message(response) -> str:
    """Format a Zenodo JSON error response as a string."""
    return f"{response.json()['message']} ({response.json()['status']})"


# Zenodo action functions


def get_deposit(
    deposit_id: int, resources: Optional[Resources] = None
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
    concept_id: Optional[int] = None, resources: Optional[Resources] = None
) -> ZenodoFunctionResponseType:
    """Create a new deposit.

    Creates a new deposit draft, possibly as a new version of an existing published
    record.

    Args:
        concept_id: An optional concept id of a published record to create a new version
            of an existing dataset.
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
    if concept_id is None:
        api = f"{zenodo_api}/deposit/depositions"
    else:
        api = f"{zenodo_api}/deposit/depositions/{concept_id}/actions/newversion"

    # Create the draft
    new_draft = requests.post(api, params=params, json={})

    # trap errors in creating the new version (not 201: created)
    if new_draft.status_code != 201:
        return {}, _zenodo_error_message(new_draft)

    if concept_id is None:
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
    metadata: dict, zenodo: dict, resources: Optional[Resources] = None
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

    zen_md["metadata"]["description"] = dataset_description(
        metadata, zenodo, render=True, resources=resources
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
    resources: Optional[Resources] = None,
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


def upload_file(
    metadata: dict,
    filepath: str,
    zenodo_filename: Optional[str] = None,
    progress_bar: bool = True,
    resources: Optional[Resources] = None,
) -> ZenodoFunctionResponseType:
    """Upload a file to Zenodo.

    Uploads the contents of a specified file to an unpublished Zenodo deposit,
    optionally using an alternative filename. If the file already exists in the deposit,
    it will be replaced.

    Args:
        metadata: The Zenodo metadata dictionary for a deposit
        filepath: The path to the file to be uploaded
        zenodo_filename: An optional alternative file name to be used on Zenodo
        progress_bar: Should the upload progress be displayed
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)
    params = zres["ztoken"]

    # Check the file and get the filename if an alternative is not provided
    filepath = os.path.abspath(filepath)
    if not (os.path.exists(filepath) and os.path.isfile(filepath)):
        raise IOError(f"The file path is either a directory or not found: {filepath} ")

    if zenodo_filename is None:
        file_name = os.path.basename(filepath)
    else:
        file_name = zenodo_filename

    # upload the file
    # - https://gist.github.com/tyhoff/b757e6af83c1fd2b7b83057adf02c139
    file_size = os.stat(filepath).st_size
    api = f"{metadata['links']['bucket']}/{file_name}"

    with open(filepath, "rb") as file_io:
        if progress_bar:
            with tqdm(
                total=file_size, unit="B", unit_scale=True, unit_divisor=1024
            ) as upload_monitor:
                # Upload the wrapped file
                wrapped_file = CallbackIOWrapper(upload_monitor.update, file_io, "read")
                fls = requests.put(api, data=wrapped_file, params=params)
        else:
            fls = requests.put(api, data=file_io, params=params)

    # trap errors in uploading file
    # - no success or mismatch in md5 checksums
    if fls.status_code != 201:
        return {}, _zenodo_error_message(fls)

    # TODO - could this be inside with above? - both are looping over the file contents
    # https://medium.com/codex/chunked-uploads-with-binary-files-in-python-f0c48e373a91
    local_hash = _compute_md5(filepath)

    if fls.json()["checksum"] != f"md5:{local_hash}":
        return {}, "Mismatch in local and uploaded MD5 hashes"
    else:
        return fls.json(), None


def discard_deposit(
    metadata: dict, resources: Optional[Resources] = None
) -> ZenodoFunctionResponseType:
    """Discard a deposit.

    Deposits can be discarded - the associated files and metadata will be deleted and
    the Zenodo ID no longer exists. Once deposits are published to records, they cannot
    be deleted via the API - contact the Zenodo team for help.

    Args:
        metadata: The Zenodo metadata dictionary for a deposit
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)
    params = zres["ztoken"]

    delete = requests.delete(metadata["links"]["self"], params=params)

    if delete.status_code == 204:
        return {"result": "success"}, None
    else:
        return {}, _zenodo_error_message(delete)


def publish_deposit(
    zenodo: dict, resources: Optional[Resources] = None
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


def delete_file(
    metadata: dict, filename: str, resources: Optional[Resources] = None
) -> ZenodoFunctionResponseType:
    """Delete an uploaded file from an unpublished Zenodo deposit.

    Args:
        metadata: The Zenodo metadata dictionary for a deposit
        filename: The file to delete from the deposit
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

    # get a dictionary of file links
    files_dict = {f["filename"]: f["links"]["self"] for f in files.json()}

    if filename not in files_dict:
        return {}, f"{filename} is not a file in the deposit"

    # get the delete link to the file and call
    delete_api = files_dict[filename]
    file_del = requests.delete(delete_api, params=params)

    if file_del.status_code != 204:
        return {}, _zenodo_error_message(file_del)
    else:
        return {"result": "success"}, None


"""
Dataset description generation (HTML and GEMINI XML)
"""


def dataset_description(
    dataset_metadata: dict,
    zenodo_metadata: dict,
    render: bool = True,
    extra: Optional[str] = None,
    resources: Optional[Resources] = None,
) -> Union[tags.div, str]:
    """Create an HTML dataset description.

    This function turns a dataset metadata JSON into html for inclusion in
    published datasets. This content is used to populate the dataset description
    section in the Zenodo metadata. Zenodo has a limited set of permitted HTML
    tags, so this is quite simple HTML.

    The available tags are: a, p, br, blockquote, strong, b, u, i, em, ul, ol,
    li, sub, sup, div, strike. Note that `<a>` is currently only available on
    Zenodo when descriptions are uploaded programmatically as a bug in their
    web interface strips links.

    The description can be modified for specific uses by including HTML via the
    extra argument. This content is inserted below the dataset description.

    Args:
        dataset_metadata: The dataset metadata
        zenodo_metadata: The Zenodo deposit metadata
        render: Should the html be returned as text or as the underlying
            dominate.tags.div object.
        extra: Additional HTML content to include in the description.
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        Either a string of rendered HTML or a dominate.tags.div object.
    """

    # zres = _resources_to_zenodo_api(resources)
    # metadata_api = zres["mdapi"]

    # PROJECT Title and authors are added by Zenodo from zenodo metadata
    # TODO - option to include here?

    desc = tags.div()

    # Dataset summary
    desc += tags.b("Description: ")
    desc += tags.p(dataset_metadata["description"].replace("\n", "</br>"))

    # Extra
    if extra is not None:
        desc += raw(extra)

    # proj_url = URL('projects', 'project_view', args=[metadata['project_id']],
    #               scheme=True, host=True)
    # desc += P(B('Project: '), 'This dataset was collected as part of the following '
    #                          'SAFE research project: ', A(B(title), _href=proj_url))
    ##

    # Funding information
    if dataset_metadata["funders"]:
        funder_info = []

        for fnd in dataset_metadata["funders"]:
            funder_details = [fnd["body"], "(", fnd["type"]]

            if fnd["ref"]:
                funder_details.append(str(fnd["ref"]))
            if fnd["url"]:
                funder_details.append(tags.a(fnd["url"], _href=fnd["url"]))

            funder_details.append(")")
            funder_info.append(tags.li(funder_details))

        desc += [
            tags.p(
                tags.b("Funding: "),
                "These data were collected as part of research funded by: ",
                tags.ul(funder_info),
            ),
            tags.p(
                "This dataset is released under the CC-BY 4.0 licence, requiring that "
                "you cite the dataset in any outputs, but has the additional condition "
                "that you acknowledge the contribution of these funders in any outputs."
            ),
        ]

    # Permits
    if dataset_metadata["permits"]:
        desc += tags.p(
            tags.b("Permits: "),
            "These data were collected under permit from the following authorities:",
            tags.ul(
                [
                    tags.li(
                        f"{pmt['authority']} ({pmt['type']} licence {pmt['number']})"
                    )
                    for pmt in dataset_metadata["permits"]
                ]
            ),
        )

    # # XML link
    # xml_url = f"{metadata_api}/xml/{zenodo['record_id']}"

    # desc += tags.p(
    #     tags.b("XML metadata: "),
    #     "GEMINI compliant metadata for this dataset is available ",
    #     tags.a("here", href=xml_url),
    # )

    # Present a description of the file or files including 'external' files
    # (data files loaded directly to Zenodo).
    ds_files = [dataset_metadata["filename"]]
    n_ds_files = 1
    ex_files = []

    if dataset_metadata["external_files"]:
        ex_files = dataset_metadata["external_files"]
        ds_files += [f["file"] for f in ex_files]
        n_ds_files += len(ex_files)

    desc += tags.p(
        tags.b("Files: "),
        f"This dataset consists of {n_ds_files} files: ",
        ", ".join(ds_files),
    )

    # Group the sheets by their 'external' file - which is None for sheets
    # in the submitted workbook - and collect them into a dictionary by source
    # file. get() is used here for older data where external was not present.

    tables_by_source = dataset_metadata["dataworksheets"]

    # Now group into a dictionary keyed by external source file - cannot sort
    # None (no comparison operators) so use a substitute
    tables_by_source.sort(key=lambda sh: sh.get("external") or False)
    tables_by_source = groupby(
        tables_by_source, key=lambda sh: sh.get("external") or False
    )
    tables_by_source = {g: list(v) for g, v in tables_by_source}

    # We've now got a set of files (worksheet + externals) and a dictionary of table
    # descriptions that might have an entry for each file.

    # Report the worksheet first
    desc += tags.p(tags.b(dataset_metadata["filename"]))

    # Report internal tables
    if False in tables_by_source:
        int_tabs = tables_by_source[False]
        desc += tags.p(
            f"This file contains dataset metadata and {len(int_tabs)} data tables:"
        )
        desc += tags.ol([tags.li(table_description(tab)) for tab in int_tabs])
    else:
        # No internal tables at all.
        desc += tags.p("This file only contains metadata for the files below")

    # Report on the other files
    for exf in ex_files:
        desc += tags.p(
            tags.b(exf["file"]), tags.p(f"Description: {exf['description']}")
        )

        if exf["file"] in tables_by_source:
            # Report table description
            ext_tabs = tables_by_source[exf["file"]]
            desc += tags.p(f"This file contains {len(ext_tabs)} data tables:")
            desc += tags.ol([tags.li(table_description(tab)) for tab in ext_tabs])

    # Add extents if populated
    if dataset_metadata["temporal_extent"] is not None:
        desc += tags.p(
            tags.b("Date range: "),
            "{0[0]} to {0[1]}".format(
                [x[:10] for x in dataset_metadata["temporal_extent"]]
            ),
        )
    if dataset_metadata["latitudinal_extent"] is not None:
        desc += tags.p(
            tags.b("Latitudinal extent: "),
            "{0[0]:.4f} to {0[1]:.4f}".format(dataset_metadata["latitudinal_extent"]),
        )
    if dataset_metadata["longitudinal_extent"] is not None:
        desc += tags.p(
            tags.b("Longitudinal extent: "),
            "{0[0]:.4f} to {0[1]:.4f}".format(dataset_metadata["longitudinal_extent"]),
        )

    # Find taxa data from each database (if they exist)
    gbif_taxon_index = dataset_metadata.get("gbif_taxa")
    ncbi_taxon_index = dataset_metadata.get("ncbi_taxa")

    # When NCBI is absent use the old format for backwards compatibility
    if gbif_taxon_index or ncbi_taxon_index:
        desc += tags.p(
            tags.b("Taxonomic coverage: "),
            tags.br(),
            "This dataset contains data associated with taxa and these have been "
            "validated against appropriate taxonomic authority databases.",
        )

    if gbif_taxon_index:
        desc += tags.p(
            tags.u("GBIF taxa details: "),
            tags.br(),
            tags.br(),
            "The following taxa were validated against the GBIF backbone dataset."
            "If a dataset uses a synonym, the accepted usage is shown followed by the "
            "dataset usage in brackets. Taxa that cannot be validated, including new "
            "species and other unknown taxa, morphospecies, functional groups and "
            "taxonomic levels not used in the GBIF backbone are shown in square "
            "brackets.",
            taxon_index_to_text(gbif_taxon_index, True, auth="GBIF"),
        )

    if ncbi_taxon_index:
        desc += tags.p(
            tags.u("NCBI taxa details: "),
            tags.br(),
            tags.br(),
            "The following taxa were validated against the NCBI taxonomy dataset."
            " If a dataset uses a synonym, the accepted usage is shown followed by the "
            "dataset usage in brackets. Taxa that cannot be validated, e.g. new or "
            "unknown species are shown in square brackets. Non-backbone taxonomic "
            "ranks (e.g. strains or subphyla) can be validated using the NCBI "
            "database. However, they will only be shown if the user explicitly "
            "provided a non-backbone taxon. When they are shown they will be "
            "accompanied by an message stating their rank.",
            taxon_index_to_text(ncbi_taxon_index, True, auth="NCBI"),
        )

    if render:
        return desc.render()
    else:
        return desc


def table_description(tab: dict) -> tags.div:
    """Convert a dict containing table contents into an HTML table.

    Function to return a description for an individual source file in a dataset.
    Typically datasets only have a single source file - the Excel workbook that
    also contains the metadata - but they may also report on external files loaded
    directly to Zenodo, and which uses the same mechanism.

    Args:
        tab: A dict describing a data table

    Returns:
        A `dominate.tags.div` instance containing an HTML description of the table
    """

    # table summary
    tab_desc = tags.div(
        tags.p(tags.b(tab["title"]), f" (described in worksheet {tab['name']})"),
        tags.p(f"Description: {tab['description']}"),
        tags.p(f"Number of fields: {tab['max_col'] - 1}"),
    )

    # The explicit n_data_row key isn't available for older records
    if "n_data_row" in tab:
        if tab["n_data_row"] == 0:
            tab_desc += tags.p(
                "Number of data rows: Unavailable (table metadata description only)."
            )
        else:
            tab_desc += tags.p(f"Number of data rows: {tab['n_data_row']}")
    else:
        tab_desc += tags.p(
            f"Number of data rows: {tab['max_row'] - len(tab['descriptors'])}"
        )

    # add fields
    tab_desc += tags.p("Fields: ")

    # fields summary
    flds = tags.ul()
    for each_fld in tab["fields"]:
        flds += tags.li(
            tags.b(each_fld["field_name"]),
            f": {each_fld['description']} (Field type: {each_fld['field_type']})",
        )

    tab_desc += flds

    return tab_desc


def generate_inspire_xml(
    dataset_metadata: dict,
    zenodo_metadata: dict,
    resources: Resources,
    lineage_statement: Optional[str] = None,
) -> bytes:
    """Convert dataset and zenodo metadata into GEMINI XML.

    Produces an INSPIRE/GEMINI formatted XML record from dataset metadata,
    and Zenodo record metadata using a template XML file. The dataset URL
    defaults to the Zenodo record but can be replaced if a separate URL (such as
    a project specific website) is used. The Gemini XML standard requires a
    statement about the lineage of a dataset - if this is not provided to this
    function it will appear as "Not provided".

    Args:
        dataset_metadata: A dictionary of the dataset metadata
        zenodo_metadata: A dictionary of the Zenodo record metadata
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.
        lineage_statement: An optional lineage statement about the data.

    Returns:
        A string containing GEMINI compliant XML.
    """

    # Get resource configuration
    # zres = _resources_to_zenodo_api(resources)

    # parse the XML template and get the namespace map
    template = Path(__file__).parent / "gemini_xml_template.xml"
    tree = etree.parse(template)
    root = tree.getroot()
    nsmap = root.nsmap

    pub_date = dt.fromisoformat(zenodo_metadata["metadata"]["publication_date"])

    # Use find and XPATH to populate the template, working through from the top of the
    # file

    # file identifier
    root.find("./gmd:fileIdentifier/gco:CharacterString", nsmap).text = "zenodo." + str(
        zenodo_metadata["id"]
    )

    # date stamp (not clear what this is - taken as publication date)
    root.find("./gmd:dateStamp/gco:DateTime", nsmap).text = pub_date.isoformat()

    # Now zoom to the data identification section
    data_id = root.find(".//gmd:MD_DataIdentification", nsmap)

    # CITATION
    citation = data_id.find("gmd:citation/gmd:CI_Citation", nsmap)
    citation.find("gmd:title/gco:CharacterString", nsmap).text = dataset_metadata[
        "title"
    ]
    citation.find(
        "gmd:date/gmd:CI_Date/gmd:date/gco:Date", nsmap
    ).text = pub_date.date().isoformat()

    # URIs - a dataset URL and the DOI
    # dataset_url = (
    #     f"{zres['mdapi']}/datasets/view_dataset"
    #     f"?id={dataset_metadata['zenodo_record_id']}"
    # )
    # citation.find(
    #     "gmd:identifier/gmd:MD_Identifier/gmd:code/gco:CharacterString", nsmap
    # ).text = dataset_url
    citation.find(
        "gmd:identifier/gmd:RS_Identifier/gmd:code/gco:CharacterString", nsmap
    ).text = zenodo_metadata["doi_url"]

    # The citation string
    authors = [au["name"] for au in dataset_metadata["authors"]]
    author_string = ", ".join(authors)
    if len(authors) > 1:
        author_string = author_string.replace(", " + authors[-1], " & " + authors[-1])

    cite_string = "{} ({}) {} [Dataset] {}".format(
        author_string,
        pub_date.year,
        dataset_metadata["title"],
        zenodo_metadata["doi_url"],
    )

    citation.find(
        "gmd:otherCitationDetails/gco:CharacterString", nsmap
    ).text = cite_string

    # ABSTRACT
    data_id.find("gmd:abstract/gco:CharacterString", nsmap).text = dataset_metadata[
        "description"
    ]

    # KEYWORDS
    # - find the container node for the free keywords
    keywords = data_id.find("./gmd:descriptiveKeywords/gmd:MD_Keywords", nsmap)
    # - get the placeholder node
    keywd_node = keywords.getchildren()[0]
    # - duplicate it if needed
    for new_keywd in range(len(dataset_metadata["keywords"]) - 1):
        keywords.append(copy.deepcopy(keywd_node))
    # populate the nodes
    for key_node, val in zip(keywords.getchildren(), dataset_metadata["keywords"]):
        key_node.find("./gco:CharacterString", nsmap).text = val

    # AUTHORS - find the point of contact with author role from the template and its
    # index using xpath() here to access full xpath predicate search.
    au_xpath = (
        "./gmd:pointOfContact[gmd:CI_ResponsibleParty/"
        "gmd:role/gmd:CI_RoleCode='author']"
    )
    au_node = data_id.xpath(au_xpath, namespaces=nsmap)[0]
    au_idx = data_id.index(au_node)

    # - duplicate it if needed into the tree
    for n in range(len(dataset_metadata["authors"]) - 1):
        data_id.insert(au_idx, copy.deepcopy(au_node))

    # now populate the author nodes, there should now be one for each author
    au_ls_xpath = (
        "./gmd:pointOfContact[gmd:CI_ResponsibleParty/"
        "gmd:role/gmd:CI_RoleCode='author']"
    )
    au_node_list = data_id.xpath(au_ls_xpath, namespaces=nsmap)

    for au_data, au_node in zip(dataset_metadata["authors"], au_node_list):
        resp_party = au_node.find("gmd:CI_ResponsibleParty", nsmap)
        resp_party.find("gmd:individualName/gco:CharacterString", nsmap).text = au_data[
            "name"
        ]
        resp_party.find(
            "gmd:organisationName/gco:CharacterString", nsmap
        ).text = au_data["affiliation"]
        contact_info = resp_party.find("gmd:contactInfo/gmd:CI_Contact", nsmap)
        email_path = (
            "gmd:address/gmd:CI_Address/gmd:electronicMailAddress/gco:CharacterString"
        )
        contact_info.find(email_path, nsmap).text = au_data["email"]

        # handle orcid resource
        orcid = contact_info.find("gmd:onlineResource", nsmap)
        if au_data["orcid"] is None:
            contact_info.remove(orcid)
        else:
            orcid.find("gmd:CI_OnlineResource/gmd:linkage/gmd:URL", nsmap).text = (
                "http://orcid.org/" + au_data["orcid"]
            )

    # CONSTRAINTS
    # update the citation information in the second md constraint
    md_path = (
        "gmd:resourceConstraints/gmd:MD_Constraints/"
        "gmd:useLimitation/gco:CharacterString"
    )
    md_constraint = data_id.find(md_path, nsmap)
    md_constraint.text += cite_string

    # embargo or not?
    embargo_path = (
        "gmd:resourceConstraints/gmd:MD_LegalConstraints/"
        "gmd:otherConstraints/gco:CharacterString"
    )
    if dataset_metadata["access"] == "embargo":
        data_id.find(embargo_path, nsmap).text = (
            f"This data is under embargo until {dataset_metadata['embargo_date']}."
            "After that date there are no restrictions to public access."
        )
    elif dataset_metadata["access"] == "restricted":
        data_id.find(embargo_path, nsmap).text = (
            "This dataset is currently not publicly available, please contact the "
            "Zenodo community owner to request access."
        )
    else:
        data_id.find(
            embargo_path, nsmap
        ).text = "There are no restrictions to public access."

    # EXTENTS
    temp_extent = root.find(".//gmd:EX_TemporalExtent", nsmap)
    temp_extent.find(".//gml:beginPosition", nsmap).text = dataset_metadata[
        "temporal_extent"
    ][0][:10]
    temp_extent.find(".//gml:endPosition", nsmap).text = dataset_metadata[
        "temporal_extent"
    ][1][:10]

    geo_extent = root.find(".//gmd:EX_GeographicBoundingBox", nsmap)
    geo_extent.find("./gmd:westBoundLongitude/gco:Decimal", nsmap).text = str(
        dataset_metadata["longitudinal_extent"][0]
    )
    geo_extent.find("./gmd:eastBoundLongitude/gco:Decimal", nsmap).text = str(
        dataset_metadata["longitudinal_extent"][1]
    )
    geo_extent.find("./gmd:southBoundLatitude/gco:Decimal", nsmap).text = str(
        dataset_metadata["latitudinal_extent"][0]
    )
    geo_extent.find("./gmd:northBoundLatitude/gco:Decimal", nsmap).text = str(
        dataset_metadata["latitudinal_extent"][1]
    )

    # Dataset transfer options: direct download and dataset view on SAFE website
    distrib = root.find("gmd:distributionInfo/gmd:MD_Distribution", nsmap)
    distrib.find(
        (
            "gmd:transferOptions[1]/gmd:MD_DigitalTransferOptions/gmd:onLine/"
            "gmd:CI_OnlineResource/gmd:linkage/gmd:URL"
        ),
        nsmap,
    ).text = zenodo_metadata["files"][0]["links"]["download"]
    distrib.find(
        (
            "gmd:transferOptions[2]/gmd:MD_DigitalTransferOptions/gmd:onLine/"
            "gmd:CI_OnlineResource/gmd:linkage/gmd:URL"
        ),
        nsmap,
    ).text += str(dataset_metadata["zenodo_record_id"])

    # LINEAGE STATEMENT
    # lineage = (
    #     "This dataset was collected as part of a research project based at The"
    #     " SAFE Project. For details of the project and data collection, see the "
    #     "methods information contained within the datafile and the project "
    #     "website: "
    # ) + URL("projects", "view_project", args=record.project_id, scheme=True,
    # host=True)

    root.find(
        (
            "gmd:dataQualityInfo/gmd:DQ_DataQuality/gmd:lineage/gmd:LI_Lineage/"
            "gmd:statement/gco:CharacterString"
        ),
        nsmap,
    ).text = lineage_statement

    # return the string contents
    return etree.tostring(tree)


def download_ris_data(
    resources: Optional[Resources] = None, ris_file: Optional[str] = None
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

    if ris_file and os.path.exists(ris_file):
        with open(ris_file, "r") as bibliography_file:
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
            raise IOError("Cannot access Zenodo API")
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
        if os.path.exists(ris_file):
            LOGGER.info(f"Appending RIS data for {len(data)} new records to {ris_file}")
            write_mode = "a"
        else:
            LOGGER.info(f"Writing RIS data for {len(data)} records to {ris_file}")
            write_mode = "w"

        with open(ris_file, write_mode) as ris_file_out:
            for this_entry in data:
                ris_file_out.write(this_entry)


def sync_local_dir(
    datadir: str,
    xlsx_only: bool = True,
    replace_modified: bool = False,
    resources: Optional[Resources] = None,
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
    def _get_file(url: str, outf: str, params: Optional[dict] = None) -> None:
        """Download a file from a URL."""
        resource = requests.get(url, params=params, stream=True)

        with open(outf, "wb") as outf_obj:
            shutil.copyfileobj(resource.raw, outf_obj)

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)
    zenodo_api = zres["zapi"]
    params = zres["ztoken"]

    # The dir argument should be an existing path
    if not (os.path.exists(datadir) and os.path.isdir(datadir)):
        raise IOError(f"{datadir} is not an existing directory")

    # Get the configured metadata api
    api = zres["mdapi"]

    # Check for an existing API url file and check it is congruent with config
    url_file = os.path.join(datadir, "url.json")

    if os.path.exists(url_file):
        with open(url_file, "r") as urlf:
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
    _get_file(f"{api}/api/index", os.path.join(datadir, "index.json"))
    _get_file(f"{api}/api/gazetteer", os.path.join(datadir, "gazetteer.geojson"))
    _get_file(
        f"{api}/api/location_aliases", os.path.join(datadir, "location_aliases.csv")
    )

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
        rec_dir = os.path.join(datadir, con_rec_id, rec_id)
        if not os.path.exists(rec_dir):
            LOGGER.info("Creating directory")
            os.makedirs(rec_dir)
        else:
            LOGGER.info("Directory found")

        # loop over the files in the record
        for this_file in dep["files"]:
            if xlsx_only and not this_file["filename"].endswith(".xlsx"):
                LOGGER.info(f"Skipping non-excel file {this_file['filename']}")
                continue

            LOGGER.info(f"Processing {this_file['filename']}")
            FORMATTER.push()

            outf = os.path.join(rec_dir, this_file["filename"])
            local_copy = os.path.exists(outf)

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
        metadata = os.path.join(rec_dir, f"{rec_id}.json")
        if os.path.exists(metadata):
            LOGGER.info("JSON Metadata found")
        else:
            LOGGER.info("Downloading JSON metadata ")
            _get_file(f"{api}/api/record/{rec_id}", metadata)

        FORMATTER.pop()

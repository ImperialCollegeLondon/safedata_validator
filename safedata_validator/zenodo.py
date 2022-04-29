import os
import datetime
from lxml import etree
from regex import P
import simplejson
import copy
from io import StringIO
import requests
# from safe_web_global_functions import safe_mailer
from itertools import groupby
from shapely import geometry
import hashlib
# from gluon.serializers import json
import rispy
import shutil

from tqdm import tqdm
from tqdm.utils import CallbackIOWrapper
from contextlib import nullcontext

from typing import Tuple, Union


from safedata_validator.resources import Resources
from safedata_validator.logger import LOGGER, FORMATTER
"""
This module provides functions to:
1. handle the publication of datasets after they have been validated 
   using safedata_validate,
2. maintain local copies of datasets in the folder structure expected 
   by the safedata R package.
3. compile a RIS format bibliographic file for published datasets.
"""

def _resources_to_zenodo_api(resources: Resources = None):
    """Private function to return the Zenodo API and access token
    from the configuration.
    """

    # Get resource configuration
    if resources is None:
        resources = Resources()

    # Get the Zenodo API 
    if resources.zenodo.use_sandbox:
        zenodo_api = resources.zenodo.zenodo_sandbox_api
        token = resources.zenodo.zenodo_sandbox_token
    else:
        zenodo_api = resources.zenodo.zenodo_api
        token = resources.zenodo.zenodo_token

    return zenodo_api, {'access_token': token}, resources.zenodo.community_name

def _compute_md5(fname:str) -> str:
    """Calculate a local file md5 hash
    """
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as fname_obj:
        for chunk in iter(lambda: fname_obj.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def _zenodo_error_message(response) -> str:
    
    return f"{response.json()['message']} ({response.json()['status']})"

"""
Zenodo action functions
"""

def get_deposit(deposit_id: int,
                resources: Resources = None
                ) -> Tuple[Union[dict, None], Union[str, None]]:
    """
    Download the metadata of a Zenodo deposit.

    Args:
        deposit_id: The Zenodo record id of an existing dataset.
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.
    
    Returns:
        A 2 tuple:
            * A dictionary containing the response content or None on error
            * An error message on failure or None on success
    """

    zenodo_api, params, _ = _resources_to_zenodo_api(resources)

    # request the deposit
    dep = requests.get(f'{zenodo_api}/deposit/depositions/{deposit_id}',
                       params=params, json={})

    # check for success and return the information.
    if dep.status_code == 200:
        return dep.json(), None
    else:
        return None, _zenodo_error_message(dep)


def create_deposit(deposit_id: int = None,
                   resources: Resources = None
                    ) -> Tuple[Union[dict, None], Union[str, None]]:
    """
    Function to create a new deposit draft, possibly as a new version of an
    existing published record

    Args:
        deposit_id: An optional id of a published record to create a new version
            of an existing dataset
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.
    Returns:
        A 2 tuple:
            * A dictionary containing the response content or None on error
            * An error message on failure or None on success
    """

    # Get resource configuration
    zenodo_api, params, _ = _resources_to_zenodo_api(resources)

    # get the correct draft api
    if deposit_id is None:
        api = f'{zenodo_api}/deposit/depositions'
    else:
        api = f'{zenodo_api}/deposit/depositions/{deposit_id}/actions/newversion'
    
    # Create the draft
    new_draft = requests.post(api, params=params, json={})

    # trap errors in creating the new version (not 201: created)
    if new_draft.status_code != 201:
        return None, _zenodo_error_message(new_draft)

    if deposit_id is None:
        return new_draft.json(), None

    # For new versions, the response is an update to the existing copy,
    # so need to separately retrieve the new draft
    api = new_draft.json()['links']['latest_draft']
    dep = requests.get(api, params=params, json={})

    # trap errors in creating the resource - successful creation of new version
    #  drafts returns 200
    if dep.status_code != 200:
        return None, _zenodo_error_message(dep)
    else:
        return dep.json(), None

# FIXME
def upload_metadata(links, token, record, zenodo_id):
    """
    Function to turn a dataset row record into a Zenodo metadata JSON and upload
    it to a deposit.

    Args:
        links: The links dictionary from a created deposit
        token: The access token to be used
        record: The database record containing the metadata to be uploaded.
        zenodo_id: The ID of the zenodo draft being created, to be used as a
            key for the GEMINI XML link.

    Returns:
        An integer indicating success (0) or failure (1) and either the
        deposit links dictionary or an error message
    """

    # extract the metadata from the record
    metadata = record.dataset_metadata['metadata']

    # basic contents
    zen_md = {
        'metadata': {
            "upload_type": "dataset",
            "publication_date": datetime.date.today().isoformat(),
            "title": metadata['title'],
            "keywords": metadata['keywords'],
            "license": 'cc-by',
            "contributors": [
                {"name": "The SAFE Project", "type": "ContactPerson",
                 "affiliation": "Imperial College London",
                 "orcid": "0000-0003-3378-2814"},
            ],
            "communities": [{"identifier": "safe"}]
        }
    }

    # set up the access rights
    if metadata['access'].lower() == 'embargo':
        zen_md['metadata']['access_right'] = 'embargoed'
        zen_md['metadata']['embargo_date'] = metadata['embargo_date']
    elif metadata['access'].lower() == 'open':
        zen_md['metadata']['access_right'] = 'open'
    elif metadata['access'].lower() == 'restricted':
        zen_md['metadata']['access_right'] = 'restricted'
        zen_md['metadata']['access_conditions'] = metadata['access_conditions']
    else:
        raise ValueError('Unknown access status')

    # set up the dataset creators - the format has already been checked and names
    # should be present and correct. Everything else is optional, so strip None
    # values and pass the rest to Zenodo
    zen_md['metadata']['creators'] = [
        {ky: auth[ky] for ky in auth if auth[ky] is not None and ky != 'email'}
        for auth in metadata['authors']]

    zen_md['metadata']['description'] = str(dataset_description(record, gemini_id=zenodo_id))

    # attach the metadata to the deposit resource
    mtd = requests.put(links['self'], params=token, data=simplejson.dumps(zen_md),
                       headers={"Content-Type": "application/json"})

    # trap errors in uploading metadata and tidy up
    if mtd.status_code != 200:
        return 1, mtd.json()
    else:
        return 0, 'success'


# FIXME
def update_published_metadata(zenodo_record_id, new_values):
    """
    Function to update the metadata on a published deposit. Currently used to
    modify the access status of deposit.

    Args:
        zenodo_record_id: The record id of the deposit to be updated
        new_values: A dictionary of new values to be substituted in to
             the existing Zenodo deposition metadata resource.
    """

    # load the correct API and token
    api, token = get_zenodo_api()

    code, dep = get_deposit(api, token, zenodo_record_id)

    if code != 0:
        return 1, dep

    links = dep['links']
    metadata = dep['metadata']

    # Unlock the published deposit for editing
    edt = requests.post(links['edit'], params=token)

    if edt.status_code != 201:
        return 1, edt.json()

    # Amend the metadata
    for key, val in new_values.items():
        if val is not None:
            metadata[key] = val
        elif key in metadata:
            metadata.pop(key)

    # If any API calls from now fail, we need to tidy up the edit
    # status of the record, or it will block subsequent attempts

    upd = requests.put(links['self'], params=token,
                       headers={"Content-Type": "application/json"},
                       data=simplejson.dumps({'metadata': metadata}))

    success_so_far = 0 if upd.status_code != 200 else 1
    ret = upd.json()

    # Republish to save the changes
    if success_so_far:
        pub = requests.post(links['publish'], params=token)
        success_so_far = 0 if pub.status_code != 202 else 1
        ret = pub.json()

    # If all steps have been successful, return a 0 code, otherwise
    # try to discard the edits and return the most recent failure
    # notice

    if success_so_far:
        return 0, ret
    else:
        dsc = requests.post(links['discard'], params=token)
        success_so_far = 0 if dsc.status_code != 201 else 1
        if not success_so_far:
            ret = dsc.json()

        return 1, ret


def upload_file(metadata: dict, 
                filepath: str, 
                zenodo_filename: str = None,
                progress_bar: bool = True,
                resources: Resources = None) -> Tuple[Union[dict, None], Union[str, None]]:
    """
    Upload the contents of a specified file to an unpublished Zenodo deposit,
    optionally using an alternative filename. If the file already exists in the
    deposit, it will be replaced.

    Args:
        metadata: The Zenodo metadata dictionary for a deposit
        file: The path to the file to be uploaded
        zenodo_filename: An optional alternative file name to be used on Zenodo
        progress_bar: Should the upload progress be displayed
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        A 2 tuple:
            * A dictionary containing the response content or None on error
            * An error message on failure or None on success
    """

    # Get resource configuration
    _, params, _ = _resources_to_zenodo_api(resources)

    # Check the file and get the filename if an alternative is not provided
    if not (os.path.exists(filepath) and os.path.isfile(filepath)):
        raise IOError(f'The file path is either a directory or not found: {filepath} ')

    if zenodo_filename is None:
        file_name = os.path.basename(filepath)
    else:
        file_name = zenodo_filename
    
    # upload the file
    file_size = os.stat(filepath).st_size
    api = f"{metadata['links']['bucket']}/{file_name}"

    with open(filepath, 'rb') as fp:

        if progress_bar:
            with tqdm(total=file_size, unit="B", unit_scale=True, unit_divisor=1024) as cm:
                # Upload the wrapped file
                wrapped_file = CallbackIOWrapper(cm.update, fp, "read")
                fls = requests.put(api, data=fp, params=params)
        else:
            fls = requests.put(api, data=fp, params=params)

    # trap errors in uploading file
    # - no success or mismatch in md5 checksums
    if fls.status_code != 200:
        return None, _zenodo_error_message(fls)

    # TODO - could this be inside with above? - both are looping over the file contents
    local_hash = _compute_md5(filepath)
    
    if fls.json()['checksum'] != f'md5:{local_hash}':
        return None, "Mismatch in local and uploaded MD5 hashes"
    else:
        return fls, None



def discard_deposit(metadata: dict,
                    resources: Resources = None
                    ) -> Tuple[Union[dict, None], Union[str, None]]:
    """
    Discard an _unpublished_ deposit.

    Args:
        metadata: The Zenodo metadata dictionary for a deposit
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        A 2 tuple:
            * A dictionary containing the response content or None on error
            * An error message on failure or None on success
    """

    # Get resource configuration
    _, params, _ = _resources_to_zenodo_api(resources)

    delete = requests.delete(metadata['links']['self'], params=params)

    if delete.status_code == 204:
        return {'result': 'success'}, None
    else:
        return None, _zenodo_error_message(delete)

# FIXME
def publish_deposit(links, token):
    """
    Function to publish a created deposit .

    Args:
        links: The links dictionary from a created deposit
        token: The access token to be used

    Returns:
        An integer indicating success (0) or failure (1) and either the
        deposit links dictionary or an error message
    """

    # publish
    pub = requests.post(links['publish'], params=token)

    # trap errors in publishing, otherwise return the publication metadata
    if pub.status_code != 202:
        return 1, pub.json()
    else:
        return 0, pub.json()


def delete_file(metadata: dict,
                filename: str,
                resources: Resources = None
                ) -> Tuple[Union[dict, None], Union[str, None]]:
    """
    Delete an uploaded file from an unpublished Zenodo deposit.

    Args:
        metadata: The Zenodo metadata dictionary for a deposit
        filename: The file to delete from the deposit
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        A 2 tuple:
            * A dictionary containing the response content or None on error
            * An error message on failure or None on success
    """

    # Get resource configuration
    _, params, _ = _resources_to_zenodo_api(resources)

    # get an up to date list of existing files (metadata
    # might be outdated)
    files = requests.get(metadata['links']['files'], params=params)

    # check the result of the files request
    if files.status_code != 200:
        # failed to get the files
        return None, _zenodo_error_message(files)

    # get a dictionary of file links
    files = {f['filename']: f['links']['self'] for f in files.json()}

    if filename not in files:
        return None, f"{filename} is not a file in the deposit"
    
    # get the delete link to the file and call
    delete_api = files[filename]
    file_del = requests.delete(delete_api, params=params)

    if file_del.status_code != 204:
        return None, _zenodo_error_message(file_del)
    else:
        return {'result': 'success'}, None



"""
Description functions
"""


def taxon_index_to_text(taxa):
    """
    Turns the taxon index for a dataset into a text representation
    of the taxonomic hierarchy used in the dataset. Takes a list
    of dicts keyed by the fields of dataset_taxa - this could come
    from db.datasets_taxa for a published dataset but also from
    the dataset_metadata of submitted datasets.
    """

    def indent(n):

        return ('&ensp;-&ensp;' * n)

    def format_name(tx):

        # format the canonical name
        if tx['taxon_rank'] in ['genus', 'species', 'subspecies']:
            return '<i>{}</i>'.format(tx['taxon_name'])
        elif tx['taxon_rank'] in ['morphospecies', 'functional group']:
            return '[{}]'.format(tx['taxon_name'])
        else:
            return tx['taxon_name']

    # Container to hold the output
    html = StringIO()

    # group by parent taxon, subsitituting 0 for None
    taxa.sort(key=lambda x: x['gbif_parent_id'] or 0)
    grouped = {k: list(v) for k, v in groupby(taxa, lambda x: x['gbif_parent_id'])}

    # start the stack with the kingdoms - these taxa will have None as a parent
    stack = [{'current': grouped[None][0], 'next': grouped[None][1:]}]

    while stack:

        # Handle the current top of the stack: format the canonical name
        current = stack[-1]['current']
        canon_name = format_name(current)

        # Look for a non-None entry in next that shares the same worksheet name
        next_ws_names = [tx['worksheet_name'] for tx in stack[-1]['next']
                         if tx['worksheet_name'] is not None]

        if current['worksheet_name'] in next_ws_names:
            # pop out the matching entry and find which is 'accepted'
            name_pair = stack[-1]['next'].pop(next_ws_names.index(current['worksheet_name']))
            if current['gbif_status'] == 'accepted':
                as_name = format_name(name_pair)
                as_status = name_pair['gbif_status']
            else:
                as_name = canon_name
                as_status = current['gbif_status']
                canon_name = format_name(name_pair)

            txt = '{} {} (as {}: {})<br>'.format(indent(len(stack)), canon_name, as_status, as_name)
        else:
            txt = '{} {} <br>'.format(indent(len(stack)), canon_name)

        html.write(txt)

        # Is this taxon a parent for other taxa - if so add that taxon to the top of
        # the stack, otherwise start looking for a next taxon to push onto the stack.
        # If there is none at the top, pop and look down.
        parent_id = current['gbif_id']
        if parent_id in grouped:
            stack.append({'current': grouped[parent_id][0], 'next': grouped[parent_id][1:]})
        else:
            while stack:
                push = stack.pop()
                if push['next']:
                    stack.append({'current': push['next'][0], 'next': push['next'][1:]})
                    break

    return XML(html.getvalue())


def dataset_description(record, gemini_id=None):
    """
    Function to turn a dataset metadata JSON into html to send to Zenodo,
    to populate the dataset view and to provide for submitted datasets.
    Zenodo has a limited set of permitted HTML tags, so this is quite simple
    HTML, but having the exact same information and layout makes sense.

    Available tags:
    a, p, br, blockquote, strong, b, u, i, em, ul, ol, li, sub, sup, div, strike.

    Note that <a> is currently only available on Zenodo when descriptions are
    uploaded programatically. A bug in their web interface strips links.

    Args:
        record: The db record for the dataset (a row from db.submitted_datasets or db.published_datasets)
        gemini_id: Should the description include a link to the GEMINI XML
            service? This isn't available on the site datasets page until a dataset
            is published as it contains links to Zenodo, but should also be included
            in the description uploaded to Zenodo.
    """

    db = current.db

    # Get shortcut to metadata and taxon_index, handling different data structures
    # published and submitted datasets
    if 'taxa' in record.dataset_metadata:
        metadata = record.dataset_metadata['metadata']
        taxon_index = record.dataset_metadata['taxa']
        taxon_index = [dict(list(zip(['worksheet_name', 'gbif_id', 'gbif_parent_id',
                                      'taxon_name', 'taxon_rank', 'gbif_status'],
                                     tx))) for tx in taxon_index]
    else:
        metadata = record.dataset_metadata
        taxon_index = record.dataset_taxa.select().as_list()

    # - get a project link back to the safe website
    qry = db((db.project_id.id == metadata['project_id']))
    proj = qry.select(
        left=db.project_id.on(db.project_id.project_details_id == db.project_details.id))
    title = proj.first().project_details.title

    # dataset summary
    desc = CAT(B('Description: '), P(XML(metadata['description'].replace('\n', '<br>'))))

    proj_url = URL('projects', 'project_view', args=[metadata['project_id']],
                   scheme=True, host=True)
    desc += P(B('Project: '), 'This dataset was collected as part of the following '
                              'SAFE research project: ', A(B(title), _href=proj_url))

    # Funding information
    if metadata['funders']:
        funder_info = []
        for fnd in metadata['funders']:
            this_funder = fnd['type']
            if fnd['ref']:
                this_funder = CAT(this_funder, ', ' + str(fnd['ref']))
            if fnd['url']:
                this_funder = CAT(this_funder, ', ', A(fnd['url'], _href=fnd['url']))

            funder_info.append(LI(CAT(fnd['body'], ' (', this_funder, ')')))

        desc += P(B('Funding: '), 'These data were collected as part of '
                                  'research funded by: ', UL(funder_info),
                  P('This dataset is released under the CC-BY 4.0 licence, requiring that '
                    'you cite the dataset in any outputs, but has the additional condition '
                    'that you acknowledge the contribution of these funders in any outputs.'))

    # Permits
    if metadata['permits']:
        desc += P(B('Permits: '), 'These data were collected under permit from the following authorities:',
                  UL([LI('{authority} ({type} licence {number})'.format(**pmt)) for pmt in metadata['permits']]))

    # Can't get the XML metadata link unless it is published, since that
    # contains references to the zenodo record
    if gemini_id is not None:
        md_url = URL('datasets', 'xml_metadata', vars={'id': gemini_id}, scheme=True, host=True)
        desc += P(B('XML metadata: '),
                  'GEMINI compliant metadata for this dataset is available ',
                  A('here', _href=md_url))

    # Present a description of the file or files including 'external' files
    # (data files loaded directly to Zenodo).
    if metadata['external_files']:
        ex_files = metadata['external_files']
        desc += P(B('Files: '), 'This dataset consists of ', len(ex_files) + 1, ' files: ',
                  ', '.join([metadata['filename']] + [f['file'] for f in ex_files]))
    else:
        ex_files = []
        desc += P(B('Files: '), 'This consists of 1 file: ', metadata['filename'])

    # Group the sheets by their 'external' file - which is None for sheets
    # in the submitted workbook - and collect them into a dictionary by source file
    tables_by_source = metadata['dataworksheets']

    # Files submitted using early versions of the dataset submission process
    # don't have external in their worksheet dictionaries (but none of those will
    # have external files).
    for tab in tables_by_source:
        if 'external' not in tab:
            tab['external'] = None

    # now group into a dictionary keyed by source file, substituting empty string
    # for None if required.
    tables_by_source.sort(key=lambda sh: sh['external'] or '')
    tables_by_source = groupby(tables_by_source, key=lambda sh: sh['external'])
    tables_by_source = {g: list(v) for g, v in tables_by_source}

    # We've now got a set of files (worksheet + externals) and a dictionary of table
    # descriptions that might have an entry for each file.

    # Report the worksheet first
    desc += P(B(metadata['filename']))

    if None in tables_by_source:
        # Report internal tables
        desc += P('This file contains dataset metadata and '
                  '{} data tables:'.format(len(tables_by_source[None])))
        table_ol = OL()
        for tab in tables_by_source[None]:
            table_ol.append(LI(table_description(tab)))

        desc += table_ol
    else:
        # No internal tables at all.
        desc += P('This file only contains metadata for the files below')

    # Report on the other files
    for exf in ex_files:
        desc += P(B(exf['file'])) + P('Description: ' + exf['description'])

        if exf['file'] in tables_by_source:
            # Report table description
            desc += P('This file contains {} data tables:'.format(len(tables_by_source[exf['file']])))
            table_ol = OL()
            for tab in tables_by_source[exf['file']]:
                table_ol.append(LI(P(table_description(tab))))

            desc += table_ol

    # Add extents if populated
    if metadata['temporal_extent'] is not None:
        desc += P(B('Date range: '),
                  '{0[0]} to {0[1]}'.format([x[:10] for x in metadata['temporal_extent']]))
    if metadata['latitudinal_extent'] is not None:
        desc += P(B('Latitudinal extent: '),
                  '{0[0]:.4f} to {0[1]:.4f}'.format(metadata['latitudinal_extent']))
    if metadata['longitudinal_extent'] is not None:
        desc += P(B('Longitudinal extent: '),
                  '{0[0]:.4f} to {0[1]:.4f}'.format(metadata['longitudinal_extent']))
    if taxon_index:
        desc += CAT(P(B('Taxonomic coverage: '), BR(),
                      ' All taxon names are validated against the GBIF backbone taxonomy. If a '
                      'dataset uses a synonym, the accepted usage is shown followed by the dataset '
                      'usage in brackets. Taxa that cannot be validated, including new species and '
                      'other unknown taxa, morphospecies, functional groups and taxonomic levels '
                      'not used in the GBIF backbone are shown in square brackets.',
                      DIV(taxon_index_to_text(taxon_index))))

    return desc


def table_description(tab):
    """
    Function to return a description for an individual source file in a dataset.
    Typically datasets only have a single source file - the Excel workbook that
    also contains the metadata - but they may also report on external files loaded
    directly to Zenodo, and which uses the same mechanism.

    Args:
        tab: A dict describing a data table

    Returns:
        A gluon object containing an HTML description of the table
    """

    # table summary
    tab_desc = CAT(P(B(tab['title']), ' (described in worksheet ', tab['name'], ')'),
                   P('Description: ', tab['description']),
                   P('Number of fields: ', tab['max_col'] - 1))

    # The explicit n_data_row key isn't available for older records
    if 'n_data_row' in tab:
        if tab['n_data_row'] == 0:
            tab_desc += P('Number of data rows: Unavailable (table metadata description only).')
        else:
            tab_desc += P('Number of data rows: {}'.format(tab['n_data_row']))
    else:
        tab_desc += P('Number of data rows: {}'.format(tab['max_row'] - len(tab['descriptors'])))

    # add fields
    tab_desc += P('Fields: ')

    # fields summary
    flds = UL()
    for each_fld in tab['fields']:
        flds.append(LI(B(each_fld['field_name']),
                       ': ', each_fld['description'],
                       ' (Field type: ', each_fld['field_type'], ')'))

    return tab_desc + flds


def generate_inspire_xml(record):
    """
    Produces an INSPIRE/GEMINI formatted XML record from a published
    dataset record, using a template XML file stored in the static files
    """

    # get the dataset and zenodo metadata
    dataset_md = record.dataset_metadata
    zenodo_md = record.zenodo_metadata

    # parse the XML template and get the namespace map
    template = os.path.join(current.request.folder, 'static', 'files', 'gemini_xml_template.xml')
    tree = etree.parse(template)
    root = tree.getroot()
    nsmap = root.nsmap

    # Use find and XPATH to populate the template, working through from the top of the file

    # file identifier
    root.find('./gmd:fileIdentifier/gco:CharacterString',
              nsmap).text = 'zenodo.' + str(record.zenodo_record_id)

    # date stamp (not clear what this is - taken as publication date)
    root.find('./gmd:dateStamp/gco:DateTime',
              nsmap).text = record.publication_date.isoformat()

    # Now zoom to the data identication section
    data_id = root.find('.//gmd:MD_DataIdentification', nsmap)

    # CITATION
    citation = data_id.find('gmd:citation/gmd:CI_Citation', nsmap)
    citation.find('gmd:title/gco:CharacterString',
                  nsmap).text = dataset_md['title']
    citation.find('gmd:date/gmd:CI_Date/gmd:date/gco:Date',
                  nsmap).text = record.publication_date.date().isoformat()

    # two identifiers - the safe project website and the DOI.
    safe_url = URL('datasets', 'view_dataset', vars={'id': record.zenodo_record_id}, scheme=True, host=True)
    citation.find('gmd:identifier/gmd:MD_Identifier/gmd:code/gco:CharacterString',
                  nsmap).text = safe_url
    citation.find('gmd:identifier/gmd:RS_Identifier/gmd:code/gco:CharacterString',
                  nsmap).text = record.zenodo_record_doi

    # The citation string
    authors = [au['name'] for au in dataset_md['authors']]
    author_string = ', '.join(authors)
    if len(authors) > 1:
        author_string = author_string.replace(', ' + authors[-1], ' & ' + authors[-1])

    cite_string = '{} ({}) {} [Dataset] {}'.format(author_string,
                                                   record.publication_date.year,
                                                   record.dataset_title,
                                                   record.zenodo_record_doi)

    citation.find('gmd:otherCitationDetails/gco:CharacterString', nsmap).text = cite_string

    # ABSTRACT
    data_id.find('gmd:abstract/gco:CharacterString', nsmap).text = dataset_md['description']

    # KEYWORDS
    # - find the container node for the free keywords
    keywords = data_id.find('./gmd:descriptiveKeywords/gmd:MD_Keywords', nsmap)
    # - get the placeholder node
    keywd_node = keywords.getchildren()[0]
    # - duplicate it if needed
    for new_keywd in range(len(dataset_md['keywords']) - 1):
        keywords.append(copy.deepcopy(keywd_node))
    # populate the nodes
    for key_node, val in zip(keywords.getchildren(), dataset_md['keywords']):
        key_node.find('./gco:CharacterString', nsmap).text = val

    # AUTHORS - find the point of contact with author role from the template and its index
    # using xpath() here to access full xpath predicate search.
    au_xpath = "./gmd:pointOfContact[gmd:CI_ResponsibleParty/gmd:role/gmd:CI_RoleCode='author']"
    au_node = data_id.xpath(au_xpath, namespaces=nsmap)[0]
    au_idx = data_id.index(au_node)

    # - duplicate it if needed into the tree
    for n in range(len(dataset_md['authors']) - 1):
        data_id.insert(au_idx, copy.deepcopy(au_node))

    # now populate the author nodes, there should now be one for each author
    au_ls_xpath = "./gmd:pointOfContact[gmd:CI_ResponsibleParty/gmd:role/gmd:CI_RoleCode='author']"
    au_node_list = data_id.xpath(au_ls_xpath, namespaces=nsmap)

    for au_data, au_node in zip(dataset_md['authors'], au_node_list):
        resp_party = au_node.find('gmd:CI_ResponsibleParty', nsmap)
        resp_party.find('gmd:individualName/gco:CharacterString',
                        nsmap).text = au_data['name']
        resp_party.find('gmd:organisationName/gco:CharacterString',
                        nsmap).text = au_data['affiliation']
        contact_info = resp_party.find('gmd:contactInfo/gmd:CI_Contact', nsmap)
        email_path = 'gmd:address/gmd:CI_Address/gmd:electronicMailAddress/gco:CharacterString'
        contact_info.find(email_path, nsmap).text = au_data['email']

        # handle orcid resource
        orcid = contact_info.find('gmd:onlineResource', nsmap)
        if au_data['orcid'] is None:
            contact_info.remove(orcid)
        else:
            orcid.find('gmd:CI_OnlineResource/gmd:linkage/gmd:URL',
                       nsmap).text = 'http://orcid.org/' + au_data['orcid']

    # CONSTRAINTS
    # update the citation information in the second md constraint
    md_path = 'gmd:resourceConstraints/gmd:MD_Constraints/gmd:useLimitation/gco:CharacterString'
    md_constraint = data_id.find(md_path, nsmap)
    md_constraint.text += cite_string

    # embargo or not?
    embargo_path = ('gmd:resourceConstraints/gmd:MD_LegalConstraints/'
                    'gmd:otherConstraints/gco:CharacterString')
    if dataset_md['access'] == 'embargo':
        data_id.find(embargo_path, nsmap).text = ('This data is under embargo until {}. After '
                                                  'that date there are no restrictions to public '
                                                  'access.').format(dataset_md['embargo_date'])
    elif dataset_md['access'] == 'closed':
        data_id.find(embargo_path, nsmap).text = ('This dataset is currently not publicly '
                                                  'available, please contact the authors to '
                                                  'request access.')
    else:
        data_id.find(embargo_path, nsmap).text = 'There are no restrictions to public access.'

    # EXTENTS
    temp_extent = root.find('.//gmd:EX_TemporalExtent', nsmap)
    temp_extent.find('.//gml:beginPosition', nsmap).text = dataset_md['temporal_extent'][0][:10]
    temp_extent.find('.//gml:endPosition', nsmap).text = dataset_md['temporal_extent'][1][:10]

    geo_extent = root.find('.//gmd:EX_GeographicBoundingBox', nsmap)
    geo_extent.find('./gmd:westBoundLongitude/gco:Decimal',
                    nsmap).text = str(dataset_md['longitudinal_extent'][0])
    geo_extent.find('./gmd:eastBoundLongitude/gco:Decimal',
                    nsmap).text = str(dataset_md['longitudinal_extent'][1])
    geo_extent.find('./gmd:southBoundLatitude/gco:Decimal',
                    nsmap).text = str(dataset_md['latitudinal_extent'][0])
    geo_extent.find('./gmd:northBoundLatitude/gco:Decimal',
                    nsmap).text = str(dataset_md['latitudinal_extent'][1])

    # Dataset transfer options: direct download and dataset view on SAFE website
    distrib = root.find('gmd:distributionInfo/gmd:MD_Distribution', nsmap)
    distrib.find(('gmd:transferOptions[1]/gmd:MD_DigitalTransferOptions/gmd:onLine/'
                  'gmd:CI_OnlineResource/gmd:linkage/gmd:URL'),
                 nsmap).text = zenodo_md['files'][0]['links']['download']
    distrib.find(('gmd:transferOptions[2]/gmd:MD_DigitalTransferOptions/gmd:onLine/'
                  'gmd:CI_OnlineResource/gmd:linkage/gmd:URL'), nsmap).text += str(record.zenodo_record_id)

    # LINEAGE STATEMENT
    lineage = ("This dataset was collected as part of a research project based at The"
               " SAFE Project. For details of the project and data collection, see the "
               "methods information contained within the datafile and the project "
               "website: ") + URL('projects', 'view_project', args=record.project_id,
                                  scheme=True, host=True)
    root.find(('gmd:dataQualityInfo/gmd:DQ_DataQuality/gmd:lineage/gmd:LI_Lineage/'
               'gmd:statement/gco:CharacterString'), nsmap).text = lineage

    # return the string contents
    return etree.tostring(tree)


def download_ris_data(resources: Resources = None, ris_file: str = None) -> list:
    """Downloads SAFE records into a RIS format bibliography file

    This function is used to maintain a bibliography file of the records
    uploaded to a safedata community on Zenodo. It accesses the Zenodo community
    specified in the resource configuration and downloads all records. It then
    optinally checks the list of downloaded DOIs against the content of an
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

    if os.path.exists(ris_file):
        with open(ris_file, 'r') as bibliography_file:
            entries = rispy.load(bibliography_file)
            for entry in entries:
                record_id = int(entry['url'].split('/')[-1])
                known_recids.append(record_id)

    # Zenodo API call to return the records associated with the SAFE community
    z_api, _, z_cname = _resources_to_zenodo_api(resources)
    api = f'{z_api}/records/?q=communities:{z_cname}'

    # Provide feedback on DOI collection
    LOGGER.info(f'Fetching record DOIs from {api}:')
    FORMATTER.push()

    # The API is paged - it contains a set of records and a link that points 
    # to the next page of records, so keep looping until there are no more next
    n_records = 0
    while True:

        # Get the data
        safe_data = requests.get(api)

        if safe_data.status_code != 200:
            raise IOError('Cannot access Zenodo API')
        else:
            # Retrieve the record data and store the DOI for each record
            safe_data = safe_data.json()
            for hit in safe_data['hits']['hits']:
                if hit['id'] not in known_recids:
                    new_doi.append(hit['doi'])
            
            # Reporting
            n_records += len(safe_data['hits']['hits'])
            LOGGER.info(f'{n_records}')

            # Update the link for the next page, unless there is no next page
            if 'next' in safe_data['links']:
                api = safe_data['links']['next']
            else:
                break

    # Use the datacite API to retrieve the citation data associated with the DOI
    # and save it out to a RIS format file
    if not new_doi:
        LOGGER.info('No new DOIs found')
        return

    # Get the DOI data    
    data = []

    FORMATTER.pop()
    LOGGER.info(f'Retrieving citation data from Datacite for {len(new_doi)} new records')
    FORMATTER.push()

    for doi in new_doi:
        
        ris_data = requests.get(f'https://data.datacite.org/application/x-research-info-systems/{doi}')
                
        if ris_data.status_code != 200:
            LOGGER.warning(f'DOI {doi} not found in datacite.org')
        else:
            # Write the response content to the data list. It comes in as byte
            # data so needs to be decoded to a string variable
            LOGGER.info(f'Retrieved citation for DOI {doi}')
            data.append(ris_data.content.decode("utf-8") + '\r\n')

    FORMATTER.pop()

    if os.path.exists(ris_file):
        LOGGER.info(f'Appending RIS data for {len(data)} new records to {ris_file}')
        write_mode = 'a'
    else:
        LOGGER.info(f'Writing RIS data for {len(data)} records to {ris_file}')
        write_mode = 'w'

    with open(ris_file, write_mode) as ris_file:
        for this_entry in data:
            ris_file.write(this_entry)




def sync_local_dir(datadir: str, api: str = None, 
                   xlsx_only: bool = True, resources: Resources = None) -> None:

    """
    The safedata R package defines a directory structure used to store metadata
    and files downloaded from a safedata community on Zenodo and from a safedata
    metadata server. This tool allows a safedata developer or community
    maintainer to create or update such a directory with _all_ of the resources
    in the Zenodo community, regardless of their public access status. This
    forms a backup (although Zenodo is heavily backed up) but also provides
    local copies of the files for testing and development of the code packages.

    You need to provide a Zenodo API token to use this script. That is obtained
    by logging into the Zenodo account managing the community and going to the
    Applications tab and creating a personal access token. These allow root
    level access to the community files so must be treated carefully!

    If this is a new data directory, you also need to provide the API url for
    the safedata metadata server.

    Args:
        datadir: The path to a local directory containing an existing safedata
            directory or an empty folder in which to create one.
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.
        api: An API from which JSON dataset metadata can be downloaded. If 
            datadir is an existing safedata directory, then the API will be read
            from `url.json`.
        xlsx_only: Should the download ignore large non-xlsx files, defaulting
            to True.
    """

    # Private helper functions
    def _get_file(url: str, outf:str, params: dict=None) -> None:
        """Download a file from a URL
        """
        resource = requests.get(url, params=params, stream=True)

        with open(outf, 'wb') as outf_obj:
            shutil.copyfileobj(resource.raw, outf_obj)

    # Get resource configuration
    zenodo_api, params, _ = _resources_to_zenodo_api(resources)

    # The dir argument should be an existing path
    if not (os.path.exists(datadir) and os.path.isdir(datadir)):
        raise IOError(f'{datadir} is not an existing directory')

    # Check for an existing API url file and then see what is provided and resolve
    url_file = os.path.join(datadir, 'url.json')

    if os.path.exists(url_file):
        with open(url_file, 'r') as urlf:
            dir_api = simplejson.load(urlf)['url'][0]
    else:
        dir_api = None

    if api is None and dir_api is None:
        raise RuntimeError('API not provided or found in directory')

    if api is not None and dir_api is not None and api != dir_api:
        raise RuntimeError('Provided api does not match existing api in directory')
    
    if api is not None and dir_api is None:
        with open(url_file, 'w') as urlf:
            simplejson.dump({'url': [api]}, urlf)
    else:
        api = dir_api
    
    # Download index files - don't bother to check for updates, this isn't
    # a frequent thing to do
    LOGGER.info("Downloading index files")
    _get_file(f'{api}/api/index', os.path.join(datadir, 'index.json'))
    _get_file(f'{api}/api/gazetteer', os.path.join(datadir, 'gazetteer.geojson'))
    _get_file(f'{api}/api/location_aliases', os.path.join(datadir, 'location_aliases.csv'))

    # Get the deposits associated with the account, which includes a list of download links
    params['page'] = 1
    deposits = []

    LOGGER.info("Scanning Zenodo deposits")
    while True:
        this_page = requests.get(zenodo_api + '/deposit/depositions',
                                 params=params,
                                 json={},
                                 headers={"Content-Type": "application/json"})

        if not this_page.ok:
            raise RuntimeError('Could not connect to Zenodo API. Invalid token?')

        if this_page.json():
            deposits += this_page.json()
            print(f" - Page {params['page']}")
            params['page'] += 1
        else:
            break

    LOGGER.info(f"Processing {len(deposits)} deposits")

    # Download the files
    for dep in deposits:

        con_rec_id = str(dep['conceptrecid'])
        rec_id = str(dep['record_id'])

        if not dep['submitted']:
            LOGGER.info(f'Unsubmitted draft {con_rec_id}/{rec_id}')
            continue

        LOGGER.info(f'Processing deposit {con_rec_id}/{rec_id}')
        FORMATTER.push()

        # Create the directory structure if needed
        rec_dir = os.path.join(datadir, con_rec_id, rec_id)
        if not os.path.exists(rec_dir):
            LOGGER.info('Creating directory')
            os.makedirs(rec_dir)
        else:
            LOGGER.info('Directory found')

        # loop over the files in the record
        for this_file in dep['files']:

            if xlsx_only and not this_file['filename'].endswith('.xlsx'):
                LOGGER.info(f"Skipping non-excel file {this_file['filename']}")
                continue

            LOGGER.info(f"Processing {this_file['filename']}")
            FORMATTER.push()

            outf = os.path.join(rec_dir, this_file['filename'])
            local_copy = os.path.exists(outf)

            if not local_copy:
                LOGGER.info("Downloading")
                _get_file(this_file['links']['download'], outf, params=params)
            elif local_copy and _compute_md5(outf) != this_file['checksum']:
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
            _get_file(f'{api}/record/{rec_id}', metadata)
        
        FORMATTER.pop()

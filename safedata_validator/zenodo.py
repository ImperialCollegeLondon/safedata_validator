import os
import datetime
from lxml import etree
import simplejson
import copy
import safedata_validator
from io import StringIO
import requests
from safe_web_global_functions import safe_mailer
from itertools import groupby
from shapely import geometry
import hashlib
from gluon.serializers import json

# The web2py HTML helpers are provided by gluon. This also provides the 'current' object, which
# provides the web2py 'request' API (note the single letter difference from the requests package!).
# The 'current' object is also extended by models/db.py to include the current 'db' DAL object
# and the 'myconf' AppConfig object so that they can accessed by this module

from gluon import *

"""
This module provides functions to handle datasets within the SAFE website. 
These are called from the datasets controller but are also needed from other locations, 
such as the scheduler, so are defined here in their own module.
"""


def verify_dataset(record_id, email=False):
    """
    Function to run safedata_validator on an uploaded file. There
    are three possible outcomes for a dataset: PASS; FAIL, if the check
    catches known formatting problems; and ERROR if check hits an exception,
    which probably means an update to the checker code to handle the new
    and exciting way of getting the file wrong.

    This assumes that the local python has safedata_validator installed and
    correctly configured with local file resources. See the documentation
    for that package for details.

    The email argument allows the function to be run by admins directly
    using the administer_datasets controller without spamming the uploader.
    Useful in error checking problem uploads.

    This function is also used as a scheduler task, so needs some extra
    care to make sure it runs in the environment of a scheduler worker
    as well as in the environment of the website.

    Args:
        record_id: The id of the record from the dataset table that is to be checked.
        email: Should the dataset uploader be emailed the outcome?

    Returns:
        A string describing the outcome of the check that gets stored in the
        scheduler results or sent back to the administer_datasets controller.
        This is primarily a user friendly bit of text for popping up in a
        website flash, not an Exception message, so those are stored elsewhere
        for an admin to look at.
    """

    # Load the host name from the configuration. When run from a controller,
    # the URL(host=TRUE) has access to the host name from requests. This isn't
    # true when it is run by a scheduler worker, which isn't operating as part
    # of the website. So rather than hardcoding, store the host name in the config
    try:
        host = current.myconf.take('host.host_name')
    except BaseException:
        raise RuntimeError('Site config does provide not the host name')

    # get the record
    record = current.db.submitted_datasets[record_id]

    # track errors to avoid hideous nested try statements
    error = False

    if record is None:
        # not a valid record? Can't email anyone so turn that off
        ret_msg = 'Verifying dataset {}: unknown record ID'.format(record_id)
        email = False
        error = True
    else:
        # otherwise, create a return dictionary for all remaining failure
        # modes (no report, but file, uploader and URL should fine) and
        # set the default outcome
        ret_dict = {'report': '',
                    'filename': record.file_name,
                    'name': record.uploader_id.first_name,
                    'dataset_url': URL('datasets', 'submitted_dataset_status',
                                       vars={'id': record.id},
                                       scheme=True, host=host)}
        outcome = 'ERROR'

    # Initialise the dataset checker:
    if not error:
        # - get paths to dataset file. Failure to find is handled by safedata_validator methods.
        fname = os.path.join(current.request.folder, 'uploads', 'submitted_datasets', record.file)
        # get the Dataset object from the file checker
        try:
            dataset = safedata_validator.Dataset(fname, verbose=False)
        except Exception as e:
            # We don't want to bail here because we might want to email the uploader,
            # but we do want to record what went wrong. We store it in the dataset record, which
            # is the only venue when run from a controller. If I could work out where the scheduler
            # run output comes from, I'd do that too.
            record.update_record(dataset_check_outcome='ERROR',
                                 dataset_check_error=repr(e))
            ret_msg = 'Verifying dataset {}: error initialising dataset checker'
            ret_msg = ret_msg.format(record.id)
            error = True

    # main processing of the dataset
    if not error:
        try:
            # load the metadata sheets
            dataset.load_summary(validate_doi=True, project_id=record.project_id)

            dataset.load_taxa()
            dataset.load_locations()

            # check the dataset worksheets, after checking there are some
            # as the metadata might only document external data files
            if dataset.dataworksheet_summaries is not None:
                for ws in dataset.dataworksheet_summaries:
                    dataset.load_data_worksheet(ws)

            # cross check the taxa and locations
            dataset.final_checks()

        except Exception as e:
            ret_msg = 'Verifying dataset {}: error running dataset checking'
            ret_msg = ret_msg.format(record.id)
            dataset_check_error = repr(e)
        else:
            if dataset.passed:
                outcome = 'PASS'
                ret_msg = 'Verifying dataset {}: dataset checking PASSED'.format(record.id)
            else:
                outcome = 'FAIL'
                ret_msg = 'Verifying dataset {}: dataset checking FAILED'.format(record.id)

            dataset_check_error = ''

        # At this point, we have a Dataset object, so can populate the record with
        # what information is available, regardless of Error, Fail or Pass
        # - The DAL handles conversion of the python structure into JSON, so
        #   can just pass the objects. Previously, used simplejson.dumps, which
        #   meant they were stored as a string, so needed reloading rather than
        #   being saved natively as JSON

        # First, need to extract the check report from the StringIO object and
        # substitute in the user filename for the local web2py filename. Also,
        # wrap it in <pre> for display purposes.
        if outcome == 'ERROR':
            report_text = ""
        else:
            # TODO - something about the safedata_validator logger setup is causing
            #        large numbers of null characters to be included in the report
            #        text, which then prevents things from being written to the DB
            #        The replace here is a hack to clear them out, but this needs
            #        to be fixed upstream.
            report_text = dataset.report().getvalue().replace('\x00', '')
            report_text = PRE(report_text.replace(fname, record.file_name))
            # Update the ret_dict to insert the report text
            ret_dict['report'] = report_text

        # get the dataset_metadata and update to store the actual filename
        dataset_metadata = dataset.export_metadata_dict()
        dataset_metadata['filename'] = record.file_name

        # In the submitted datasets table, taxa, locations and metadata are
        # all stored in dataset_metadata. When a dataset is published, taxa
        # and locations are separated out into the dataset_* tables.
        record.update_record(dataset_check_outcome=outcome,
                             dataset_check_report=report_text,
                             dataset_check_error=dataset_check_error,
                             dataset_title=dataset.title,
                             dataset_metadata={'metadata': dataset_metadata,
                                               'taxa': dataset.taxon_index,
                                               'locations': dataset.location_index})

    # notify the user
    if email:
        opts = {'PASS': ['Dataset passed checks', 'dataset_check_pass.html'],
                'FAIL': ['Dataset failed checks', 'dataset_check_fail.html'],
                'ERROR': ['Error checking dataset', 'dataset_check_error.html']}

        safe_mailer(to=record.uploader_id.email,
                    cc=['data@safeproject.net'],
                    reply_to='data@safeproject.net',
                    cc_info=False,
                    subject=opts[outcome][0],
                    template=opts[outcome][1],
                    template_dict=ret_dict)

    # A task run by a worker does not automatically commit changes, so
    # save any by changes before ending
    current.db.commit()

    return ret_msg


def get_zenodo_api():
    """
    Function to provide the zenodo API endpoint and the access token
    from the site config. The config specifies whether the testing
    sandbox or the live site is to be used.
    """

    try:
        sandbox = int(current.myconf.take('zenodo.use_sandbox'))
    except BaseException:
        raise RuntimeError('Site config does not provide zenodo.use_sandbox')

    if sandbox:
        try:
            token = {'access_token': current.myconf.take('zenodo.sandbox_access_token')}
        except BaseException:
            raise RuntimeError('Site config does not provide zenodo.sandbox_access_token')

        api = 'https://sandbox.zenodo.org/api/'
    else:
        try:
            token = {'access_token': current.myconf.take('zenodo.access_token')}
        except BaseException:
            raise RuntimeError('Site config does not provide zenodo.access_token')

        api = 'https://zenodo.org/api/'

    return api, token


def submit_dataset_to_zenodo(record_id, deposit_id=None):
    """
    Function that attempts to publish a dataset record to Zenodo and
    to populate the published_datasets table with result of that attempt.
    This handles the logic of selecting which method to use: create excel,
    update excel or adopt external.

    Args:
        record_id: The id of the dataset table record to be submitted
        deposit_id: An integer giving the id of an existing Zenodo deposit to adopt
            using this dataset record.
    Returns:
        A string describing the outcome.
    """

    # get the current db
    db = current.db

    # check the record exists and hasn't already been submitted
    record = db.submitted_datasets[record_id]

    if record is None:
        return 'Publishing dataset: unknown record ID {}'.format(record_id)
    elif record.dataset_check_outcome != 'PASS':
        return 'Publishing dataset: record ID {} has not passed format checking'.format(record_id)
    elif record.zenodo_submission_status == 'ZEN_PASS':
        return 'Publishing dataset: record ID {} already published'.format(record_id)

    # load the correct API and token
    api, token = get_zenodo_api()

    # There are then four possible options of things that could be published:
    # 1) a brand new excel-only dataset,
    # 2) an update to an existing excel-only dataset,
    # 3) a brand new dataset with external files and
    # 4) an update to an existing dataset with online files.

    metadata = record.dataset_metadata['metadata']

    # external_files contains an empty list or a list of dictionaries
    if metadata['external_files']:
        code, links, response = adopt_external_zenodo(api, token, record, deposit_id)
        external = True
    else:
        if record.concept_id is None:
            code, links, response = create_excel_zenodo(api, token, record)
        else:
            code, links, response = update_excel_zenodo(api, token, record)
        external = False

    if code > 0:
        # There has been a problem. If this is an internal Excel file only, then try
        # and delete the failed deposit and update the record
        if links is not None and not external:
            # This can fail and leave a hanging deposit, but we won't let that stop the function
            _, _ = delete_deposit(links, token)

        # update the record
        record.update_record(zenodo_submission_status='ZEN_FAIL',
                             zenodo_submission_date=datetime.datetime.now(),
                             zenodo_error=response)
        return "Failed to publish record"
    else:

        # Set the most recent flag for existing published versions to False
        if record.concept_id is not None:
            db(db.published_datasets.zenodo_concept_id == record.concept_id
               ).update(most_recent=False)

        # remove the dataset metadata from the zenodo response, since the
        # contents is information we have already, so we can store the rest
        del response['metadata']

        # Now create a published datasets entry
        published_record = db.published_datasets.insert(
            uploader_id=record.uploader_id,
            upload_datetime=record.upload_datetime,
            submission_id=record.id,
            dataset_title=metadata['title'],
            dataset_access=metadata['access'],
            dataset_embargo=metadata['embargo_date'],
            dataset_conditions=metadata['access_conditions'],
            dataset_description=metadata['description'],
            dataset_metadata=record.dataset_metadata['metadata'],
            temporal_extent_start=metadata['temporal_extent'][0],
            temporal_extent_end=metadata['temporal_extent'][1],
            geographic_extent=geometry.box(metadata['longitudinal_extent'][0],
                                           metadata['latitudinal_extent'][0],
                                           metadata['longitudinal_extent'][1],
                                           metadata['latitudinal_extent'][1]).wkt,
            publication_date=datetime.datetime.now(),
            most_recent=True,
            zenodo_record_id=response['record_id'],
            zenodo_record_doi=response['doi_url'],
            zenodo_record_badge=response['links']['badge'],
            zenodo_concept_id=response['conceptrecid'],
            zenodo_concept_doi=response['links']['conceptdoi'],
            zenodo_concept_badge=response['links']['conceptbadge'],
            zenodo_metadata=response)

        # add an associated project link to the dataset concept
        # if one does not already exist
        project_link = db((db.project_datasets.project_id == record.project_id) &
                          (db.project_datasets.concept_id == response['conceptrecid'])
                          ).select()

        if not project_link:
            db.project_datasets.insert(project_id=record.project_id,
                                       concept_id=response['conceptrecid'],
                                       user_id=record.uploader_id,
                                       date_added=datetime.date.today())

        # update to include the UTM 50N bbox geometry
        db(db.published_datasets.id == published_record).update(
            geographic_extent_utm50n=db.published_datasets.geographic_extent.st_transform(32650))

        # populate index tables
        # A) Taxa
        taxa = record.dataset_metadata['taxa']
        taxa = [dict(list(zip(['dataset_id', 'worksheet_name', 'gbif_id', 'gbif_parent_id',
                               'taxon_name', 'taxon_rank', 'gbif_status'],
                              [published_record] + tx))) for tx in taxa]

        db.dataset_taxa.bulk_insert(taxa)

        # B) Files, using the Zenodo response
        files = response['files']
        for each_file in files:
            each_file['dataset_id'] = published_record
            each_file['download_link'] = each_file['links']['download']
            each_file['file_zenodo_id'] = each_file.pop('id')

        db.dataset_files.bulk_insert(files)

        # C) Locations
        locations = record.dataset_metadata['locations']
        locations = [dict(list(zip(['dataset_id', 'name', 'new_location', 'type', 'wkt_wgs84'],
                                   [published_record] + loc))) for loc in locations]

        new_locs = db.dataset_locations.bulk_insert(locations)

        # update the UTM 50 N geometry where possible
        db(db.dataset_locations.id.belongs(new_locs)).update(
            wkt_utm50n=db.dataset_locations.wkt_wgs84.st_transform(32650))

        # D) Dataworksheets and fields
        for data in metadata['dataworksheets']:

            worksheet_id = db.dataset_worksheets.insert(dataset_id=published_record, **data)

            for fld in data['fields']:
                fld['dataset_id'] = published_record
                fld['worksheet_id'] = worksheet_id

            db.dataset_fields.bulk_insert(data['fields'])

        # E) Authors
        for auth in metadata['authors']:
            db.dataset_authors.insert(dataset_id=published_record, **auth)

        # F) Funders
        if metadata['funders'] is not None:
            for fndr in metadata['funders']:
                db.dataset_funders.insert(dataset_id=published_record, **fndr)

        # G) Permits
        if metadata['permits'] is not None:
            for perm in metadata['permits']:
                db.dataset_permits.insert(dataset_id=published_record, **perm)

        # K) Keywords
        if metadata['keywords'] is not None:
            for kywd in metadata['keywords']:
                db.dataset_keywords.insert(dataset_id=published_record, keyword=kywd)

        # remove the dataset from the submitted_datasets table
        record.delete_record()

        # Flush the cached index of published datasets
        current.cache.ram('index', None)

        return "Published dataset to {}".format(response['doi_url'])


def create_excel_zenodo(api, token, record):
    """
    A function to work through the Zenodo API steps to publish an new Excel only dataset.
    It works through the publication steps as long as each step keeps returning a zero
    success code, otherwise we get to the end with the most recent failure

    Args:
        api: The API URL to use: sandbox or main site
        token: A dictionary containing the key 'access_token'
        record: The dataset row for the record to publish
    Returns:
        i) An integer code indicating success or failure,
        ii) the links object for the deposit - which is needed to delete
            partially created deposits and
        iii) A response object from Zenodo - which will contain either a
            failure message or the publication details.
    """

    # create the new deposit
    code, response, zenodo_id = create_deposit(api, token)

    # upload the record metadata
    if code == 0:
        # store previous response containing the links dictionary
        links = response
        code, response = upload_metadata(links, token, record, zenodo_id)
    else:
        links = None

    # upload the file
    if code == 0:
        code, response = upload_file(links, token, record)

    # publish the deposit
    if code == 0:
        code, response = publish_deposit(links, token)

    # Return what we've got
    if code > 0:
        return 1, links, response
    else:
        return 0, links, response


def update_excel_zenodo(api, token, record):
    """
    A function to work through the Zenodo API steps to publish a new version of an Excel
    only dataset. It works through the publication steps as long as each step keeps returning
    a zero success code, otherwise we get to the end with the most recent failure

    Args:
        api: The API URL to use: sandbox or main site
        token: A dictionary containing the key 'access_token'
        record: The dataset row for the record to publish
    Returns:
        i) An integer code indicating success or failure,
        ii) the links object for the deposit - which is needed to delete
            partially created deposits and
        iii) A response object from Zenodo - which will contain either a
            failure message or the publication details.
    """
    # get a new draft of the existing record from the zenodo id of the
    # most recent published record (can't create a draft from a concept id)
    most_recent = current.db((current.db.published_datasets.zenodo_concept_id == record.concept_id) &
                             (current.db.published_datasets.most_recent == True)
                             ).select().first()

    code, response, zenodo_id = create_deposit_draft(api, token, most_recent.zenodo_record_id)

    # upload the record metadata
    if code == 0:
        # store previous response containing the links dictionary
        links = response
        code, response = upload_metadata(links, token, record, zenodo_id)
    else:
        links = None

    # delete the existing file
    if code == 0:
        code, response = delete_previous_file(links, token)

    # upload the new file
    if code == 0:
        code, response = upload_file(links, token, record)

    # publish the deposit
    if code == 0:
        code, response = publish_deposit(links, token)

    # Return what we've got
    if code > 0:
        return 1, links, response
    else:
        return 0, links, response


def adopt_external_zenodo(api, token, record, deposit_id):
    """
    A function to work through the Zenodo API steps to publish a dataset that adopts
    external files in an existing deposit. It works through the publication steps as
    long as each step keeps returning a zero success code, otherwise we get to the end
    with the most recent failure

    Args:
        api: The API URL to use: sandbox or main site
        token: A dictionary containing the key 'access_token'
        record: The dataset row for the record to publish
        deposit_id: An integer giving the id of an existing Zenodo deposit to adopt
            using this dataset record.
    Returns:
        i) An integer code indicating success or failure,
        ii) the links object for the deposit - which is needed to delete
            partially created deposits and
        iii) A response object from Zenodo - which will contain either a
            failure message or the publication details.
    """

    # get the deposit
    code, response = get_deposit(api, token, deposit_id)

    # upload the record metadata
    if code == 0:
        # store previous response containing the links dictionary
        # and the list of remote files
        remote_files = response['files']
        links = response['links']
        code, response = upload_metadata(links, token, record, deposit_id)

        # If we got a deposit, check the files found in the deposit match
        # with the external files specified in the record metadata.
        remote_filenames = {rfile['filename'] for rfile in remote_files}
        external_files = set([r['file'] for r in record.dataset_metadata['metadata']['external_files']])

        if not remote_filenames == external_files:
            code = 1
            response = "Files in deposit do not match external files listed in Excel file"
            links = None
    else:
        links = None

    # Upload the Excel file - the expectation here is that the Excel file
    # associated with previous drafts is deleted as part of the manual file
    # update process, so we only have to upload the one submitted to the website
    if code == 0:
        code, response = upload_file(links, token, record)

    # publish the deposit
    if code == 0:
        code, response = publish_deposit(links, token)

    # Return what we've got
    if code > 0:
        return 1, links, response
    else:
        return 0, links, response


"""
Zenodo action functions
"""


def get_deposit(api, token, deposit_id):
    # request the deposit
    dep = requests.get(api + 'deposit/depositions/{}'.format(deposit_id), params=token, json={},
                       headers={"Content-Type": "application/json"})

    # check for success and return the information.
    if dep.status_code != 200:
        return 1, dep.json()
    else:
        return 0, dep.json()


def create_deposit(api, token):
    """
    Function to create a new deposit
    Args:
        api: The api URL to be used (standard or sandbox)
        token: The access token to be usedz
    Returns:
        An integer indicating success (0) or failure (1) and either the
        deposit links dictionary or an error message
    """

    # get a new deposit resource
    dep = requests.post(api + '/deposit/depositions', params=token, json={},
                        headers={"Content-Type": "application/json"})

    # trap errors in creating the resource - successful creation of new deposits returns 201
    if dep.status_code != 201:
        return 1, dep.json(), None
    else:
        return 0, dep.json()['links'], dep.json()['id']


def create_deposit_draft(api, token, deposit_id):
    """
    Function to create a new draft of an existing published record
    Args:
        api: The api URL to be used (standard or sandbox)
        token: The access token to be used
        deposit_id: The id of the published record
    Returns:
        An integer indicating success (0) or failure (1) and either the
        deposit links dictionary for the new draft or an error message
    """

    # get the draft api
    new_draft = requests.post(api + '/deposit/depositions/{}/actions/newversion'.format(deposit_id),
                              params=token, json={},
                              headers={"Content-Type": "application/json"})

    # trap errors in creating the new version
    if new_draft.status_code != 201:
        return 1, new_draft.json(), None

    # now get the newly created version
    api = new_draft.json()['links']['latest_draft']
    dep = requests.get(api, params=token, json={},
                       headers={"Content-Type": "application/json"})

    # trap errors in creating the resource - successful creation of new version
    #  drafts returns 200
    if dep.status_code != 200:
        return 1, dep.json(), None
    else:
        return 0, dep.json()['links'], dep.json()['id']


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

    # If all steps have been succesful, return a 0 code, otherwise
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


# # Untested snippet from: https://gist.github.com/tyhoff/b757e6af83c1fd2b7b83057adf02c139
# # possibly of use when and if this module gets built into safedata_validator,
# # to provide user feedback on uploading bulk data to drafts.
#
# from tqdm import tqdm
# from tqdm.utils import CallbackIOWrapper
#
# file_path = os.path.abspath(__file__)
# upload_url = https://some-bucket.s3.amazonaws.com
#
# file_size = os.stat(file_path).st_size
# with open(file_path, "rb") as f:
#     with tqdm(total=file_size, unit="B", unit_scale=True, unit_divisor=1024) as t:
#         wrapped_file = CallbackIOWrapper(t.update, f, "read")
#         requests.put(upload_url, data=wrapped_file)

def upload_file(links, token, record):
    """
    Function to upload the Excel datafile submitted for a record to Zenodo deposit.

    Args:
        links: The links dictionary from a created deposit
        token: The access token to be used
        record: The database record containing the metadata to be uploaded.

    Returns:
        An integer indicating success (0) or failure (1) and either the
        deposit links dictionary or an error message
    """

    # upload the new file
    bucket_url = links['bucket']
    fname = os.path.join(current.request.folder, 'uploads', 'submitted_datasets', record.file)

    with open(fname, 'rb') as fp:
        fls = requests.put(bucket_url + '/' + record.file_name,
                           data=fp,
                           params=token)

    # trap errors in uploading file
    # - no success or mismatch in md5 checksums
    if fls.status_code != 200:
        return 1, fls.json()
    elif fls.json()['checksum'] != 'md5:' + record.file_hash:
        return 1, "Mismatch in local and uploaded MD5 hashes"
    else:
        return 0, 'success'


def delete_deposit(links, token):
    """
    Function to delete an (unpublished) partially created deposit if the publication
    process fails.

    Args:
        links: The links dictionary from a created deposit
        token: The access token to be used

    Returns:
        An integer indicating success (0) or failure (1) and either the
        deposit links dictionary or an error message
    """

    delete = requests.delete(links['self'], params=token)

    if delete.status_code != 204:
        return 1, delete.json()
    else:
        return 0, 'success'


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


def delete_previous_file(links, token):
    """
    Function to delete a previously uploaded file from a new version of a deposit,
    prior to replacing it with an updated one.

    Args:
        links: The links dictionary from a created deposit
        token: The access token to be used

    Returns:
        An integer indicating success (0) or failure (1) and either the
        deposit links dictionary or an error message
    """

    # get the existing files
    files = requests.get(links['files'], params=token)

    # check the result of the files request
    if files.status_code != 200:
        # failed to get the files
        return 1, files.json()
    elif len(files.json()) != 1:
        # multiple files
        return 1, files.json()

    # get the delete link to the file and call
    delete_api = files.json()[0]['links']['self']
    file_del = requests.delete(delete_api, params=token)

    if file_del.status_code != 204:
        return 1, file_del.json()
    else:
        return 0, 'success'


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


# API search functions

def dataset_query_to_json(qry, most_recent=False, ids=None,
                          fields=[('published_datasets', 'zenodo_concept_id'),
                                  ('published_datasets', 'zenodo_record_id'),
                                  ('published_datasets', 'dataset_title')]):
    """
    Shared function to take a Query including rows in db.published datasets
    and return a standardised set of attributes and a count.
    """

    db = current.db

    if most_recent:
        qry &= (db.published_datasets.most_recent == True)

    if ids is not None:
        qry &= (db.published_datasets.zenodo_record_id.belongs(ids))

    # Turn fields argument into fields references and select
    fields = [db[t][f] for t, f in fields]
    rows = db(qry).select(*fields, distinct=True)

    return {'count': len(rows), 'entries': rows}


def dataset_taxon_search(gbif_id=None, name=None, rank=None):
    """Search for datasets by taxon information

    Examples:
        /api/search/taxa?name=Formicidae
        /api/search/taxa?gbif_id=4342
        /api/search/taxa?rank=Family

    Args:
        gbif_id (int): A GBIF taxon id.
        name (str): A scientific name
        rank (str): A taxonomic rank. Note that GBIF only provides
            kingdom, phylum, order, class, family, genus and species.
    """

    db = current.db
    qry = (db.published_datasets.id == db.dataset_taxa.dataset_id)

    if gbif_id is not None:
        qry &= (db.dataset_taxa.gbif_id == gbif_id)

    if name is not None:
        qry &= (db.dataset_taxa.taxon_name == name)

    if rank is not None:
        qry &= (db.dataset_taxa.taxon_rank == rank.lower())

    return qry


def dataset_author_search(name=None):
    """Search for datasets by author name

    Examples:
        /api/search/authors?name=Wilk

    Args:
        name (str): An author name or part of a name
    """

    db = current.db
    qry = (db.published_datasets.id == db.dataset_authors.dataset_id)

    if name is not None:
        qry &= (db.dataset_authors.name.contains(name.lower()))

    return qry


def dataset_locations_search(name=None):
    """Search for datasets with data at a named location

    Examples:
        /api/search/location?name=A_1

    Args:
        name (str): A location name
    """

    db = current.db
    qry = (db.published_datasets.id == db.dataset_locations.dataset_id)

    if name is not None:
        qry &= (db.dataset_locations.name.contains(name.lower()))

    return qry


def dataset_date_search(date=None, match_type='intersect'):
    """Search for datasets by temporal extent

    Examples:
        /api/search/dates?date=2014-06-12
        /api/search/dates?date=2014-06-12,2015-06-12
        /api/search/dates?date=2014-06-12,2015-06-12&match_type=contains
        /api/search/dates?date=2014-06-12,2015-06-12&match_type=within

    Args:
        date (str): A string containing one or two (comma separated) dates in ISO format (2019-06-12)
        match_type (str): One of 'intersect', 'contain' and 'within' to match the provided dates to
            the temporal extents of datasets. The 'contain' option returns datasets that span a date
            range and 'within' returns datasets that fall within that date range.
    """

    db = current.db

    # parse the data query: can be a single iso date or a comma separated pair.
    try:
        date_vals = date.split(',')
        date_vals = [datetime.datetime.strptime(dt, '%Y-%m-%d').date() for dt in date_vals]
    except ValueError:
        return {'error': 400, 'message': 'Could not parse dates: {}'.format(date)}

    # Check the number of dates and enforce order
    if len(date_vals) > 2:
        return {'error': 400, 'message': 'date contains more than two values: {}'.format(date)}
    elif len(date_vals) == 2:
        date_vals.sort()
    else:
        # make singular dates a 'range' to use same predicate mechanisms
        date_vals += date_vals

    # check the dates against the temporal extents of the datasets
    if match_type == 'intersect':
        qry = (~ ((db.published_datasets.temporal_extent_start >= date_vals[1]) |
                  (db.published_datasets.temporal_extent_end <= date_vals[0])))
    elif match_type == 'contain':
        qry = ((db.published_datasets.temporal_extent_start >= date_vals[0]) &
               (db.published_datasets.temporal_extent_end <= date_vals[1]))
    elif match_type == 'within':
        qry = ((db.published_datasets.temporal_extent_start <= date_vals[0]) &
               (db.published_datasets.temporal_extent_end >= date_vals[1]))
    else:
        return {'error': 400, 'message': 'Unknown date match type: {}'.format(match_type)}

    return qry


def dataset_field_search(text=None, ftype=None):
    """Search for datasets by data field information

    Examples:
        /api/search/fields?text=temperature
        /api/search/fields?ftype=numeric
        /api/search/fields?text=temperature&ftype=numeric

    Args:
        text (str): A string to look for within the field name and description.
        ftype (str): A field type to match.
    """

    db = current.db
    qry = (db.published_datasets.id == db.dataset_fields.dataset_id)

    if text is not None:
        qry &= ((db.dataset_fields.field_name.contains(text)) |
                (db.dataset_fields.description.contains(text)))

    if ftype is not None:
        qry &= (db.dataset_fields.field_type.ilike(ftype))

    return qry


def dataset_text_search(text=None):
    """Search for datasets by free text search

    Examples:
        /api/search/text?text=humus

    Args:
        text (str): A string to look within dataset, worksheet and field
        descriptions and titles and in dataset keywords.
    """

    db = current.db

    qry = ((db.published_datasets.id == db.dataset_fields.dataset_id) &
           (db.published_datasets.id == db.dataset_worksheets.dataset_id) &
           (db.published_datasets.id == db.dataset_keywords.dataset_id) &
           ((db.dataset_fields.field_name.contains(text)) |
            (db.dataset_fields.description.contains(text)) |
            (db.dataset_worksheets.title.contains(text)) |
            (db.dataset_worksheets.description.contains(text)) |
            (db.published_datasets.dataset_title.contains(text)) |
            (db.published_datasets.dataset_description.contains(text)) |
            (db.dataset_keywords.keyword.contains(text))
            ))

    return qry


def dataset_parse_spatial(wkt=None, location=None):
    """
    Shared function to parse query geometry options - either a location or a WKT - and get
    # the query geometry as a UTM 50N geometry.
    """
    db = current.db

    if (location is not None) and (wkt is not None):
        return {'error': 400, 'message': 'Provide a location name or a WKT geometry, not both'}
    elif location is None and wkt is None:
        return {'error': 400, 'message': 'Provide either a location name or a WKT geometry'}
    elif location is not None:
        gazetteer_location = gazetteer_location = db(db.gazetteer.location == location).select(
            db.gazetteer.wkt_utm50n.st_aswkb().with_alias('wkb')).first()
        if gazetteer_location is not None:
            query_geom = gazetteer_location.wkb
        else:
            return {'error': 400, 'message': "Unknown location"}
    elif wkt is not None:
        # Validate the geometry - there isn't currently a validator, so use the DB (?shapely)
        # i) Does the WKT parse correctly to a WKB string, assuming WGS84 for now
        try:
            query_geom = db.executesql("select st_geomfromtext('{}', 4326);".format(wkt))[0][0]
        except (db._adapter.driver.ProgrammingError, db._adapter.driver.InternalError):
            return {'error': 400, 'message': "Could not parse WKT geometry"}

        # ii) Do the coordinates seem like lat long?
        is_lat_long = db.executesql("SELECT ST_XMin(g) >= -180 AND ST_XMax(g) <= 180 AND "
                                    "   ST_YMin(g) >= -90 AND ST_YMax(g) <= 90 "
                                    "   FROM (SELECT st_geomfromwkb(decode('{0}', 'hex')) "
                                    "AS g) AS sel ;".format(query_geom))[0][0]
        if not is_lat_long:
            return {'error': 400, 'message': "WKT geometry coordinates not as lat/long"}

        # iii) Convert to UTM50N
        query_geom = \
        db.executesql("SELECT st_transform(st_geomfromwkb(decode('{0}', 'hex')), 32650);".format(query_geom))[0][0]

    return query_geom


def dataset_spatial_search(wkt=None, location=None, distance=0):
    """Spatial search for sampling locations.

    This endpoint can search for datasets using either a user-provided geometry or the geometry
    of a named location from the SAFE gazetteer. The sampling locations provided in each dataset
    are tested to see if they intersect the search geometry and a buffer distance can also be
    provided to search around the query geometry.

    Note that this endpoint will not retrieve datasets that have not provided sampling locations
    or use new locations that are missing coordinate information. The bounding box endpoint
    uses the dataset geographic extent, which is provided for all datasets.

    Examples:
        /api/search/spatial?location=A_1
        /api/search/spatial?location=A_1&distance=50
        /api/search/spatial?wkt=Point(116.5 4.75)
        /api/search/spatial?wkt=Point(116.5 4.75)&distance=50000
        /api/search/spatial?wkt=Polygon((110 0, 110 10,120 10,120 0,110 0))

    Args:
        wkt (str): A well-known text geometry. This is assumed to use latitude and longitude
            coordinates in WGS84 (EPSG:4326).
        location (str): A location name used to select a query geometry from the SAFE gazetteer.
        distance (float): A search distance in metres. All geometries are converted to the
            UTM 50N projection to provide appopriate distance searching.
    """

    db = current.db

    # validate the query geometry options and report back if there is an error
    query_geom = dataset_parse_spatial(wkt, location)
    if isinstance(query_geom, dict):
        return query_geom

    qry = ((db.published_datasets.id == db.dataset_locations.dataset_id) &
           (db.dataset_locations.name == db.gazetteer.location) &
           ((db.gazetteer.wkt_utm50n.st_distance(query_geom) <= distance) |
            (db.dataset_locations.wkt_utm50n.st_distance(query_geom) <= distance)))

    return qry


def dataset_spatial_bbox_search(wkt=None, location=None, match_type='intersect', distance=None):
    """Spatial search for dataset bounding boxes

    The endpoint can search for datasets using either a user-provided geometry or a location
    name from the SAFE gazetteer. This endpoint uses only the dataset bounding box, which is
    provided for all datasets, rather than sampling location information which may not be
    recorded for some datasets.

    Examples:
        /api/search/bbox?wkt=Polygon((110 0, 110 10,120 10,120 0,110 0))
        /api/search/bbox?wkt=Polygon((116 4.5,116 5,117 5,117 4.5,116 4.5))
        /api/search/bbox?wkt=Polygon((116 4.5,116 5,117 5,117 4.5,116 4.5))&match_type=contain
        /api/search/bbox?wkt=Point(116.5 4.75)&match_type=within

    Args:
        wkt (str): A well-known text geometry. This is assumed to use latitude and longitude
            coordinates in WGS84 (EPSG:4326).
        location (str): A location name used to select a query geometry from the SAFE gazetteer.
        match_type (str): One of 'intersect', 'contain' and 'within' to match the provided geometry
            to the geographic extents of datasets. The 'contain' option returns datasets that
            completely cover the query geometry and 'within' returns datasets that fall entirely
            within the query geometry.

    """

    db = current.db

    # validate the query geometry options and report back if there is an error
    query_geom = dataset_parse_spatial(wkt, location)
    if isinstance(query_geom, dict):
        return query_geom

    # validate the match type
    if match_type not in ['intersect', 'contain', 'within', 'distance']:
        return {'error': 400, 'message': "Unknown spatial match type: {}".format(match_type)}

        # Query the geographic extents with the appropriate predicate
    if match_type == 'intersect':
        qry = (db.published_datasets.geographic_extent_utm50n.st_intersects(query_geom))
    elif match_type == 'contain':
        qry = (db.published_datasets.geographic_extent_utm50n.st_contains(query_geom))
    elif match_type == 'within':
        qry = (db.published_datasets.geographic_extent_utm50n.st_within(query_geom))
    elif match_type == 'distance':
        qry = (db.published_datasets.geographic_extent_utm50n.st_distance(query_geom) <= distance)

    return qry


def get_index():
    """
    Function to generate a JSON string containing the formatted contents of the dataset files
    table. This is used as the core index of the safedata package, so is cached in a file
    like format along with MD5 hashes of the index and other files to speed up checking for
    updates and refreshing the index.
    """

    # version of the data contained in the dataset description
    db = current.db
    qry = (db.published_datasets.id == db.dataset_files.dataset_id)
    val = dataset_query_to_json(qry,
                                fields=[('published_datasets', 'publication_date'),
                                        ('published_datasets', 'zenodo_concept_id'),
                                        ('published_datasets', 'zenodo_record_id'),
                                        ('published_datasets', 'dataset_access'),
                                        ('published_datasets', 'dataset_embargo'),
                                        ('published_datasets', 'dataset_title'),
                                        ('published_datasets', 'most_recent'),
                                        ('dataset_files', 'checksum'),
                                        ('dataset_files', 'filename'),
                                        ('dataset_files', 'filesize')])

    # repackage the db output into a single dictionary per file
    entries = val['entries'].as_list()
    [r['published_datasets'].update(r.pop('dataset_files')) for r in entries]
    val['entries'] = [r['published_datasets'] for r in entries]

    # Find the hashes
    index_hash = hashlib.md5(json(val).encode('utf-8')).hexdigest()

    # Use the file hash of the static gazetteer geojson
    gazetteer_file = os.path.join(current.request.folder, 'static', 'files', 'gis', 'gazetteer.geojson')
    with open(gazetteer_file) as f:
        gazetteer_hash = hashlib.md5(f.read().encode('utf-8')).hexdigest()

    # Use the file hash of the static locations alias csv
    location_aliases_file = os.path.join(current.request.folder, 'static', 'files', 'gis', 'location_aliases.csv')
    with open(location_aliases_file) as f:
        location_aliases_hash = hashlib.md5(f.read().encode('utf-8')).hexdigest()

    return dict(hashes=dict(index=index_hash,
                            gazetteer=gazetteer_hash,
                            location_aliases=location_aliases_hash),
                index=val)

    # return the dictionary - this will be json serialised by the API when it is returned
    return val

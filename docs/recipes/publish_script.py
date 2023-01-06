"""Example script to publish a dataset using sdv from within Python."""

import simplejson

from safedata_validator.field import Dataset
from safedata_validator.logger import LOGGER
from safedata_validator.resources import Resources
from safedata_validator.zenodo import post_metadata

resources = Resources("safedata_validator_local_test.cfg")
data_xlsx = "Test_format_good_NCBI.xlsx"

ds = Dataset(resources)
ds.load_from_workbook(data_xlsx)
ds.to_json()


# TODO - run the deposit creation and publication steps
zenodo_json = "zenodo_1143714.json"
data_json = "Test_format_good_NCBI.json"


with open(zenodo_json) as zn_json:
    znjson = simplejson.load(zn_json)

with open(data_json) as ds_json:
    dsjson = simplejson.load(ds_json)

# Run the post metadata function
response, error = post_metadata(zenodo=znjson, metadata=dsjson, resources=resources)

# Report on the outcome.
if error is not None:
    LOGGER.error(f"Failed to publish deposit: {error}")
else:
    LOGGER.info(f"Published to: {response['links']['record']}")

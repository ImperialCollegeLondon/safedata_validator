"""Example script to publish a dataset using sdv from within Python."""

import simplejson

from safedata_validator.field import Dataset
from safedata_validator.logger import LOGGER
from safedata_validator.resources import Resources
from safedata_validator.zenodo import post_metadata

# Local paths to the configuration file and the dataset to be validated
config_path = "safedata_validator_local_test.cfg"
data_xlsx = "Test_format_good_NCBI.xlsx"

# Create a Resources object from the config file and then create a
# dataset instance using those validation resources
resources = Resources(config_path)
ds = Dataset(resources)

# Load the dataset from the Excel workbook, which validates the content
ds.load_from_workbook(data_xlsx)

# Export the validate dataset metadata to JSON, creating the file
# Test_format_good_NCBI.json
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
    assert error is not None
    LOGGER.info(f"Published to: {response['links']['record']}")

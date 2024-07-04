"""Example script to validate a dataset using safedata_validator from within Python."""

import simplejson

from safedata_validator.field import Dataset
from safedata_validator.resources import Resources

# Local paths to the configuration file and the dataset to be validated
config_path = "config.cfg"
dataset = "SAFE_dataset.xlsx"

# Create a Resources object from the config file and then create a dataset instance
# using those validation resources
resources = Resources(config_path)
ds = Dataset(resources)

# Load the dataset from the Excel workbook, which validates the content
ds.load_from_workbook(dataset)

# Export the validated metadata to a file
with open("SAFE_dataset.json") as json_file:
    simplejson.dump(ds.to_json(), json_file)

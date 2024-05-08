"""Example script to publish a dataset using safedata_validator from within Python."""

import simplejson

from safedata_validator.field import Dataset
from safedata_validator.logger import LOGGER
from safedata_validator.resources import Resources
from safedata_validator.server import post_metadata
from safedata_validator.zenodo import (
    create_deposit,
    generate_inspire_xml,
    publish_deposit,
    upload_file,
    upload_metadata,
)

# Local paths to the configuration file and the dataset to be validated
config_path = "config.cfg"
dataset = "SAFE_dataset.xlsx"
extra_file = "Supplementary_files.zip"
xml_file = "SAFE_dataset_GEMINI.xml"

# Create a Resources object from the config file and then create a dataset instance
# using those validation resources
resources = Resources(config_path)
ds = Dataset(resources)

# Load the dataset from the Excel workbook, which validates the content
ds.load_from_workbook(dataset)

# Extract the validated dataset metadata
data_metadata = simplejson.loads(ds.to_json())

# Create the new deposit to publish the dataset
zenodo_metadata, error = create_deposit(resources=resources)

# Monitor the success of individual steps
all_good = error is None

# Generate XML
xml_content = generate_inspire_xml(
    dataset_metadata=data_metadata, zenodo_metadata=zenodo_metadata, resources=resources
)
with open(xml_file, "w") as xml_out:
    xml_out.write(xml_content)

# Post the files
for file in [dataset, extra_file, xml_file]:
    if all_good:
        file_upload_response, error = upload_file(
            metadata=zenodo_metadata, filepath=file, resources=resources
        )
        all_good = error is None

# Post the metadata
if all_good:
    md_upload_response, error = upload_metadata(
        metadata=data_metadata, zenodo=zenodo_metadata, resources=resources
    )
    all_good = error is None

# Publish the deposit
if all_good:
    publish_response, error = publish_deposit(
        zenodo=zenodo_metadata, resources=resources
    )
    all_good = error is None

# Show the new publication
publish_response["links"]["html"]

# Post the dataset metadata to the safedata server
if all_good:
    response, error = post_metadata(
        zenodo=zenodo_metadata, metadata=data_metadata, resources=resources
    )

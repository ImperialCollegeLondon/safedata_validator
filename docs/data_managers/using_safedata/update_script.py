"""Python script to publish a new version of a dataset using safedata_validator."""

from pathlib import Path

import simplejson

from safedata_validator.resources import Resources
from safedata_validator.server import post_metadata
from safedata_validator.zenodo import (
    create_deposit,
    delete_files,
    generate_inspire_xml,
    publish_deposit,
    upload_files,
    upload_metadata,
)

# Local paths to the files to be published
dataset = "Example.xlsx"
metadata_path = "Example.json"
extra_file = "Supplementary_files.zip"
xml_file = "Example_GEMINI.xml"

# Create a Resources object from a configuration file in a standard location
resources = Resources()

# Extract the validated dataset metadata
with open(metadata_path) as md_json:
    data_metadata = simplejson.load(md_json)

# Create a new version of an existing dataset using the record ID of the most recent
# version
zenodo_metadata, error = create_deposit(
    new_version=1143714,
    resources=resources,
)

# Monitor the success of individual steps
all_good = error is None

# Generate XML
xml_content = generate_inspire_xml(
    dataset_metadata=data_metadata, zenodo_metadata=zenodo_metadata, resources=resources
)
with open(xml_file, "w") as xml_out:
    xml_out.write(xml_content)

# Get the names of the existing files from the JSON metadata of the deposit and delete
# them all before uploading the provided versions. Note that the publish_dataset
# function does this in a much more sophisticated way.
existing_online_files = [p["key"] for p in zenodo_metadata["files"]]
if all_good:
    file_delete_response, error = delete_files(
        metadata=zenodo_metadata, filenames=existing_online_files, resources=resources
    )
    all_good = error is None


# Post the files
files = [Path(f) for f in (dataset, extra_file, xml_file)]
if all_good:
    file_upload_response, error = upload_files(
        metadata=zenodo_metadata, filepaths=files, resources=resources
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

# Post the dataset metadata to the safedata server
if all_good:
    response, error = post_metadata(
        zenodo=zenodo_metadata, metadata=data_metadata, resources=resources
    )

"""Python script to publish a dataset using safedata_validator.

Note that this is essentially just a slimmed down version of the publish_dataset
function underlying the `safedata_zenodo publish_dataset` command line tool.
"""

from pathlib import Path

import simplejson

from safedata_validator.resources import Resources
from safedata_validator.zenodo import (
    ZenodoResources,
    create_deposit,
    generate_inspire_xml,
    publish_deposit,
    upload_files,
    upload_metadata,
)

# Local paths to the dataset file
dataset = "Example.xlsx"
metadata_path = "Example.json"
extra_file = "Supplementary_files.zip"
xml_file = "Example_GEMINI.xml"

# Create a Resources object from a configuration file in a standard location and convert
# to the Zenodo specific resource class
resources = Resources()
zenodo_resources = ZenodoResources(resources)

# Extract the validated dataset metadata
with open(metadata_path) as md_json:
    data_metadata = simplejson.load(md_json)

# Create the new deposit to publish the dataset
create_response = create_deposit(zen_res=zenodo_resources)

# Bail if unsuccessful
if not create_response.ok:
    raise RuntimeError(create_response.error_message)

# Extract the Zenodo metadata from the response
zenodo_metadata = create_response.json_data

# Generate XML
xml_content = generate_inspire_xml(
    dataset_metadata=data_metadata, zenodo_metadata=zenodo_metadata, resources=resources
)
with open(xml_file, "w") as xml_out:
    xml_out.write(xml_content)

# Post the files
files = [Path(f) for f in (dataset, extra_file, xml_file)]
file_upload_response = upload_files(
    zenodo=zenodo_metadata, filepaths=files, zen_res=zenodo_resources
)

if not file_upload_response.ok:
    raise RuntimeError(file_upload_response.error_message)

# Post the metadata
md_upload_response = upload_metadata(
    metadata=data_metadata, zenodo=zenodo_metadata, zen_res=zenodo_resources
)

if not md_upload_response.ok:
    raise RuntimeError(md_upload_response.error_message)

# Publish the deposit
publish_response = publish_deposit(zenodo=zenodo_metadata, zen_res=zenodo_resources)

if not publish_response.ok:
    raise RuntimeError(publish_response.error_message)

# Show the new publication
print(publish_response.json_data["links"]["html"])

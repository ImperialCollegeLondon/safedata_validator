#!/bin/sh

# The code below assumes that the safedata_validator tools find a configuration
# file in one of the standard locations. If not, the path can be specified 
# explicitly using, the -r flag:
#   -r /path/to/safedata_validator_local_test.cfg 

# Validate the dataset. If successful this step will export a dataset metadata
# file called SAFE_dataset.json
safedata_validate SAFE_dataset.xlsx

# Publish the dataset to Zenodo
# 1) Create a new deposit as a new version of an existing record. Again this will
#    generate a deposit metadata file called something like zenodo_1156212.json
safedata_zenodo create_deposit -c 1143713

# 2) Generate a GEMINI XML metadata file for the deposit
safedata_zenodo generate_xml zenodo_1156212.json SAFE_dataset.json 1156212_GEMINI.xml

# 3) Upload the dataset file, external files named in the dataset summary and the XML
#    metadata. This uses the zenodo metadata file to confirm the upload destination.
safedata_zenodo upload_file zenodo_1156212.json SAFE_dataset.xlsx
safedata_zenodo upload_file zenodo_1156212.json Supplementary_files.zip
safedata_zenodo upload_file zenodo_1156212.json 1156212_GEMINI.xml

# 3) Update the Zenodo deposit webpage - this populates the deposit description
#    on Zenodo from the dataset metadata
safedata_zenodo upload_metadata zenodo_1156212.json SAFE_dataset.json

# 4) Finally, publish the deposit to create the final record and DOI
safedata_zenodo publish_deposit zenodo_1156212.json

# Publish the dataset and zenodo metadata to the metadata server
safedata_metadata post_metadata zenodo_1156212.json SAFE_dataset.json
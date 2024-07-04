#!/bin/sh

# The code below assumes that:
# * the safedata_validator tools find a configuration file in one of the 
#   standard locations. If not, the path can be specified explicitly 
#   using, the -r flag:
#     -r /path/to/safedata_validator_local_test.cfg 
# * the Example.xlsx file has been sucessfully validated, generating the
#   Example.json metadata file

# Publish the dataset to Zenodo
# 1) Create a new deposit, which will generate a deposit metadata file called
#    something like zenodo_1143714.json
safedata_zenodo create_deposit

# 2) Generate a GEMINI XML metadata file for the deposit
safedata_zenodo generate_xml zenodo_1143714.json Example.json 1143714_GEMINI.xml

# 3) Upload the dataset file, external files named in the dataset summary and the XML
#    metadata. This uses the zenodo metadata file to confirm the upload destination.
safedata_zenodo upload_file zenodo_1143714.json Example.xlsx
safedata_zenodo upload_file zenodo_1143714.json Supplementary_files.zip
safedata_zenodo upload_file zenodo_1143714.json 1143714_GEMINI.xml

# 4) Update the Zenodo deposit webpage - this populates the deposit description
#    on Zenodo from the dataset metadata
safedata_zenodo upload_metadata zenodo_1143714.json Example.json

# 5) Finally, publish the deposit to create the final record and DOI
safedata_zenodo publish_deposit zenodo_1143714.json
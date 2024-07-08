#!/bin/sh

# The code below assumes that:
# * the safedata_validator tools find a configuration file in one of the 
#   standard locations. If not, the path can be specified explicitly 
#   using, the -r flag:
#     -r /path/to/safedata_validator_local_test.cfg 
# * the Example.xlsx file has been sucessfully validated, generating the
#   Example.json metadata file


# Publish the dataset to Zenodo
# 1) Create a new deposit as a new version of the most recent version of an existing
#    record. Again this will generate a deposit metadata file called something like
#    zenodo_1156212.json
safedata_zenodo create_deposit --new-version 1143714

# 2) Generate a GEMINI XML metadata file for the deposit
safedata_zenodo generate_xml zenodo_1156212.json Example.json 1156212_GEMINI.xml

# 3) Delete the existing files on the deposit. This uses the zenodo metadata file to
#    confirm the upload destination. 
#
#    Note that in this example, we do not check to see if any of the files are
#    identical. The command below simple deletes all of the files and then reuploads the
#    provided versions. The publish_dataset subcommand handles this in a much more
#    sophisticated way.
safedata_zenodo delete_files zenodo_1156212.json \
    Example.xlsx  Supplementary_files.zip 1143714_GEMINI.xml

# 4) Upload the dataset file, external files named in the dataset summary and the new
#    XML metadata. This uses the zenodo metadata file to confirm the upload destination.
safedata_zenodo upload_files zenodo_1156212.json \
    Example.xlsx  Supplementary_files.zip 1156212_GEMINI.xml

# 3) Update the Zenodo deposit webpage - this populates the deposit description
#    on Zenodo from the dataset metadata
safedata_zenodo upload_metadata zenodo_1156212.json Example.json

# 4) Finally, publish the deposit to create the final record and DOI
safedata_zenodo publish_deposit zenodo_1156212.json

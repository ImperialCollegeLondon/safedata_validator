#!/bin/sh

# The code below assumes that the safedata_validator tools find a configuration
# file in one of the standard locations. If not, the path can be specified 
# explicitly using, the -r flag:
#   -r /path/to/safedata_validator_local_test.cfg 

# Validate the dataset. If successful this step will export a dataset metadata
# file called SAFE_dataset.json
safedata_validate SAFE_dataset.xlsx
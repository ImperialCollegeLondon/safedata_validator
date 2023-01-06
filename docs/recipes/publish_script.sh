#!/bin/sh

safedata_validate -r safedata_validator_local_test.cfg Test_format_good_NCBI.xlsx


safedata_zenodo -r safedata_validator_local_test.cfg  create_deposit
safedata_zenodo -r safedata_validator_local_test.cfg  upload_file zenodo_1143714.json Test_format_good_NCBI.xlsx
safedata_zenodo -r safedata_validator_local_test.cfg  upload_metadata zenodo_1143714.json Test_format_good_NCBI.json
safedata_zenodo -r safedata_validator_local_test.cfg  publish_deposit zenodo_1143714.json


safedata_server -r safedata_validator_local_test.cfg post_metadata zenodo_1143714.json Test_format_good_NCBI.json
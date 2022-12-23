# This script echoes the command line usage for script tools to files in this
# directory that can then be transcluded into markdown files using the python
# markdown extension markdown-include. The syntax is:
#
# {!path/to/file!}
#
# Note that the path is from the root directory running the mkdocs build/serve
# process, not relative to the markdown file in which the included file is 
# referenced.


echo "cl_prompt $ safedata_validate -h" > safedata_validate.txt
safedata_validate -h >> safedata_validate.txt


echo "cl_prompt $ safedata_zenodo -h" > safedata_zenodo_top.txt
safedata_zenodo -h >> safedata_zenodo_top.txt

subcommands=(create_deposit get_deposit discard_deposit publish_deposit \
             upload_file delete_file upload_metadata amend_metadata \
             sync_local_dir maintain_ris generate_html generate_xml show_config)

for subc in "${subcommands[@]}";
do
    echo $subc
    echo "cl_prompt $ safedata_zenodo $subc -h" > safedata_zenodo_$subc.txt
    safedata_zenodo $subc -h >> safedata_zenodo_$subc.txt
done

echo "cl_prompt $ safedata_server -h" > safedata_server_top.txt
safedata_server -h >> safedata_server_top.txt

subcommands=(post_metadata update_gazetteer)

for subc in "${subcommands[@]}";
do
    echo $subc
    echo "cl_prompt $ safedata_server $subc -h" > safedata_server_$subc.txt
    safedata_server $subc -h >> safedata_server_$subc.txt
done

echo "cl_prompt $ safedata_build_local_gbif -h" > safedata_build_local_gbif.txt
safedata_build_local_gbif -h >> safedata_build_local_gbif.txt

echo "cl_prompt $ safedata_build_local_ncbi -h" > safedata_build_local_ncbi.txt
safedata_build_local_ncbi -h >> safedata_build_local_ncbi.txt

# This script echoes the command line usage for script tools to files in this
# directory that can then be transcluded into markdown files using the 
# mkdocs plugin include-markdown. The syntax is:
#
# {%
# include "path/to/file.txt"
# %}
#
# Note that the path can be relative to the docs directory or relative to the markdown
# file in which the included file is referenced.


echo "cl_prompt $ safedata_validate -h" > safedata_validate.txt
safedata_validate -h >> safedata_validate.txt


echo "cl_prompt $ safedata_zenodo -h" > safedata_zenodo_top.txt
safedata_zenodo -h >> safedata_zenodo_top.txt

subcommands=(create_deposit get_deposit discard_deposit publish_deposit \
             upload_file delete_file upload_metadata amend_metadata \
             sync_local_dir maintain_ris generate_html generate_xml)

for subc in "${subcommands[@]}";
do
    echo $subc
    echo "cl_prompt $ safedata_zenodo $subc -h" > safedata_zenodo_$subc.txt
    safedata_zenodo $subc -h >> safedata_zenodo_$subc.txt
done

echo "cl_prompt $ safedata_metadata -h" > safedata_metadata_top.txt
safedata_metadata -h >> safedata_metadata_top.txt

subcommands=(post_metadata update_resources)

for subc in "${subcommands[@]}";
do
    echo $subc
    echo "cl_prompt $ safedata_metadata $subc -h" > safedata_metadata_$subc.txt
    safedata_metadata $subc -h >> safedata_metadata_$subc.txt
done

echo "cl_prompt $ safedata_build_local_gbif -h" > safedata_build_local_gbif.txt
safedata_build_local_gbif -h >> safedata_build_local_gbif.txt

echo "cl_prompt $ safedata_build_local_ncbi -h" > safedata_build_local_ncbi.txt
safedata_build_local_ncbi -h >> safedata_build_local_ncbi.txt

# Validating a dataset

<!-- markdownlint-disable MD046 MD033 -->

The description below assumes that a user has provided a SAFE formatted dataset and an
additional ZIP file containing additional files:

* `Example.xlsx`
* `Supplementary_files.zip`

The validation process is then typically to use the
[`safedata_validate`](../command_line_tools/safedata_validate.md) command line tool to
run the validation, as shown in the tab below. You can also use the `safedata_validator`
package from within Python to do this, but you would only typically do this if you were
building this functionality into another application.

=== "bash"

    ```sh
    {%
    include "data_managers/using_safedata/validate_script.sh"
    %}
    ```

=== "python"

    ```python
    {%
    include "data_managers/using_safedata/validate_script.py"
    %}
    ```

Running validation checks the contents of the data file against the
[required data format](../../data_providers/data_format/overview.md) and generates a
report on the validation process. This report details the checks being carried out,
highlights any issues with the dataset and generates a final pass or fail verdict.
The notes in the report for any warning and errors raised should help resolve the
problems and datasets must not be published until they pass validation!

The report details from a successful validation can be seen below:

```text
{%
include "data_providers/data_format/Example.log"
%}
```

Successful validation will also generate a JSON file containing the detailed metadata
for the dataset. This file is used to provide information for the publication process
and also the detailed metadata to be uploaded to the metadata server. This file is very
long and is not really intended to be read, but the contents can be shown below.

<details><summary>JSON metadata contents</summary>

```json
{%
include "data_providers/data_format/Example.json"
%}
```

</details>

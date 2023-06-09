"""Test an entry point.

This script is used to debug the entry point function imported and run below. This
function expecting to be provided with values via command line arguments and argparse,
and so debugging requires a way to pass those arguments in to the function.

Within Visual Studio Code, a new debug configuration can be added to .vscode/launch.json
containing the required arguments, which will then be passed into the script as it runs.
The JSON to do that looks like this:

```json
    {
        "name": "_safedata_validator_cli debug",
        "type": "python",
        "request": "launch",
        "program": "${file}",
        "console": "integratedTerminal",
        "args": [
            "-r",
            "local/config_files/live_safeweb_gbif_2016_07_25.cfg",
            "../safedata_directory/1198522/1237724/template_SinghRamesh.xlsx"
        ]
    }
```

From the RUN AND DEBUG pane:

* make sure you have _this_ file as the active file to be debugged,
* select the configuration name above from the dropdown,
* start the debugger **from the RUN AND DEBUG window**.

This will then run the entry point function below and provide it with the required
arguments given in the 'args' data above. Breakpoints in the code can then be used to
interact with the running code and allow step through.
"""

from safedata_validator.entry_points import _safedata_validator_cli

_safedata_validator_cli()

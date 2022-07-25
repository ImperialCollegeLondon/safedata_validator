"""Validate a safedata directory.

This script runs the validation on a set of safedata formatted files stored in the
directory structure used in the `safedata` R package. This revalidation is extremely
useful for debugging.

It would be good to track the memory use in loading datasets, to understand how chunk
size affects the memory profile and run time. However @profile only records the change
in memory usage from running each step, not the memory use while running it, and I can't
get the time based profile API to run.
"""

import csv
import os
from datetime import datetime
from glob import glob

from safedata_validator.field import Dataset
from safedata_validator.logger import CONSOLE_HANDLER, COUNTER_HANDLER, LOG
from safedata_validator.resources import Resources

safe_dir = "/Users/dorme/Research/SAFE/Database/safedata_directory"
resource_file = (
    "/Users/dorme/Research/SAFE/Database/safedata_validator_package/local_db/"
    "safedata_validator.cfg"
)

# Get the xlsx files in the directory - these are always nested 2 directories deep:
# concept_id/deposit_id/filename.xlsx

test_files = glob("*/*/*.xlsx", recursive=True, root_dir=safe_dir)
resources = Resources(resource_file)

COUNTER_HANDLER.setLevel("ERROR")
CONSOLE_HANDLER.setLevel(100)  # Muted

files = []
errors = []


# @profile
def process_file(idx: int, this_file: str, resources: Resources):
    """Validate a single file."""

    # get file name parts
    file_parts = this_file.split(os.path.sep)

    # Record the start of processing
    time_in = datetime.now()

    # Create the dataset from this file
    ds = Dataset(resources=resources)

    try:
        full_path = os.path.join(safe_dir, this_file)
        ds.load_from_workbook(full_path, console_log=False)

        # # This won't run - should get memory usage while loading is going on, but
        # # throws an error
        # mem = memory_usage((ds.load_from_workbook,
        #                     tuple(),
        #                     {'filename': this_file,
        #                         'console_log': False}), interval=0.05, timeout=1)

        # Record completion time
        time_out = datetime.now()

        # Record the errors
        if ds.n_errors:
            # split the log entries on the error sign, dropping the blank entry before
            # the first error sign.
            error_log = LOG.getvalue().split("!")[1:]
            for each_err in error_log:
                errors.append((this_file, each_err))

        # Report and store processing time
        filesize = os.path.getsize(full_path)
        proc_time = (time_out - time_in).total_seconds()

        print(
            f"{idx:3} {(filesize/ 10**6):7.2f} MB {proc_time:7.3f} "
            f"{ds.n_errors:3} {this_file} "
        )

        files.append([this_file, filesize, proc_time, ds.n_errors, ds.error_breakdown])

        # output json
        jsonf = f"{file_parts[0]}_{file_parts[1]}.json"
        with open(jsonf, "w") as jsonout:
            jsonout.write(ds.to_json())

    except Exception as e:
        print(f"{idx:3} {this_file} failed {str(e)}")


# Loop over the files
for idx, this_file in enumerate(test_files):

    process_file(idx, this_file, resources)

# Save a file by file summary of processing time and file size
with open("testing_gbif20xx.csv", "w") as csvfile:
    csvwriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)

    csvwriter.writerow(
        (
            "filepath",
            "filesize",
            "proc_time",
            "n_errors",
            "summary_errors",
            "loc_errors",
            "gbif_errors",
            "ncbi_errors",
            "data_errors",
        )
    )
    for rw in files:

        # Pop off the error breakdown and flatten
        erb = rw.pop()
        erb = list(erb.values())
        erb = erb[0:5] + [sum([c for n, c in erb[5]])]
        rw_unpack = rw[0:3] + erb

        csvwriter.writerow(rw_unpack)

# Save a table of individual errors by file
with open("errors_gbif20xx.csv", "w") as csvfile:
    csvwriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)

    csvwriter.writerow(("filepath", "error_msg"))

    for fname, each_err in errors:
        csvwriter.writerow((fname, each_err.strip()))

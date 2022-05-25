"""
This script runs the validation on a set of safedata formatted files stored in
the directory structure used in the `safedata` R package. This revalidation is
extremely useful for debugging.
"""

import os
from glob import glob
import humanize
import csv
from datetime import datetime

# HACK - Problem here with circular imports - need to load locations first
#        but more broadly need to fix the structure!
from safedata_validator.locations import Locations
from safedata_validator.field import Dataset
from safedata_validator.logger import COUNTER_HANDLER, CONSOLE_HANDLER, LOG
from safedata_validator.resources import Resources

test_files = glob('/Users/dorme/Research/SAFE/Database/safedata_directory/*/*/*.xlsx', recursive=True)
resources = Resources('/Users/dorme/Research/SAFE/Database/safedata_validator_package/local_db/safedata_validator.cfg')

COUNTER_HANDLER.setLevel('ERROR')
CONSOLE_HANDLER.setLevel(100) # Muted

files = []
errors = {}

# Loop over the files
for idx, this_file in enumerate(test_files):
    
    # Record the start of processing
    time_in = datetime.now()

    # Create the dataset from this file
    ds = Dataset(resources=resources)
    ds.load_from_workbook(this_file, console_log=False)

    # Record completion time
    time_out = datetime.now()

    # Check the errors
    if ds.n_errors:
        errors[this_file] = LOG.getvalue()

    # Report and store processing time
    filesize = os.path.getsize(this_file)
    proc_time = (time_out - time_in).total_seconds()

    print(f'{idx} {humanize.naturalsize(filesize)} {proc_time} {this_file}')

    files.append((this_file, filesize, proc_time))



with open('testing.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for rw in files:
        spamwriter.writerow(rw)


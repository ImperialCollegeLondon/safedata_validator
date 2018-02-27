import os
from safe_dataset_checker import *


loc_json = 'SAFE_locations.json'
gbif_db = 'backbone-current-simple.sqlite'


# part by part
ds = Dataset('Test_format_good.xlsx', gbif_database=gbif_db)
ds.load_summary()
ds.load_locations(locations_json=loc_json)
ds.load_taxa()
ds.load_data_worksheet(ds.dataworksheet_summaries[0])
ds.load_data_worksheet(ds.dataworksheet_summaries[1])
ds.final_checks()


ds = Dataset('Test_format_bad.xlsx', gbif_database=gbif_db)
ds.load_summary()
ds.load_locations(locations_json=loc_json)
ds.load_taxa()
ds.load_data_worksheet(ds.dataworksheet_summaries[0])
ds.load_data_worksheet(ds.dataworksheet_summaries[1])
ds.final_checks()

# Checking logging modes

check_file('Test_format_good.xlsx', locations_json=loc_json, gbif_database=gbif_db)
sys.stderr.flush()
sys.stdout.flush()
sys.stdout.write('#\n# Checking file remotely and silently\n#\n')
check_file('Test_format_good.xlsx', locations_json=loc_json, verbose=False)

check_file('Test_format_bad.xlsx', locations_json=loc_json, gbif_database=gbif_db)
sys.stderr.flush()
sys.stdout.flush()
sys.stdout.write('#\n# Checking file remotely and silently\n#\n')
check_file('Test_format_bad.xlsx', locations_json=loc_json, verbose=False)


# Checking against a batch of files

datasets = os.walk('datasets')
results = {}

for pth, subdirs, files in datasets:
    if '/' in pth:
        xlfiles = [os.path.join(pth, fl) for fl in files if fl.endswith('.xlsx')]
        for xlf in xlfiles:
            results[xlf] = check_file(xlf, locations_json=loc_json, gbif_database=gbif_db, verbose=False)



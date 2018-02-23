import sys
from safe_dataset_checker import *


loc_json = 'SAFE_locations.json'
gbif_db = 'backbone-current-simple.sqlite'


check_file('Test_format_good.xlsx', locations_json=loc_json, gbif_database=gbif_db)
sys.stderr.flush()
sys.stdout.write('#\n# Checking file remotely and silently\n#\n')
check_file('Test_format_good.xlsx', locations_json=loc_json, verbose=False)

check_file('Test_format_bad.xlsx', locations_json=loc_json, gbif_database=gbif_db)
sys.stderr.flush()
sys.stdout.write('#\n# Checking file remotely and silently\n#\n')
check_file('Test_format_bad.xlsx', locations_json=loc_json, verbose=False)

import os
from safe_dataset_checker import *

# ./safe_dataset_checker.py -g backbone-current-simple.sqlite -l SAFE_locations.json Test_format_good.xlsx

loc_json = 'SAFE_locations.json'
gbif_db = 'backbone-current-simple.sqlite'
fname = 'Test_format_good.xlsx'
fname = '/Users/dorme/Downloads/Pillay_R_et_al_Dryobalanops_lanceolata_AllData.xlsx'

# part by part
ds = Dataset(fname, gbif_database=gbif_db)
ds.load_summary()
ds.load_locations(locations_json=loc_json)
ds.load_taxa()

for meta in ds.dataworksheet_summaries:
    ds.load_data_worksheet(meta)

ds.final_checks()


ds = Dataset('Test_format_bad.xlsx', gbif_database=gbif_db)
ds.load_summary()
ds.load_locations(locations_json=loc_json)
ds.load_taxa()
ds.load_data_worksheet(ds.dataworksheet_summaries[0])
ds.load_data_worksheet(ds.dataworksheet_summaries[1])
ds.load_data_worksheet(ds.dataworksheet_summaries[2])
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
xlfiles = []
for pth, subdirs, files in datasets:
    if '/' in pth:
        xlfiles.extend([os.path.join(pth, fl) for fl in files
                        if fl.endswith('.xlsx') and not fl.startswith('~$')])

for xlf in xlfiles:
    print xlf
    res = check_file(xlf, locations_json=loc_json, gbif_database=gbif_db, verbose=False)
    results[xlf] = (res.passed, res.report().getvalue())


for k, v in results.iteritems():
    if not v[0]:
        print k
        print v[1]



print(results['datasets/6/datasets.file.ba4584070320ed17.74656d706c6174655f47726179416e74436f6d7065746974696f6e2e786c7378.xlsx'].passed)
print(results['datasets/6/datasets.file.ba4584070320ed17.74656d706c6174655f47726179416e74436f6d7065746974696f6e2e786c7378.xlsx'].report().getvalue())
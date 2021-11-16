# A basic script for me to test how taxon works and to test that my new functions also work

# Import necessary files
from safedata_validator.dataset import Dataset

# File name of excel file
fname = "example_inputs/Test_format_good.xlsx"

# initialise the dataset object
dataset = Dataset(fname, verbose=True, gbif_database=None, locations=None)

# Command to run the taxa checks
dataset.load_taxa()

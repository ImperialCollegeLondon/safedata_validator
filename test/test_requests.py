import requests

"""This test is to double check that the fake filesystem correctly handles
using the requests package.
"""


URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=science%5bjournal%5d+AND+breast+cancer+AND+2008%5bpdat%5d"

def test_requests():

    x = requests.get(URL)

def test_requests_fake(config_filesystem):

    x = requests.get(URL)

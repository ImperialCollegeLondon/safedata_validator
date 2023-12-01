"""This test checks that the fake filesystem correctly handles the requests package."""

import requests  # type: ignore

URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&"
    "term=science%5bjournal%5d+AND+breast+cancer+AND+2008%5bpdat%5d"
)


def test_requests():

    requests.get(URL)


def test_requests_fake(config_filesystem):

    requests.get(URL)

"""This test checks that the fake filesystem correctly handles the requests package."""

import requests  # type: ignore

URL = "https://zenodo.org/communities/safe/"


def test_requests():
    requests.get(URL)


def test_requests_fake(config_filesystem):
    requests.get(URL)

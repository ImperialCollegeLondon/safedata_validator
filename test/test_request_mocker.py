"""Prototype requests mocking for tests."""

import requests


def test_mocking(mocked_requests):
    """Prototype test.

    Note that this will fail if the SDV_DO_NOT_MOCK_REQUESTS environment variable is
    set.
    """
    assert requests.get("http://google.com").text == "data"

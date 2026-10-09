"""Test the zenodo module functionality.

This is currently extremely limited - needs a mocked API to Zenodo to do actual testing
of most of the functionality. At the moment, the only thing that is tested is that the
function used to merge Zenodo metadata into dataset metadata for publication to a
safedata_server instance works as intended.
"""

import json
from contextlib import contextmanager

from .conftest import FIXTURE_FILES


@contextmanager
def does_not_raise():
    yield


def test_merge_metadata(user_config_file):
    """Test the merge_metadata function."""

    from safedata_validator.models import UploadMetadata
    from safedata_validator.zenodo import merge_metadata

    # Need to create fake files to be used for the validation outputs
    user_config_file.create_file("/tmp/publishable_metadata.json")

    dataset_metadata = json.load(open(FIXTURE_FILES.rf.good_seq_file_dataset_json))
    zenodo_metadata = json.load(open(FIXTURE_FILES.rf.good_seq_file_zenodo_json))

    _ = merge_metadata(
        path="/tmp/publishable_metadata.json",
        dataset_metadata=dataset_metadata,
        zenodo_metadata=zenodo_metadata,
    )

    # Check the output JSON data can be loaded and validated by the pydantic model
    UploadMetadata.model_validate(json.load(open("/tmp/publishable_metadata.json")))

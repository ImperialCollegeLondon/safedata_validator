"""This module provides functions to interact with the metadata server web application.

1. send new dataset metadata to the server
2. update the resources on the server to match the local versions.
"""  # D415

from __future__ import annotations

import json
from collections import namedtuple
from dataclasses import dataclass, field
from hashlib import md5
from pathlib import Path

import requests  # type: ignore
from pydantic import ValidationError

from safedata_validator.logger import LOGGER
from safedata_validator.models import UploadMetadata
from safedata_validator.resources import Resources


@dataclass
class MetadataResources:
    """Packaging for Metadata resources.

    This dataclass is used to package the Metadata server specific elements of the
    configuration.
    """

    resources: Resources
    """A safedata_validator resources instance."""
    api: str = field(init=False)
    """The configured Zenodo API to be used."""
    token: dict[str, str] = field(init=False)
    """A dictionary providing the authentication token for the API."""

    def __post_init__(self) -> None:
        """Populate the post init attributes."""

        # Get the appropriate API and token
        self.api = self.resources.metadata.api
        self.headers = {"Authorization": f"Token {self.resources.metadata.token}"}
        self.ssl_verify = self.resources.metadata.ssl_verify


def post_metadata(metadata_file: Path, server_resources: MetadataResources) -> bool:
    """Post the dataset metadata to the metadata server.

    Args:
        metadata_file: The path to the published metadata for a dataset
        server_resources: The server resources to be used.

    Returns:
        False on failure or True on success
    """

    try:
        with open(metadata_file) as fobj:
            _ = UploadMetadata.model_validate(json.load(fobj))
    except IsADirectoryError:
        LOGGER.error(f"Metadata file is a directory: {metadata_file}")
        return False
    except FileNotFoundError:
        LOGGER.error(f"Metadata file not found: {metadata_file}")
        return False
    except (UnicodeDecodeError, json.JSONDecodeError):
        LOGGER.error(f"Could not parse JSON file: {metadata_file}")
        return False
    except ValidationError as excep:
        LOGGER.error(f"Validation found errors in: {metadata_file}")
        LOGGER.error(str(excep))
        return False

    # post the metadata to the server
    response = requests.post(
        f"{server_resources.api}/api/datasets/upload/",
        headers=server_resources.headers,
        files={"file": open(metadata_file, "rb")},
        verify=server_resources.ssl_verify,
    )

    # Check what is in the response is received from the server
    response_data = response.json()
    if not response.ok:
        LOGGER.error(f"Failed to post metadata: {response_data['detail']}")
        return False

    LOGGER.info(f"Metadata posted: {server_resources.api}{response_data['url']}")
    return True


ResourceData = namedtuple(
    "ResourceData", ["label", "md5_key", "upload_api", "local_path"]
)
"""Named tuple for resource data."""


def update_resources(server_resources: MetadataResources) -> bool:
    """Update the resources on the metadata server.

    The metadata server provides the gazetteer, location aliases and any project IDs as
    part of the safedata R package workflow. The web server also uses those resources
    internally to provide information. This function is used to post the current
    resources to an API on the server that is used to refresh those reseources.

    Args:
        server_resources: The server resources to be used.

    Returns:
        False on failure or True on success
    """

    # Setup endpoints to upload files
    endpoints: list[ResourceData] = [
        ResourceData(
            label="Gazetteer",
            md5_key="gazetteer",
            upload_api="api/gazetteer/upload/",
            local_path=server_resources.resources.gaz_path,
        ),
        ResourceData(
            label="Location aliases",
            md5_key="location_aliases",
            upload_api="api/gazetteer/aliases/upload/",
            local_path=server_resources.resources.localias_path,
        ),
    ]

    # Add the project database if provided
    if server_resources.resources.project_database is not None:
        LOGGER.warning("Metadata server does not currently handle project database")

        # endpoints.append(
        #     ResourceData(
        #         label="Project Database",
        #         md5_key="projects",
        #         upload_api="projects",
        #         local_path=server_resources.resources.project_database,
        #     )
        # )

    success = True

    # Get the md5 digest of the remote file
    remote_md5_dict = requests.get(
        f"{server_resources.api}/resources/hashes.json"
    ).json()

    # Post changed files
    # - Note that these files have been validated by the Resources  class
    #   so no further validation here.
    for resource_data in endpoints:
        # Get the md5 digest for the local file
        local_md5 = md5(open(resource_data.local_path, "rb").read()).hexdigest()

        # Update only if the digests differ
        if local_md5 == remote_md5_dict[resource_data.md5_key]:
            LOGGER.info(f"{resource_data.label} up to date")
        else:
            LOGGER.info(f"{resource_data.label} update starting...")
            # post the resource files to the server
            response = requests.post(
                f"{server_resources.api}/{resource_data.upload_api}",
                headers=server_resources.headers,
                files={"file": open(resource_data.local_path, "rb")},
            )

            if not response.ok:
                LOGGER.error(f"{resource_data.label} update failed:")
                LOGGER.error(response.json()["detail"])
                success = False
            else:
                LOGGER.info(f"{resource_data.label} updated")

    if not success:
        LOGGER.error("Failed to update all resources.")
    else:
        LOGGER.info("Resources updated")

    return success

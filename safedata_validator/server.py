"""This module provides functions to interact with the metadata server web application.

1. send new dataset metadata to the server
2. update the resources on the server to match the local versions.
"""  # D415

from __future__ import annotations

import json
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

    # Setup endpoints to upload files to
    endpoints: list[tuple[str, str, str]] = [
        ("Gazetteer", "gazetteer", server_resources.resources.gaz_path),
        (
            "Location aliases",
            "gazetteer/aliases",
            server_resources.resources.localias_path,
        ),
    ]

    # Add the project database if provided
    if server_resources.resources.project_database is not None:
        endpoints.append(
            (
                "Project database",
                "projects",
                server_resources.resources.project_database,
            )
        )

    success = True

    # Post changed files - note that these files have been validated by the Resources
    # class so no further validation here.
    for name, endpoint, file in endpoints:
        # Get the md5 digest for the local file
        local_md5 = md5(open(file, "rb").read()).hexdigest()
        # Get the md5 digest of the remote file
        get_remote = requests.get(f"{server_resources.api}/api/{endpoint}/hash/")
        remote_md5 = get_remote.json()["md5"]

        # Update only if the digests differ
        if local_md5 == remote_md5:
            LOGGER.info(f"{name} up to date")
        else:
            LOGGER.info(f"{name} update starting...")
            # post the resource files to the server
            response = requests.post(
                f"{server_resources.api}/api/{endpoint}/upload/",
                headers=server_resources.headers,
                files={"file": open(file, "rb")},
            )

            if not response.ok:
                LOGGER.error(f"{name} update failed:")
                LOGGER.error(response.json()["detail"])
                success = False
            else:
                LOGGER.info(f"{name} updated")

    if not success:
        LOGGER.error("Failed to update all resources.")
    else:
        LOGGER.info("Resources updated")

    return success

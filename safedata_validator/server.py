"""This module provides functions to interact with the metadata server web application.

1. send new dataset metadata to the server
2. update the resources on the server to match the local versions.
"""  # D415

from __future__ import annotations

from dataclasses import InitVar, dataclass, field

import requests  # type: ignore

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
        self.token = {"access_token": self.resources.metadata.token}
        self.ssl_verify = self.resources.metadata.ssl_verify


@dataclass
class MetadataResponse:
    """Metadata server response processor.

    This dataclass is a processor around `requests.Response` objects from calls to a
    metadata server. If the response is successful, it parses the returned data payload;
    otherwise it formats as much information as possible into an error message.
    """

    response: InitVar[requests.Response]
    """The incoming response from a Zenodo API call."""
    ok: bool = field(init=False)
    """Was the response ok."""
    status_code: int = field(init=False)
    """The status code returned by the response."""
    json_data: dict = field(init=False, default_factory=lambda: dict())
    """The JSON data payload from a successful response."""
    error_message: str | None = field(init=False, default=None)
    """A formatted error message from a failed response."""

    def __post_init__(self, response: requests.Response) -> None:
        """Populate the ZenodoResponse object."""
        # Basic status
        self.ok = response.ok
        self.status_code = response.status_code
        # Now either populate json data or the error message
        if self.ok:
            self.json_data = response.json()
        else:
            self.error_message = response.text


def post_metadata(
    metadata: dict, zenodo: dict, server_resources: MetadataResources
) -> MetadataResponse:
    """Post the dataset metadata and zenodo metadata to the metadata server.

    Args:
        metadata: The dataset metadata dictionary for a dataset
        zenodo: The dataset metadata dictionary for a deposit
        server_resources: The server resources to be used.

    Returns:
        See [here][safedata_validator.server.MetadataResources].
    """

    # Get payload
    payload = {"metadata": metadata, "zenodo": zenodo}

    # post the metadata to the server
    return MetadataResponse(
        requests.post(
            f"{server_resources.api}/post_metadata",
            params=server_resources.token,
            json=payload,
            verify=server_resources.ssl_verify,
        )
    )


def update_resources(server_resources: MetadataResources) -> MetadataResponse:
    """Update the resources on the metadata server.

    The metadata server provides the gazetteer, location aliases and any project IDs as
    part of the safedata R package workflow. The web server also uses those resources
    internally to provide information. This function is used to post the current
    resources to an API on the server that is used to refresh those reseources.

    Args:
        server_resources: The server resources to be used.

    Returns:
        See [here][safedata_validator.server.MetadataResources].
    """

    # Get payload
    files = {
        "gazetteer": open(server_resources.resources.gaz_path, "rb"),
        "location_aliases": open(server_resources.resources.localias_path, "rb"),
    }

    if server_resources.resources.project_database is not None:
        files["project_database"] = open(
            server_resources.resources.project_database, "rb"
        )

    # post the resource files to the server
    return MetadataResponse(
        requests.post(
            f"{server_resources.api}/update_resources",
            params=server_resources.token,
            files=files,
        )
    )

"""This module provides functions to interact with the metadata server web application.

1. send new dataset metadata to the server
2. update the resources on the server to match the local versions.
"""  # D415

from __future__ import annotations

import requests  # type: ignore

from safedata_validator.resources import Resources
from safedata_validator.zenodo import (
    ZenodoFunctionResponseType,
    _resources_to_zenodo_api,
)


def post_metadata(
    metadata: dict, zenodo: dict, resources: Resources | None = None
) -> ZenodoFunctionResponseType:
    """Post the dataset metadata and zenodo metadata to the metadata server.

    Args:
        metadata: The dataset metadata dictionary for a dataset
        zenodo: The dataset metadata dictionary for a deposit
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)

    payload = {"metadata": metadata, "zenodo": zenodo}

    # post the metadata to the server
    mtd = requests.post(
        f"{zres['mdapi']}/post_metadata",
        params={"token": zres["mdtoken"]},
        json=payload,
        verify=zres["mdssl"],
    )

    # trap errors in uploading metadata and tidy up
    if mtd.status_code != 201:
        return {}, mtd.text
    else:
        return mtd.json(), None


def update_resources(resources: Resources) -> ZenodoFunctionResponseType:
    """Update the resources on the metadata server.

    The metadata server provides the gazetteer, location aliases and any project IDs as
    part of the safedata R package workflow. The web server also uses those resources
    internally to provide information. This function is used to post the current
    resources to an API on the server that is used to refresh those reseources.

    Args:
        resources: The safedata_validator resource configuration to be used. If
            none is provided, the standard locations are checked.

    Returns:
        See [here][safedata_validator.zenodo.ZenodoFunctionResponseType].
    """

    # Get resource configuration
    zres = _resources_to_zenodo_api(resources)

    # Get payload
    files = {
        "gazetteer": open(resources.gaz_path, "rb"),
        "location_aliases": open(resources.localias_path, "rb"),
    }

    if resources.project_database is not None:
        files["project_database"] = open(resources.project_database, "rb")

    # post the resource files to the server
    response = requests.post(
        f"{zres['mdapi']}/update_resources",
        params={"token": zres["mdtoken"]},
        files=files,
    )

    # Trap errors in uploading resources and tidy up
    if response.status_code != 201:
        return {}, response.text
    else:
        return {}, None

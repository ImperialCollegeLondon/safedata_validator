"""This module provides pydantic models for validating the different JSON files used in
``safedata_validator``. There are three main models:

* The :class:`Summary` model validates the metadata stored in the Summary sheet.
* The :class:`Zenodo` model validates the metadata provided from Zenodo for published
  datasets.
* The :class:`UploadMetadata` model validates the combination of the two - a ``Summary``
  model that also includes the published Zenodo metadata.

These models currently do not attempt to do the fine scale validation carried out in the
{mod}`safedata_validator.summary` module. It is possible that the functionality in that
module could be replaced by extending the validation on these models.
"""  # noqa: D205, D415

from __future__ import annotations

from datetime import date
from typing import Any, Literal

from pydantic import BaseModel


class Author(BaseModel):
    """Metadata summary author block."""

    name: str
    affiliation: str | None = None  # Affiliation not mandatory
    email: str | None = None  # Email not mandatory
    orcid: str | None = None  # ORCID not mandatory


class ExternalFile(BaseModel):
    """Metadata summary external file block."""

    file: str
    description: str


class Field(BaseModel):
    """Metadata data worksheet field header metadata."""

    field_name: str
    description: str
    field_type: str
    units: str | None
    method: str | None
    levels: str | None
    taxon_field: str | None
    taxon_name: str | None
    interaction_field: str | None
    interaction_name: str | None
    range: str | None
    col_idx: int


class Dataworksheet(BaseModel):
    """Metadata data worksheet summary metadata."""

    taxa_fields: list
    max_row: int
    max_col: int
    name: str
    title: str
    description: str
    descriptors: list[str]
    external: str | None = None  # Not mandatory
    fields: list[Field]
    field_name_row: int
    n_data_row: int


class Funders(BaseModel):
    """Metadata summary funders block."""

    body: str
    type: str
    ref: str | int | None = None
    url: str | None = None


class Permits(BaseModel):
    """Metadata summary permits block."""

    type: str
    authority: str
    number: str


class Locations(BaseModel):
    """Metadata locations."""

    name: str
    new_location: bool
    wkt_wgs84: str | None = None


class GBIFTaxa(BaseModel):
    """Metadata GBIF Taxa."""

    worksheet_name: str | None  # None for parent taxa
    taxon_id: int
    parent_id: int | None  # None for root taxon.
    taxon_name: str
    taxon_rank: str
    taxon_status: str


class Summary(BaseModel):
    """Metadata summary model."""

    project_ids: list[int]
    title: str
    description: str
    authors: list[Author]
    filename: str
    external_files: list[ExternalFile] | None = None
    access: Literal["Open", "Embargo", "Restricted", "open", "embargo", "restricted"]
    embargo_date: date | None = None
    access_conditions: str | None = None
    funders: list[Funders] | None = None
    permits: list[Permits] | None = None
    keywords: list[str]
    dataworksheets: list[Dataworksheet]
    gbif_timestamp: str
    gbif_taxa: list[GBIFTaxa] | None = None
    sequenced_taxa: dict[str, Any]  # Needs work
    locations: list[Locations] | None = None
    validator_version: str
    temporal_extent: list[str]
    latitudinal_extent: list[float]
    longitudinal_extent: list[float]


class Creator(BaseModel):
    """Zenodo creator metadata."""

    name: str
    affiliation: str | None = None
    orcid: str | None = None


class Contributor(BaseModel):
    """Zenodo contributor metadata."""

    name: str
    affiliation: str
    type: str
    orcid: str


class Community(BaseModel):
    """Zenodo community metadata."""

    identifier: str


class PrereserveDoi(BaseModel):
    """Zenodo prereserve DOI metadata."""

    doi: str
    recid: int


class ZenodoMetadata(BaseModel):
    """Zenodo metadata section."""

    title: str
    doi: str
    publication_date: str
    description: str
    access_right: str
    creators: list[Creator]
    contributors: list[Contributor] | None = None
    keywords: list[str] | None = None
    license: str | None = None  # Missing on restricted records
    imprint_publisher: str
    communities: list[Community]
    upload_type: str
    prereserve_doi: PrereserveDoi


class Links(BaseModel):
    """Zenodo deposit links."""

    self: str
    html: str
    doi: str
    parent_doi: str | None = None  # Not sure why missing in some zenodo metadata
    badge: str
    conceptbadge: str | None = None  # Not sure why missing in some zenodo metadata
    files: str
    bucket: str
    latest_draft: str
    latest_draft_html: str
    publish: str
    edit: str
    discard: str
    newversion: str
    record: str
    record_html: str
    latest: str
    latest_html: str


class FileLinks(BaseModel):
    """Zenodo file links metadata."""

    self: str
    download: str


class File(BaseModel):
    """Zenodo file metadata."""

    id: str
    filename: str
    filesize: int
    checksum: str
    links: FileLinks


class Zenodo(BaseModel):
    """Zenodo metadata model."""

    created: str
    modified: str
    id: int
    conceptrecid: str
    doi: str
    conceptdoi: str
    doi_url: str
    metadata: ZenodoMetadata
    title: str
    links: Links
    record_id: int
    owner: int
    files: list[File]
    state: str
    submitted: bool


class UploadMetadata(Summary):
    """Complete upload metadata containing summary and published metadata."""

    zenodo: Zenodo

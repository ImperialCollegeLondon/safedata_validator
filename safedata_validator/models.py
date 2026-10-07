"""This module provides pydantic models for validating the different JSON files used in
``safedata_validator``. There are three main models:

* The :class:`Summary` model validates the metadata stored in the Summary sheet.
* The :class:`Zenodo` model validates the metadata provided from Zenodo for published
  datasets.
* The :class:`UploadMetadata` model validates the combination of the two - a ``Summary``
  model that also includes the published Zenodo metadata.

These models currently do not attempt to do the fine scale validation carried out in the
{mod}`safedata_validator.summary` module. It is possible that the functionality in that
module could be replaced by extending the validation on these models and that could then
become the basis of a `summary.toml` file for use with a CSV based deposit system.
"""  # noqa: D205, D415

# NOTES
#  - current metadata fails on setting an EmailStr validator (badly formed emails)

from __future__ import annotations

from datetime import date
from typing import Any, Literal

from pydantic import AnyHttpUrl, BaseModel, ConfigDict, model_validator


class ForbiddenExtra(BaseModel):
    """Base class to enforce core configuration settings."""

    model_config = ConfigDict(extra="forbid")


class Author(ForbiddenExtra):
    """Metadata summary author block."""

    name: str
    affiliation: str | None = None  # Affiliation not mandatory
    email: str | None  # Email not mandatory
    # Should be this, but identifies failures
    # email: EmailStr | None = None  # Email not mandatory
    orcid: str | None = None  # ORCID not mandatory


class ExternalFile(ForbiddenExtra):
    """Metadata summary external file block."""

    file: str
    description: str


class Field(ForbiddenExtra):
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


class Dataworksheet(ForbiddenExtra):
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


class Funders(ForbiddenExtra):
    """Metadata summary funders block."""

    body: str
    type: str
    ref: str | int | None = None
    url: AnyHttpUrl | None = None


class Permits(ForbiddenExtra):
    """Metadata summary permits block."""

    type: str
    authority: str
    number: str


class Locations(ForbiddenExtra):
    """Metadata locations."""

    name: str
    new_location: bool
    wkt_wgs84: str | None = None


class GBIFTaxa(ForbiddenExtra):
    """Metadata GBIF Taxa."""

    worksheet_name: str | None  # None for parent taxa
    taxon_id: int
    parent_id: int | None  # None for root taxon.
    taxon_name: str
    taxon_rank: str
    taxon_status: str


class Summary(ForbiddenExtra):
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
    temporal_extent: list[date]
    latitudinal_extent: list[float]
    longitudinal_extent: list[float]


class Creator(ForbiddenExtra):
    """Zenodo creator metadata."""

    name: str
    affiliation: str | None = None
    orcid: str | None = None


class Contributor(ForbiddenExtra):
    """Zenodo contributor metadata."""

    name: str
    affiliation: str
    type: str
    orcid: str


class Community(ForbiddenExtra):
    """Zenodo community metadata."""

    identifier: str


class PrereserveDoi(ForbiddenExtra):
    """Zenodo prereserve DOI metadata."""

    doi: str
    recid: int


class RelatedIdentifier(ForbiddenExtra):
    """Zenodo related identifier metadata."""

    identifier: str
    relation: str
    resource_type: str
    scheme: str


class Grants(ForbiddenExtra):
    """Zenodo grants metadata."""

    id: str


class ZenodoMetadata(ForbiddenExtra):
    """Zenodo metadata section."""

    title: str
    doi: str
    publication_date: date
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
    # These last three are only found in a handful of existing records.
    related_identifiers: list[RelatedIdentifier] | None = None
    grants: list[Grants] | None = None
    version: str | None = None


class Links(ForbiddenExtra):
    """Zenodo deposit links."""

    self: AnyHttpUrl
    html: AnyHttpUrl
    doi: AnyHttpUrl
    # Not sure why missing in some zenodo metadata
    parent_doi: AnyHttpUrl | None = None
    badge: AnyHttpUrl
    # Not sure why missing in some zenodo metadata
    conceptbadge: AnyHttpUrl | None = None
    files: AnyHttpUrl
    bucket: AnyHttpUrl
    thumb250: AnyHttpUrl | None = None
    thumbs: dict[str, AnyHttpUrl] | None = None
    latest_draft: AnyHttpUrl
    latest_draft_html: AnyHttpUrl
    publish: AnyHttpUrl
    edit: AnyHttpUrl
    discard: AnyHttpUrl
    newversion: AnyHttpUrl
    record: AnyHttpUrl
    record_html: AnyHttpUrl
    latest: AnyHttpUrl
    latest_html: AnyHttpUrl


class FileLinks(ForbiddenExtra):
    """Zenodo file links metadata."""

    self: AnyHttpUrl
    download: AnyHttpUrl


class File(ForbiddenExtra):
    """Zenodo file metadata."""

    id: str
    filename: str
    filesize: int
    checksum: str
    links: FileLinks


class Zenodo(ForbiddenExtra):
    """Zenodo metadata model."""

    created: str
    modified: str
    id: int
    conceptrecid: str
    doi: str
    conceptdoi: str
    doi_url: AnyHttpUrl
    metadata: ZenodoMetadata
    title: str
    links: Links
    record_id: int
    owner: int
    files: list[File]
    state: str
    submitted: bool

    @model_validator(mode="after")
    def check_published(self) -> Zenodo:
        """Check the record is published."""

        if (self.state != "done") or not self.submitted:
            raise ValueError(
                "Zenodo metadata do not show completed publication of the dataset"
            )

        return self


# Enforce
# - submitted must be true and state 'done'


class UploadMetadata(Summary):
    """Complete upload metadata containing summary and published metadata."""

    zenodo: Zenodo

# The `taxa` submodule

::: safedata_validator.taxa
    options:
        show_source: true
        group_by_category: false
        members: false
        show_root_heading: false
        show_root_toc_entry: false

## GBIF validation

::: safedata_validator.taxa.GBIFTaxon
    options:
        show_root_heading: True
        header_level: 3

::: safedata_validator.taxa.GBIFValidator
    options:
        show_root_heading: True
        header_level: 3
        members: true
        filters: ["!__del__"]

::: safedata_validator.taxa.GBIFError
    options:
        show_root_heading: True
        header_level: 3

## NCBI Validation

::: safedata_validator.taxa.NCBITaxon
    options:
        show_root_heading: True
        header_level: 3

::: safedata_validator.taxa.NCBIValidator
    options:
        show_root_heading: True
        header_level: 3
        members: true
        filters: ["!__del__"]

::: safedata_validator.taxa.NCBIError
    options:
        show_root_heading: True
        header_level: 3

## Taxon worksheet classes

::: safedata_validator.taxa.GBIFTaxa
    options:
        show_root_heading: True
        header_level: 3

::: safedata_validator.taxa.NCBITaxa
    options:
        show_root_heading: True
        header_level: 3

::: safedata_validator.taxa.Taxa
    options:
        show_root_heading: True
        header_level: 3

## Helper functions

::: safedata_validator.taxa.taxon_index_to_text
    options:
        show_root_heading: True
        header_level: 3

::: safedata_validator.taxa.taxa_strip
    options:
        show_root_heading: True
        header_level: 3

::: safedata_validator.taxa.construct_bi_or_tri
    options:
        show_root_heading: True
        header_level: 3

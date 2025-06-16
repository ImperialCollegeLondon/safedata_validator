# Sequenced taxa worksheets

<!-- markdownlint-disable MD033 -->
<style>

/*fixing cell widths so everything lines up and adding borders*/
table {
  table-layout: fixed;
}

tbody td {
  width: 14em;
  min-width: 14em;
  max-width: 14em;
  border: 1px solid lightgrey;
}

thead th {
  width: 14em;
  min-width: 14em;
  max-width: 14em;
  border: 1px solid lightgrey;
}
</style>
<!-- markdownlint-enable MD033 -->

Like the [GBIFTaxa worksheet](./gbif_taxa.md), a sequenced taxa worksheet is used to
record taxonomic details of organisms referred to in the Data worksheets, but is
intended for use with taxa identified through sequencing rather than field observation.
Like the GBIFTaxa worksheet, the taxonomic information is used to:

* Cross-check the listed taxa against taxon fields in the [data worksheets](./data.md)
  to confirm that there is a complete list of the taxa recorded in the datasets.
* Generate an hierarchical taxon index that can be used to search datasets for taxa of
  interest.

However, unlike the GBIF taxa worksheets, `safedata_validator` does not validate the
provided taxa against the underlying bioinformatic datasets. Latin binomial names from
field observations can be matched against any broad taxonomic database - we use GBIF
because it is well-curated and updated - but taxonomic ids from sequencing can only
really be validated against the specific bioinformatic database and version used by the
researchers. This is beyond the scope of `safedata_validator`, so no validation is
performed. We instead expect researchers to provide simple information in the Summary
data about the bioinformatic dataset used to generate sequenced taxa id, which could
be used to revisit taxon identification.

A sequenced taxa worksheet is used to provide a table of taxonomic ranks for matched
sequences. This is a common export format from bioinformatics workflows and allows
`safedata_validator` to provide an hierarchical taxon index for these taxa as well as
providing a list of sequenced taxon names to be checked against the taxon names provided
in [data worksheets](./data.md).

Multiple sequenced taxa worksheets can be provided. The names of these sheets and the
details (name, version, etc) of the reference databases used to generate them should be
provided in the Summary metadata. If taxa are used anywhere in the dataset at least one
worksheet of this type or the GBIFTaxa worksheet must be included. It is also an option
to provide both a GBIFTaxa worksheet and one or more sequenced taxa worksheets, where
taxa identified through both sequencing and field observation are being reported.

We would also encourage you (where possible) to include the raw sequencing data used to
generate your taxonomies, as it improves replicability to have the raw genomic
information in addition to your taxonomic assignment of sequences presented in the
sequenced taxa worksheets. Sequencing data should be included as [external
files](other_formats.md), and you should note in the comments which sequence worksheet
it provides raw data for.

## Taxon table layout

The table format looks like this:

<!-- markdownlint-disable MD013 -->
{{ read_excel('Example.xlsx', sheet_name = 'Sequenced', keep_default_na = False, tablefmt = 'github') }}
<!-- markdownlint-enable MD013 -->

The table must contain column headers in the **first row** of the worksheet. The Name
column is mandatory and must contain a local name for **all** of the taxa that you are
going to use in the rest of the dataset, aside those that are described in another
sequenced taxa worksheet or in a GBIFTaxa worksheet. The other columns are to help us
maintain a taxonomic index for the taxa used in your datasets.

* **Name**: This column is mandatory and these are the names that you are going to use
  in your data worksheet.  These can be abbreviations or codes: if you want to use
  `E_coli_FaMsAh` in your data worksheets, rather than typing out `Escherichia coli
  FaMsAh gene for 16S rRNA, LC842012` every time, then that is fine.

  Please note that different worksheet names need to be used if a taxon is detected in
  multiple ways and hence represented in [multiple taxa
  worksheets](./taxa.md#duplicate-taxa).

* **Ranks**: Here the column name (e.g. Phylum) provides a **taxonomic rank**, and the
  row entries provide the relevant names for this rank. A top level rank must be
  provided, this is one of `domain`, `superkingdom` or `kingdom`. If either `domain` or
  `superkingdom` is provided `kingdom` can also be provided. The other thing ranks that
  can be provided are the standard backbone ranks (`phylum` down to `species`). It is up
  to you which of these you provide but if a rank is provided every higher backbone rank
  must be provided, e.g. if `order` is provided as a rank `class` and `phylum` must also
  be provided. Non-backbone ranks can be provided (e.g. `subspecies` or `strain`) but
  they are treated as additional information and are therefore not validated.

    !!! Note

        Missing rank entries are generally completely fine, e.g. if you have genera that
        haven't been assigned to a family level no entry has to be provided for the
        rank. However, an entry for highest taxonomic rank has to be provided for every
        row.

    Names can be provided in plain text, or alternatively in a commonly used notation,
    where the rank is indicated by a lower case first letter and the name follows after
    two underscores (e.g. `k__Bacteria` for Kingdom Bacteria). Notation of this type
    should be placed in the correct rank columns, and validation is carried out to check
    that the rank implied by the notation matches the column rank.

    Entries for species rank should not be provided as binomials! Instead, the species
    name and the genus name should be provided separately, and the validator will
    construct the binomial automatically for the searchable metadata. Standard taxonomic
    tags (e.g. "candidatus") are fine to include as part of the name, however they are
    removed from the searchable metadata.
  
* **Comments and other fields**: These fields are obviously optional. If you do have
  particular notes that you want to make - new species notes and the like - then these
  can be very useful for future researchers trying to place taxa. Equally if you want to
  record further information about the taxon rows, you can add additional fields as long
  as they do not duplicate any of the field names mentioned above.

## My data is not sequencing data

You should record this data using GBIF format on a GBIFTaxa worksheet instead.

## My data doesn't contain taxa

Fine. You can omit both the GBIFTaxa worksheet and the sequenced taxa worksheets!

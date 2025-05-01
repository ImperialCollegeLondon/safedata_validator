# The SeqTaxa worksheet

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

This worksheet plays a similar role to the [GBIFTaxa worksheet](./gbif_taxa.md), that is
recording the taxonomic information for organisms referred to in the Data worksheets.
The key difference is that this taxonomic information is not validated against a
reference database. The reason for this is that software to generate taxonomies from
sequencing data use a wide range of different reference databases, and it is not
appropriate for a piece of validation software like `safedata_validator` to dictate
which reference databases can be used. This approach is only suitable for sequencing
data, and this worksheet is recommended for this class of data. If taxa are used
anywhere in the dataset either this worksheet or the GBIFTaxa worksheet must be
included. It is also an option to provide both a GBIFTaxa worksheet and a SeqTaxa
worksheet, e.g. in cases where both sequencing and observational data are being
reported.

TODO - Once multiple sheets are implemented I need to describe this + talk about the
required sheet metadata

## Taxon table layout

The table format looks like this:

TODO - Once I've settled I've added the metadata this should be populated from the
template

The table must contain column headers in the **first row** of the worksheet. The Name
column is mandatory and must contain a local name for **all** of the taxa that you are
going to use in the rest of the dataset, aside those that are already described on a
GBIFTaxa worksheet.

If both a SeqTaxa and a GBIFTaxa worksheet are provided the same
taxa can be included in both, e.g. a species found both by observation and eDNA
sequencing. However, to avoid confusion these should be given different names, i.e.
`Vulpes_obs` and `Vulpes_seq` for observed and sequenced instances of `Vulpes`,
respectively. Names cannot be duplicated either within a SeqTaxa worksheet or from a
GBIFTaxa worksheet (when one exists)! Note that these can be abbreviations or codes:
if you want to use `Crbe` in your data worksheets, rather than typing out
`Crematogaster borneensis` every time, then that is fine.

!!! Note

    These are the names that you are going to use in your data worksheet. The
    other columns are to help us validate the taxonomy of your names.

* **Name**: This column can be optionally used to note that a row contains a new taxon
  that is not expected to be present in the NCBI database. The taxon will be included as
  a new taxon as a child of the next taxonomic rank.

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
        rank. However, an entry for at least one of the top level ranks has to be
        provided for every row.

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
  
* **Comments and other fields**: These fields are obviously optional. If you
  do have particular notes that you want to make - explaining disagreements with NCBI
  taxonomy, new species notes and the like - then these can be very useful for future
  researchers trying to place taxa. Equally if you want to record further information
  about NCBI taxon rows, you can add additional fields as long as they do not duplicate
  any of the field names mentioned above.

## My data is not sequencing data

You should record this data using GBIF format on a GBIFTaxa worksheet instead.

## My data doesn't contain taxa

Fine. You can omit either or both of the GBIFTaxa and SeqTaxa worksheets!

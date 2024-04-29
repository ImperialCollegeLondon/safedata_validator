# The NCBITaxa worksheet

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
The key difference is that taxonomic information is recorded using NCBI taxonomy rather
than GBIF taxonomy. This taxonomy better matches the output of sequencing data, and so
this worksheet is recommended for this class of data. If taxa are used anywhere in the
dataset either this worksheet or the GBIFTaxa worksheet must be included. It is also an
option to provide both a GBIFTaxa worksheet and a NCBITaxa worksheet, e.g. in cases
where both sequencing and observational data are being reported.

## NCBI Taxon validation

In order to help keep the taxonomy as clean as possible and to allow us to index the
taxonomic coverage of datasets, we will check all taxon names in NCBITaxa worksheet
against the NCBI taxonomy database. If you want to check your taxon names and ranks,
then the search engine is here:

[https://www.ncbi.nlm.nih.gov/taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)

No online taxonomy is ever going to be 100% up to date (or 100% agree with your
taxonomic usage!) but the NCBI backbone has extremely good taxonomic coverage.

## Taxon table layout

The table format looks like this:

<!-- markdownlint-disable MD013 -->
{{ read_excel('Example.xlsx', sheet_name = 'NCBITaxa', keep_default_na = False, tablefmt = 'github') }}
<!-- markdownlint-enable MD013 -->

The table must contain column headers in the **first row** of the worksheet. The Name
column is mandatory and must contain a local name for **all** of the taxa that you are
going to use in the rest of the dataset, aside those that are already described on a
GBIFTaxa worksheet.

If both a NCBITaxa and a GBIFTaxa worksheet are provided the same
taxa can be included in both, e.g. a species found both by observation and eDNA
sequencing. However, to avoid confusion these should be given different names, i.e.
`Vulpes_obs` and `Vulpes_seq` for observed and sequenced instances of `Vulpes`,
respectively. Names cannot be duplicated either within a NCBITaxa worksheet or from a
GBIFTaxa worksheet (when one exists)! Note that these can be abbreviations or codes:
if you want to use `Crbe` in your data worksheets, rather than typing out
`Crematogaster borneensis` every time, then that is fine.

!!! Note

    These are the names that you are going to use in your data worksheet. The
    other columns are to help us validate the taxonomy of your names.

* **New**: This column can be optionally used to note that a row contains a new taxon
  that is not expected to be present in the NCBI database. The taxon will be included as
  a new taxon as a child of the next taxonomic rank.

* **Ranks**: Here the column name (e.g. Phylum) provides a **taxonomic rank**, and the
  row entries provide the relevant names for this rank. In contrast to GBIF, which only
  uses a small set of backbone ranks, the NCBI database also includes a large number of
  intermediate ranks (e.g. subphylum, strain). Any of these ranks may be included as
  headers in the worksheet, with the exception of `clade` and `no rank` as these ranks
  can be duplicated within a taxon hierarchy.
  
    You are only required to provide the taxonomic name for the specific rank that you
  are trying to match and can leave other fields empty. However it is probably more
  useful to provide a more complete taxonomy! If you do provide higher taxonomic
  information then it must be congruent with the hierarchy for the focal taxon. For
  example, specifying Family Anatidae is sufficient to identify waterfowl, but providing
  Order Carnivora (rather than Order Anseriformes) would result in an error.

    Names can be provided in plain text, or alternatively in a commonly used notation,
  where the rank is indicated by a lower case first letter and the name follows after
  two underscores (e.g. `k__Bacteria` for Kingdom Bacteria). Notation of this type
  should be placed in the correct rank columns, and validation is carried out to check
  that the rank implied by the notation matches the column rank.

    Two special cases are that NCBI outputs typically separate out the components of
  binomial and trinomial names: for example, they might return `g__Escherichia` and
  `s__coli`. In order to be able to match _complete_ species and subspecies names
  against the database, you must provide field information for genus, species and
  subspecies ranks. This information is used to assemble complete names for validation
  against the NCBI names. Note that if you have already compiled complete names, so that
  your genus field contains `Escherichia` and your species field contains `Escherichia
  coli` then this will also be accepted, as long as the parts are compatible.
  
    !!! Note

        Missing rank entries are completely fine, e.g. leaving out phylum information
        for some taxa but providing it for others. However, sufficient information must
        be provided to unambiguously identify each taxon.

* **Comments and other fields**: These fields are obviously optional. If you
  do have particular notes that you want to make - explaining disagreements with NCBI
  taxonomy, new species notes and the like - then these can be very useful for future
  researchers trying to place taxa. Equally if you want to record further information
  about NCBI taxon rows, you can add additional fields as long as they do not duplicate
  any of the field names mentioned above.

## Common issues

### Kingdom or Superkingdom

 NCBI defines a taxonomic rank above kingdom, which it terms a superkingdom. Bacteria is
 defined as a superkingdom as is Eukaryota. Within the NCBI database Bacterial taxa have
 no kingdom defined, but Eukaryotes generally have kingdom information provided. As
 such, we allow Bacteria (and Archaea) to be entered as **either** a superkingdom **or**
 a kingdom, but Eukaryota can **only** be entered as a superkingdom. It is, however, a
 perfectly valid option to not enter superkingdom information at all and just enter
 kingdom information for Eukaryotes (e.g. Fungi, Metazoa, etc) and phylum information
 for Prokaryotes (Bacteria and Archaea).

### Non-canon NCBI taxonomy

 The NCBI taxonomy database is regularly updated (particularly for microbial taxa). This
 means that taxa names can become synonymised or superseded. Generally, the current
 canonical name for taxa can be found based on superseded information and this will
 generate a warning like this:

    ? Non-canon name Enterococcus coli at rank species: synonym for Escherichia coli

 Then, both the superseded taxa information and the up to date information will be added
 to the taxon index. If you don't want both to be recorded, simply replace the
 superseded taxa information with the up to date information.

### Non-backbone ranks

 We consider backbone ranks to be those of the GBIF backbone (e.g. kingdom, phylum, ...,
 subspecies) with the addition of superkingdom. It is also fine to include non-backbone
 ranks such as strain or superorder. However, when the lineage of each taxa is found,
 only backbone ranks will be included, i.e. non-backbone ranks will only be recorded if
 they are the least nested taxonomic level for a specific taxon.

## My data is not sequencing data, and is hard to convert to NCBI taxonomy

You should record this data using GBIF format on a Taxa worksheet instead.

## My data doesn't contain taxa

Fine. You can omit either or both of the GBIFTaxa and NCBITaxa worksheets!

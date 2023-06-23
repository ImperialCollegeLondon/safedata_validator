# The NCBITaxa worksheet

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

| Name | NCBI ID | Superkingdom | Kingdom | Phylum | Class | Comments |
|---|---|---|---|---|---| --- |
| G_proteobacteria | 1236 | Bacteria |  | Proteobacteria | Gammaproteobacteria |  |
| E_mycetes | 147545 | Eukaryota | Fungi | Ascomycota | Eurotiomycetes |  |
| Dinophyceae |  | Eukaryota |  |  | Dinophyceae |  |
| Acidobact |  |  | k__Bacteria | p__Acidobacteria | c__Acidobacteriia |  |  |

The table must contain column headers in the **first row** of the worksheet. The first
two columns (Name, NCBI ID) are mandatory and contain the following:

* **Name**: This column must contain a local name for **all** of the taxa that you are
  going to use in the rest of the dataset, aside those that are already described on a
  GBIFTaxa worksheet. If both a NCBITaxa and a GBIFTaxa worksheet are provided the same
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

* **NCBI ID**: This column is required but providing actual entries is optional, and can
  be filled in with the ID that NCBI assigns to the taxon. This helps to ensure that the
  correct match is found in the NCBI database.

* **Ranks**: Here the column name (e.g. Phylum) provides a **taxonomic rank**, and the
  row entries provide the relevant names for this rank. At least two columns using GBIF
  backbone ranks must be provided. Non-backbone ranks (e.g. subphylum, strain) can also
  be provided, so long as they are defined in the NCBI database. The ranks should be
  listed in descending order from left to right. It is fine to skip specific backbone
  ranks entirely. However, if species is provided as a rank, genus must also be
  provided, and the same goes for subspecies and species. Species information can either
  be provided as a binomial or as a single name. Likewise subspecies information can
  either be provided as the full trinomial or as a single name. Names can be provided in
  plain text, or alternatively in a commonly used notation, where the rank is indicated
  by a lower case first letter and the name follows after two underscores (e.g.
  k__Bacteria for Kingdom Bacteria). Notation of this type should be placed in the
  correct rank columns, and validation is carried out to check that the rank implied by
  the notation matches the column rank.

!!! Note

    Missing rank entries are completely fine, e.g. leaving out
    phylum information for some taxa but providing it for others.
    However, sufficient information must be provided to unambiguously
    identify each taxon.

* **Comments**: This is entirely optional - if you have a fairly standard set of taxa
  with no serious issues then you can leave it out entirely or it can be empty. If you
  do have particular notes that you want to make - explaining disagreements with NCBI
  taxonomy, new species notes and the like - then these can be very useful for future
  researchers trying to place taxa.

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

### Superseded NCBI taxonomy

 The NCBI taxonomy database is regularly updated (particularly for microbial taxa). This
 means that both taxa names and NCBI IDs can become superseded. Generally, the up to
 date entries can be found based on superseded information, when this is the case a
 warning will be issued in this format:

    ? Bacillus foraminis not accepted usage should be Mesobacillus foraminis instead

 Then, both the superseded taxa information and the up to date information will be
 recorded. If you don't want both to be recorded, simply replace the superseded taxa
 information with the up to date information.

### Non-backbone ranks

 We consider backbone ranks to be those of the GBIF backbone (e.g. kingdom, phylum, ...,
 subspecies) with the addition of superkingdom. It is also fine to include non-backbone
 ranks such as strain or superorder. However, when the lineage of each taxa is found,
 only backbone ranks will be included, i.e. non-backbone ranks will only be recorded if
 they are the lowest taxonomic level for a specific taxon. Clades should only be
 included if they are the lowest known rank for a specific taxa. If they are used they
 should all be included as a single column, as repeated column names results in an
 error. Only the ordering of backbone ranks are checked so clades can be placed out of
 order to accomplish this.

### New and unrecognised taxa

 If a taxon is new or not recognised by NCBI (and you're sure it isn't just a typo!) it
 can still be recorded. This requires that the entry for the next lowest rank in the
 taxon's hierarchy is recognised by NCBI. If so the worksheet name, rank and parent taxa
 details of the original taxa are recorded. When this happens a warning of the following
 type will be generated:

    ? My new strain not registered with NCBI, but higher level taxon Escherichia coli is

## My data is not sequencing data, and is hard to convert to NCBI taxonomy

You should record this data using GBIF format on a Taxa worksheet instead.

## My data doesn't contain taxa

Fine. You can omit both the GBIFTaxa and NCBITaxa worksheet!

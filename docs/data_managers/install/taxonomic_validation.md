# Taxonomic validation process

This document provides an overview of the validation process used for GBIF and NCBI
taxa.

## GBIF validation

GBIF validation searches the GBIF backbone taxonomy database for a taxon that matches
both the taxon scientific name and taxonomic rank provided in a dataset. If a match is
found then the GBIF database will supply one of the following status codes:

```txt
accepted
doubtful 
homotypic synonym
synonym
heterotypic synonym
proparte synonym
misapplied
```

The first two options (accepted and doubtful) are treated as canonical matches, but the
remaining status codes will provide a link to an **accepted** taxonomic usage. The
validation process therefore extracts three things from the backbone database:

1. The status of the name and rank provided by the user.
2. The accepted usage (which might be the same).
3. The backbone taxonomic hierarchy for the taxa, so we can index at higher taxonomic
  levels.

However, this is complicated as there is no guarantee that taxa will always hook into
the backbone at the next deeper taxonomic level. For example, the accepted taxon
_Wanosuchus atresus_ is only hooked in at order level: it is accepted as a species, but
its parent is Crocodylia as the genus is doubtful. Similarly, _Goniopholis tenuidens_ is
a synonym at species level but again has Crocodylia as a parent (and is considered a
synonym for the family Goniopholidae).

In some cases, taxa can hook in at taxonomic levels _more nested_ than their own: the
species _Brittonastrum greenei_ is a synonym of the **subspecies** _Agastache
pallidiflora pallidiflora_.

The `parent_taxon_level_analysis.R` file in this repository contains some code
to check this:

* All **accepted taxa** map to a more nested parent but 5% map to a more nested parent
  more than one step up the hierarchy. The table below shows child taxon level as rows
  and parent taxon level as columns.

<!-- markdownlint-disable MD013 -->
|           | kingdom| phylum| class| order| family|   genus| species| subspecies| variety| form|
|:----------|-------:|------:|-----:|-----:|------:|-------:|-------:|----------:|-------:|----:|
|kingdom    |       0|      0|     0|     0|      0|       0|       0|          0|       0|    0|
|phylum     |     100|      0|     0|     0|      0|       0|       0|          0|       0|    0|
|class      |       5|    316|     0|     0|      0|       0|       0|          0|       0|    0|
|order      |       7|     45|  1327|     0|      0|       0|       0|          0|       0|    0|
|family     |    2191|   1339|  4267| 14423|      0|       0|       0|          0|       0|    0|
|genus      |    3427|   4985|  5584|  6260| 220735|       0|       0|          0|       0|    0|
|species    |    1567|    706|  1529|   696|   8944| 2449414|       0|          0|       0|    0|
|subspecies |      41|      7|     3|     2|    832|     268|  200902|          0|       0|    0|
|variety    |      53|     10|     0|    26|   2661|      50|   82914|         32|       0|    0|
|form       |      12|      4|     0|     4|    815|      18|   19272|          0|      56|    0|
<!-- markdownlint-enable MD013 -->

* Only 77% of **unaccepted taxa** map to a parent at the next most nested
  taxonomic level and 4.5% map to a parent at the same or a less nested level,
  as in the example above.

<!-- markdownlint-disable MD013 -->
|           | kingdom| phylum| class| order| family|   genus| species| subspecies| variety| form|
|:----------|-------:|------:|-----:|-----:|------:|-------:|-------:|----------:|-------:|----:|
|kingdom    |       0|      0|     0|     0|      0|       0|       0|          0|       0|    0|
|phylum     |      22|      8|     0|     0|      0|       0|       0|          0|       0|    0|
|class      |       0|     14|     1|     0|      0|       0|       0|          0|       0|    0|
|order      |       0|      5|    32|     0|      0|       0|       0|          0|       0|    0|
|family     |      21|    157|   481|  3599|      0|       0|       0|          0|       0|    0|
|genus      |    8555|  24242| 25055| 31010| 185911|       0|       0|          0|       0|    0|
|species    |      64|     24|   173|   405|   2142| 1886329|  121225|         84|       5|    0|
|subspecies |       3|      0|     1|     0|    151|   77512|   26266|         13|       0|    0|
|variety    |       2|      1|     0|     2|    367|  212954|   50062|         47|       4|    0|
|form       |       0|      0|     0|     0|    128|   48126|   10449|          3|       2|    0|
<!-- markdownlint-enable MD013 -->

### GBIF validation process

Validation against the local GBIF database works by using an initial SQL query using the
provided name and rank:

```SQL
SELECT * FROM backbone
    WHERE canonicalName = 'XX'
    AND taxonRank = 'YY';
```

This will return all exact matches for the combination, which can include rows for
multiple taxa with different description dates and authors, references and taxonomic
statuses. The `safedata_validator` package examines the returned rows and tries to
identify a single accepted, doubtful or non-canonical usage _in that order of
preference_.

The backbone table row for the identified taxon provides both `acceptedNameUsageID` and
`parentNameUsageID` fields that can be easily used in subsequent searches to get the
accepted name and hook the taxon into the higher taxonomy. Although canonical higher
taxon names are provided, their taxon IDs are not, so higher taxa are saved into a set
containing unique pairs of names and ranks and then added to the index when all taxa
have been processed.

### Problems

Rare edge cases include taxon names with two equally approved usages: for
example, the genus _Morus_ is an accepted usage for both mulberries and gannets.
This kind of problem is described in the JSON `"note"` field and provided to the
user. In these rare cases, an accepted usage would require a GBIF taxon ID to be
provided to discriminate between them.

## NCBI validation

The basic idea is to search taxon names in the NCBI database to look for matches, if
a match is found then the supplied rank and (optionally) NCBI taxa ID are used to
verify whether that this match is as expected. When searching taxa names that are
no longer considered to be the proper scientific name for a taxon, the validation
automatically maps onto the new taxonomic name. In cases where a two taxon are
considered to be equivalent NCBI removes the superseded taxon, and then records
the ID of this taxon in a separate table. This means that our validation method
can account for the entry of outdated taxa (and their IDs), which pass validation
with a warning being provided to the user.

We want three things from validation:

* To know whether the name provided by the user refers to an existing taxon.
* To know whether or not this taxon is superseded.
* The backbone taxonomic hierarchy for the taxa, so we can index at higher
  taxonomic levels. Taxa do not always neatly hook into the backbone taxonomic
  rank immediately above. However, this isn't too much of a problem in the NCBI
  case as parent IDs are always supplied so a full taxon hierarchy can always be
  constructed, which can then be pruned to form a backbone hierarchy.

Unlike NCBI does not provide status codes to taxa. As describes above it just
merges outdated taxa names and taxa IDs with up to date names and IDs. Our
validation proceed makes a note of where this happens. However, we do record
status codes so that GBIF and NCBI taxonomic information can be stored in the
same format. The three status codes we define are as follows:

```txt
accepted
merged
user
```

The first of these (accepted) is for taxa found using a currently valid name or
ID, the second of these (merged) is for valid taxa found using either a previously
valid name or ID that has been merged with another name ID. Finally, 'user' is used
for taxa that are not found in the NCBI database but have valid parent information.
This is particularly useful for potentially novel taxa, i.e. ones defined by the
user.

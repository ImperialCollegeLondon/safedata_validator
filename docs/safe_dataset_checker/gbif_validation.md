# GBIF validation notes

The basic idea is to search GBIF with a taxon name and rank to find if there is a match and then the
taxonomic status of the match. GBIF has the following status codes:

    doubtful, accepted, homotypic synonym, synonym,
    heterotypic synonym, proparte synonym, misapplied

Misapplications and synonyms always have a suggested **accepted** usage, doubtful taxa never do. However,
doubtful taxa typically still have a taxonomic hierarchy and can be encountered pretty much anywhere on
the rank walk to the root.

We want three things from validation:

* The status of the name provided by the user.
* The accepted usage (which might be the same).
* The backbone taxonomic hierarchy for the taxa, so we can index at higher taxonomic levels. However, this is complicated as there is no guarantee that taxa will always hook into the backbone at the next level above their own.

For example, the accepted taxon _Wanosuchus atresus_ is only hooked in at order level: it is accepted as
a species, but its parent is Crocodylia as the genus is doubtful. Similarly, _Goniopholis tenuidens_ is a
synonym at species level but again has Crocodylia as a parent (and is considered a synonym for the family
Goniopholidae).

  In some cases, taxa can hook in at levels below their own: the species _Brittonastrum greenei_ is a
  synonym of the **subspecies** _Agastache pallidiflora pallidiflora_.

  The `parent_taxon_level_analysis.R` file in this repository contains some code to check this:
  
  * All **accepted taxa** map to a more nested parent but 5% map to a more nested parent more than one
    step up the hierarchy. The table below shows child taxon level as rows and parent taxon level as
    columns.

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

  * Only 77% of **unaccepted taxa** map to a parent at the next most nested taxonomic level and 4.5% map
    to a parent at the same or a less nested level, as in the example above.

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

The two different APIs for offline and online usage work to try and get the same answer.

## Online use

The GBIF API provides the `species/match` endpoint, which allows searching for names against the GBIF
backbone taxonomy. More details are provided at the [GBIF developer
site](https://www.gbif.org/developer/species), but the endpoint is handy as it accepts variables to turn
off fuzzy matching (`strict=true`) and to check the name is at the stated rank (e.g. `rank=species`). For
example:

[http://api.gbif.org/v1/species/match?name=Panthera%20leo&strict=true&rank=species](http://api.gbif.org/v1/species/match?name=Panthera%20leo&strict=true&rank=species)

The JSON data returned looks like this:
```{json}
{"usageKey":5219404,
 "scientificName":"Panthera leo (Linnaeus, 1758)",
 "canonicalName":"Panthera leo",
 "rank":"SPECIES",
 "status":"ACCEPTED",
 "confidence":99,
 "matchType":"EXACT",
 "kingdom":"Animalia",
 "phylum":"Chordata",
 "order":"Carnivora",
 "family":"Felidae",
 "genus":"Panthera",
 "species":"Panthera leo",
 "kingdomKey":1,
 "phylumKey":44,
 "classKey":359,
 "orderKey":732,
 "familyKey":9703,
 "genusKey":2435194,
 "speciesKey":5219404,
 "synonym":false,
 "class":"Mammalia"}
```

With strict searching, here we have a single exact single match (`"matchType":"EXACT"`) of the name
_Panthera leo_ at species level. If no exact match is found at the rank then instead we would get
`"matchType":"NONE"` and possibly some `"notes"` with further details. Note that the `match` endpoint
checks names against multiple fields but with strict matching this is almost certainly going to be to the
canonical name field (scientific name minus authorship).

A search can bring in multiple matches, where a single canonical name is used in multiple synonyms: for
example _Zenicomus photuroides_ has one accepted use as _Zenicomus photuroides_ Thomson, 1868 and about
10 synonyms with different authors.

The GBIF `species/match` API automatically looks for a single 'best' match - it will return the single
taxon with that canonical name over a set of synonyms: the JSON data will contain `"status": "ACCEPTED"`.
The `verbose=true` argument to the API can used to also return those alternative matches but we don't use
that information here as we don't expect authorship details. If there are only non-accepted usages, then
the `status` will be something else: for example `SYNONYM` or `MISAPPLIED`.

For accepted taxa, this information is enough to get the information needed: we've got the accepted usage
and a set of backbone taxa. The only things to be wary of are potential missing levels and the fact that
taxon keys only go up to species level. In the example above the `usageKey` matches the `speciesKey`, but
the backbone includes subspecies, variety and form levels. We can hook subspecies in, since the
`speciesKey` is the least nested parent possible but variety and form do not have enough information to
build the full chain without further queries. These levels are not supported by the checker.

When GBIF returns a synonym, there is more to be done. The JSON data from the `species/match` endpoint
contains the `usageKey` for the search name. That _probably_ means that the highest taxon key provided is
the accepted taxon and the next highest is the parent. However, the program currently checks this using
the `species/{usageKey}` endpoint, which explictly contains an `acceptedKey` and a `parentKey` which
ought to match the previous information. For example:

[http://api.gbif.org/v1/species/8548608](http://api.gbif.org/v1/species/8548608)

## Offline use

Validation against the local GBIF database works slightly differently as it is a straight match against
the canonical name with no preference filtering. The equivalent database query is:

```{SQL}
SELECT * FROM backbone
    WHERE canonicalName = 'XX'
    AND taxonRank = 'YY';
```

This will return all exact matches, including the less favoured alternatives that GBIF omits, so that
preference is reimplemented by the program. The exact algorithm used by the API is unknown, so the
program just looks for a single accepted use and then falls back to other statuses. Note that currently
the SQL query is case sensitive where the API query is not.

In many ways, this is easier to use because the database rows contain both `acceptedNameUsageID` and
`parentNameUsageID` fields that can be easily used to get the accepted name and hook the taxon into the
higher taxonomy. Although canonical higher taxon names are provided, their taxon IDs are not, so higher
taxa are saved into a set containing unique pairs of names and ranks and then added to the index when all
taxa have been processed.

## Problems

Rare edge cases include taxon names with two equally approved usages: for example, the genus _Morus_ is
an accepted usage for both mulberries and gannets. This kind of problem is described in the JSON `"note"`
field and provided to the user. In these rare cases, an accepted usage would require a GBIF taxon ID to
be provided to discriminate between them.


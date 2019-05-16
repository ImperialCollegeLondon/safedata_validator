
# The Taxa worksheet

Many datasets will involve data taken from organisms, whether that is a count of the number of
individuals or measurement of a trait such as body length. In order to help us keep track of
taxa, all datasets using taxa **must** contain a Taxa spreadsheet, providing taxonomic
information.

Note that you must only provide details for taxa actually used in the data worksheets. This
ensures that the taxonomic index for a dataset is accurate and also double checks that it the
omission of a taxon from the data worksheets is not an error.

# Taxon validation

In order to help keep the taxonomy as clean as possible and to allow us to index the taxonomic
coverage of datasets, we will check all taxon names in Taxa worksheet against the GBIF backbone
taxonomy database. If you want to check your taxon names and ranks, then the search engine is
here:

[https://www.gbif.org/species/search](https://www.gbif.org/species/search)

No online taxonomy is ever going to be 100% up to date (or 100% agree with your taxonomic
usage!) but the GBIF backbone has very good taxonomic coverage and is well curated.

# Taxon table layout

The table format looks like this:

| Name | Taxon name | Taxon type | Taxon ID | Parent name | Parent type | Parent ID | Comments |
|---|---|---|---|---|---|---|---|
| Crematogaster borneensis | Crematogaster borneensis | Species |   |   |   |   |   |
| Dolichoderus sp. | Dolichoderus | Genus |   |   |   |   |   |
| Gannets	   | Morus | Genus | 2480962 |   |   |   |   |
| Morphospecies 1 | NA | Morphospecies |   | Formicidae | Family |   |   |

The table must contain column headers in the **first row** of the worksheet. The headers must
include:

 * **Name**: This column must contain a local name for **all** of the taxa that you are going to use in the rest of the dataset. You cannot have duplicated names! Note that these can be abbreviations or codes: if you want to use `Crbe` in your data worksheets, rather than typing out `Crematogaster borneensis` every time, then that is fine.

!!! Note

    These are the names that you are going to use in your data worksheet. The 
	other columns are to help us validate the taxonomy of your names.

 * **Taxon name**: This column must contain the scientific name of the taxon, which will be used for taxon validation via GBIF. Note that this should not include any taxon authority, so _Panthera tigris_ not _Panthera tigris_ (Linnaeus, 1758).

 * **Taxon type**: This column must provide the **taxonomic type** of the named taxon, which is usually the taxonomic **rank**. For example, the taxon _Pongo pygmaeus_ would be of type `Species` and the taxon Formicidae would be of type `Family`. There are two additional possible values which are described in more detail below: `Morphospecies` and `Functional group`.

 * **Taxon ID**: This column is optional - you don't have to include it and you don't have to fill in every row if it is included. It is **only needed** if the taxon name and rank are ambiguous.  
     For example, the genus _Morus_ can refer to either mulberries or gannets. In these rare cases, you will need to look up the taxon in GBIF and provide the GBIF ID number. The example in the table allows us to confirm that you mean this _Morus_: [https://www.gbif.org/species/2480962](https://www.gbif.org/species/2480962).

There are three types of taxa, described below, that can't be validated directly against GBIF.
In these cases, the next three columns (**Parent name**, **Parent type** and **Parent
ID**) are used to provide a parent taxon that can be validated and provides a taxonomic hook
to allow us to place the taxon in the backbone taxonomy. They are used in exactly the same way
as the Taxon columns: the only restriction is that the Parent type must be one of seven
core ranks used in the GBIF backbone. Again **Parent ID** is only needed to resolve ambiguous 
taxon names.

The Parent columns only need to be filled in for taxa falling in one of the following
groups.

## Comments

The comments field is entirely optional - if you have a fairly standard set of taxa with no
 serious issues then you can leave it out entirely or it can be empty. If you do have particular
 notes that you want to make - strong disagreements with GBIF taxonomy, new species notes and the 
 like - then these can be very useful for future researchers trying to place taxa. 

# New and unrecognized taxa

If a taxon is new or not recognized by GBIF (and you're sure you're right!) then provide a
parent name and type to allow us to hook the taxon into the index. For example, before _Pongo
tapanuliensis_ was recognised as a species, you would need to provide _Pongo_ as a parent
name of type 'genus' allows us to place the new taxon.

Note that you should still provide the name and taxonomic level as usual, and the parent 
information will allow us to place that taxon in the tree. A comment is useful in these cases
but not mandatory.

| Name | Taxon name | Taxon type | Taxon ID | Parent name | Parent type | Parent ID | Comments |
|---|---|---|---|---|---|---|---|
| lost_orang | Pongo tapanuliensis| Species |   |  Pongo| Genus   |   |  New species |

You would then get a message in the file validation report saying:

    - Row 1 (lost_orang): not found in GBIF but has valid parent information


# Morphospecies and Functional groups

For morphospecies and functional groups, the taxon name is the label to be used in the dataset.
Set the Scientific name to be 'NA' - it cannot be blank - and then specify the taxon type as
'Morphospecies' or 'Functional group'.

Now you need to provide a parent taxon and type. The level of taxonomic certainty for
morphospecies and functional groups is quite variable, but we'd like the finest taxonomic level
you can provide. As an example, in the table above, 'Morphospecies #1' is simply identified as
being an ant.

# Less common taxonomic levels

The GBIF backbone taxonomy only includes the following eight major levels: Kingdom, Phylum,
Order, Class, Family, Genus, Species and Subspecies. If you need to use taxa defined at any
intermediate levels, then again provide a parent taxon and type. For example, if you were
counting bees and only identifying to tribe level (Bombini, Euglossini, etc.) then the parent
family Apidae would allow us to hook the taxa into the backbone taxonomy. The subfamily Apinae
would be more precise, but subfamily isn't one of the backbone taxonomic levels.

# My data doesn't contain taxa

Fine. You can omit the Taxa worksheet!
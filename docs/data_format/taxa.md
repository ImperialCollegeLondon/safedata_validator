
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

| Name | Taxon name | Taxon type | Taxon ID | Ignore ID | Parent name | Parent type | Parent ID | Comments |
|---|---|---|---|---|---|---|---|---|
| crbo | Crematogaster borneensis | Species |   |   |   |   |   |   |
| dolic_sp | Dolichoderus | Genus |   |   |   |   |   |   |
| gannets  | Morus | Genus | 2480962 |   |   |   |   |   |
| lost_orang | Pongo tapanuliensis | Species |   |   | Pongo | Genus  |   |  New species |
| morpho1 | NA | Morphospecies |   |   | Formicidae | Family |   |   |
| bombines | Bombini | Tribe |   |   | Apidae | Family |   |   |
| micr_hid | Microcopris hidakai | Species |   | 1090433 | Microcopris | Genus |   |   |

The table must contain column headers in the **first row** of the worksheet. The first three columns (Name, Taxon name and Taxon type) are mandatory and contain the following:

 * **Name**: This column must contain a local name for **all** of the taxa that you are going to use in the rest of the dataset. You cannot have duplicated names! Note that these can be abbreviations or codes: if you want to use `Crbe` in your data worksheets, rather than typing out `Crematogaster borneensis` every time, then that is fine.

!!! Note

    These are the names that you are going to use in your data worksheet. The
	other columns are to help us validate the taxonomy of your names.

 * **Taxon name**: This column must contain the scientific name of the taxon, which will be used for taxon validation via GBIF. Note that this should not include any taxon authority, so _Panthera tigris_ not _Panthera tigris_ (Linnaeus, 1758).

 * **Taxon type**: This column must provide the **taxonomic type** of the named taxon, which is usually the taxonomic **rank**. For example, the taxon _Pongo pygmaeus_ would be of type `Species` and the taxon Formicidae would be of type `Family`. If the taxon type is not one of the core GBIF backbone levels ( Kingdom, Phylum, Order, Class, Family, Genus, Species and Subspecies), or is one of the special values `Morphospecies` or `Functional group`, then you will need to provide a **parent taxon** (see below).

* **Comments**: This is entirely optional - if you have a fairly standard set of taxa with no  serious issues then you can leave it out entirely or it can be empty. If you do have particular notes that you want to make - explaining disagreements with GBIF taxonomy, new species notes and the  like - then these can be very useful for future researchers trying to place taxa.

The other fields are optional and are used to handle the following taxonomic issues. You only need to put anything in these optional fields for the rows they apply to: leave them blank for all other rows as in the example table above.

## Ambiguous taxon names

Some taxon names map to more than one taxon.  For example, the genus _Morus_ can refer to either mulberries or gannets. The checking code will raise an error if it encounters an ambiiguous name.  In these rare cases, you will need to look up the taxon in GBIF and provide the GBIF ID number in the **Taxon ID** field.

The example in the table allows us to confirm that you mean this _Morus_: [https://www.gbif.org/species/2480962](https://www.gbif.org/species/2480962).


## Incorrect GBIF taxonomy

 In general, we follow the GBIF 'accepted usage' of a taxon name and this includes following when GBIF say a taxon is a synonym. We include your taxon name in our taxonomic indices, but under the taxonomic hierarchy suggested by GBIF.  

 In the example in the table, GBIF asserts that _Microcopris hidakai_ is a synonym of [_Onthophagus hidakai_](https://www.gbif.org/species/1090433) and so our taxonomic index will include both those names under the parent genus _Onthophagus_. The taxon validation of your dataset does warn about this in the output:

    ? Row 7 (micr_hid): considered a homotypic_synonym of Onthophagus hidakai (1090433) in GBIF backbone

If you strongly disagree with that accepted usage, then you can use the **Ignore ID** field to explicitly ignore the GBIF match. You will need to insert the GBIF ID number of the accepted usage (reported in the warning message) in this field. You will then need to provide a **parent taxon** (see below)

# Parent taxon details

If a taxon name cannot be matched against GBIF,  then you will need to provide details of a parent taxon. This is a taxon that can be validated in GBIF and provides a hook to allow us to place the child taxon in the backbone taxonomy. Note that you should still provide the taxon information as usual.

In these cases, the next three columns (**Parent name**, **Parent type** and **Parent
ID**) are used to provide a parent taxon. They are used in exactly the same way
as the Taxon columns: the only restriction is that the Parent type must be one of
core ranks used in the GBIF backbone: Kingdom, Phylum, Order, Class, Family, Genus, Species and Subspecies. Again **Parent ID** is only needed to resolve ambiguous
taxon names.

The table shows examples for the following cases:

## New and unrecognized taxa

If a taxon is new or not recognized by GBIF (and you're sure it isn't just a typo!) then provide a
parent name and type to allow us to hook the taxon into the index. For example, before _Pongo
tapanuliensis_ was recognised as a species, you would need to provide _Pongo_ as a parent
name of type 'genus',

You would then get a message in the file validation report saying:

    - Row 4 (lost_orang): not found in GBIF but has valid parent information


## Morphospecies and Functional groups

For morphospecies and functional groups, the taxon name is the label to be used in the dataset.
Set the Scientific name to be 'NA' - it cannot be blank - and then specify the taxon type as
'Morphospecies' or 'Functional group'.

Now you need to provide a parent taxon and type. The level of taxonomic certainty for
morphospecies and functional groups is quite variable, but we'd like the finest taxonomic level
you can provide. As an example, in the table above, 'Morphospecies #1' is simply identified as
being an ant. The validation report will show:

    - Row 5 (morpho1): Morphospecies with valid parent information

## Less common taxonomic levels

If you want to use taxa defined at any intermediate taxonomic levels that are not included in the GBIF backbone, then you will again need to provide a parent taxon. In the example in the table, taxa identified to tribe level (Bombini) are hooked in as being in the family Apidae. The subfamily Apinae would be more precise, but subfamily isn't one of the backbone taxonomic levels.

## Ignored taxa

If you have set an **Ignore ID** value for a taxon, then that is explicitly rejecting the parent taxon that we would naturally use from GBIF and so you will need to provide a replacement. The example in the table is explicitly saying that _Microcopris hidakai_ is not a synonym but is a separate taxon that belongs in _Microcopris_.

# My data doesn't contain taxa

Fine. You can omit the Taxa worksheet!

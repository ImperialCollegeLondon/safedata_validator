# The GBIFTaxa worksheet

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

This worksheet exists to record taxonomic in information for organisms referred to in
the Data worksheets. This taxonomic information should be recorded in a format that can
be validated against the GBIF backbone taxonomy. Generally the GBIF backbone taxonomy is
most suitable for taxa discovered via direct observation (rather than sequencing), so
this worksheet is most appropriate for storing the details of directly observed taxa. If
taxa are used anywhere in the dataset either this worksheet or the [NCBITaxa
worksheet](./ncbi_taxa.md) must be included. It is also an option to provide both a
GBIFTaxa worksheet and a NCBITaxa worksheet, e.g. in cases where both sequencing and
direct observational data are being reported.

## GBIF Taxon validation

In order to help keep the taxonomy as clean as possible and to allow us to index the
taxonomic coverage of datasets, we will check all taxon names in GBIFTaxa worksheet
against the GBIF backbone taxonomy database. If you want to check your taxon names and
ranks, then the search engine is here:

[https://www.gbif.org/species/search](https://www.gbif.org/species/search)

No online taxonomy is ever going to be 100% up to date (or 100% agree with your taxonomic
usage!) but the GBIF backbone has very good taxonomic coverage and is well curated.

## Taxon table layout

The table format looks like this:

<!-- markdownlint-disable MD013 -->
{{ read_excel('Example.xlsx', sheet_name = 'GBIFTaxa', keep_default_na = False, tablefmt = 'github') }}
<!-- markdownlint-enable MD013 -->

The table must contain column headers in the **first row** of the worksheet. The first
three columns (Name, Taxon name and Taxon type) are mandatory and contain the following:

- **Name**: This column must contain a local name for **all** of the taxa that you are
  going to use in the rest of the dataset. You cannot have duplicated names! Note that
  these can be abbreviations or codes: if you want to use `Crbe` in your data
  worksheets, rather than typing out `Crematogaster borneensis` every time, then that is
  fine.

!!! Note
    These are the names that you are going to use in your data worksheet. The
    other columns are to help us validate the taxonomy of your names.

- **Taxon name**: This column must contain the scientific name of the taxon, which will
  be used for taxon validation via GBIF. Note that this should not include any taxon
  authority, so _Panthera tigris_ not _Panthera tigris_ (Linnaeus, 1758).

- **Taxon type**: This column must provide the **taxonomic type** of the named taxon,
  which is usually the taxonomic **rank**. For example, the taxon _Pongo pygmaeus_ would
  be of type `Species` and the taxon Formicidae would be of type `Family`. If the taxon
  type is not one of the core GBIF backbone ranks ( Kingdom, Phylum, Order, Class,
  Family, Genus, Species and Subspecies), or is one of the special values
  `Morphospecies` or `Functional group`, then you will need to provide a **parent
  taxon** (see below).

- **Comments**: This is entirely optional - if you have a fairly standard set of taxa
  with no  serious issues then you can leave it out entirely or it can be empty. If you
  do have particular notes that you want to make - explaining disagreements with GBIF
  taxonomy, new species notes and the  like - then these can be very useful for future
  researchers trying to place taxa.

The other fields are optional and are used to handle the following taxonomic issues. You
only need to put anything in these optional fields for the rows they apply to: leave
them blank for all other rows as in the example table above.

### Ambiguous taxon names

Some taxon names map to more than one taxon.  For example, the genus _Morus_ can refer
to either mulberries or gannets. The checking code will raise an error if it encounters
an ambiguous name.  In these rare cases, you will need to look up the taxon in GBIF and
provide the GBIF ID number in the **Taxon ID** field.

The example in the table allows us to confirm that you mean this _Morus_:
[https://www.gbif.org/species/2480962](https://www.gbif.org/species/2480962).

### Incorrect GBIF taxonomy

In general, we follow the GBIF 'accepted usage' of a taxon name and this includes
following when GBIF say a taxon is a synonym. We include your taxon name in our
taxonomic indices, but under the taxonomic hierarchy suggested by GBIF.

If you strongly disagree with that accepted usage, then you can use the **Ignore ID**
field to explicitly ignore the GBIF match. You will need to insert the GBIF ID number of
the accepted usage (reported in the warning message) in this field. You will then need
to provide a **parent taxon** (see below)

In the example in the table, a (fictional) recently discovered species of gannet
([_Morus_](https://www.gbif.org/species/2480962)) has been included. Unfortunately, this
species shares a scientific name ([_Morus rubra_](https://www.gbif.org/species/5361886))
with a species of mulberry ([_Morus_](https://www.gbif.org/species/2984545)). This
results in information about the mulberry species being incorrectly included in the
taxonomic hierarchy. To stop this from happening, information on the correct parent
genus must be included, and the incorrect species match must be explicitly ignored using
the **Ignore ID** field.

The example in the table is made up, but similar situations have arisen in previous
datasets published using the `safedata` system. That said conflicts of this type are
rare, so you are unlikely to have to ever use **Ignore ID**.

## Parent taxon details

If a taxon name cannot be matched against GBIF,  then you will need to provide details
of a parent taxon. This is a taxon that can be validated in GBIF and provides a hook to
allow us to place the child taxon in the backbone taxonomy. Note that you should still
provide the taxon information as usual.

In these cases, the next three columns (**Parent name**, **Parent type** and **Parent
ID**) are used to provide a parent taxon. They are used in exactly the same way as the
Taxon columns: the only restriction is that the Parent type must be one of core ranks
used in the GBIF backbone: Kingdom, Phylum, Order, Class, Family, Genus, Species and
Subspecies. Again **Parent ID** is only needed to resolve ambiguous taxon names.

The table shows examples for the following cases:

### New and unrecognized taxa

If a taxon is new or not recognized by GBIF (and you're sure it isn't just a typo!) then
provide a parent name and type to allow us to hook the taxon into the index. For
example, before _Pongo tapanuliensis_ was recognised as a species, you would need to
provide _Pongo_ as a parent name of type 'genus',

You would then get a message in the file validation report saying:

> \- Row 4 (lost_orang): not found in GBIF but has valid parent information

### Morphospecies and Functional groups

For morphospecies and functional groups, the **Taxon type** should be specified as
'Morphospecies' or 'Functional group'. The **Name** field should be filled in with
whatever name is used to label it in the dataset. In this case, the **Taxon name** field
is not checked or used anywhere in the validation process, so can be filled in with
whatever you wish. However, this the field still **has** to be filled out. If you do not
wish to provide a name here, we recommend just filling it in with 'NA'.

Now you need to provide a parent taxon and type. The level of taxonomic certainty for
morphospecies and functional groups is quite variable, but we'd like the finest
taxonomic level you can provide. As an example, in the table above, 'Morphospecies #1'
is simply identified as being an ant. The validation report will show:

> \- Row 5 (morpho1): Morphospecies with valid parent information

### Less common taxonomic levels

If you want to use taxa defined at any intermediate taxonomic levels that are not
included in the GBIF backbone, then you will again need to provide a parent taxon. In
the example in the table, taxa identified to tribe level (Bombini) are hooked in as
being in the family Apidae. The subfamily Apinae would be more precise, but subfamily
isn't one of the backbone taxonomic levels.

### Ignored taxa

If you have set an **Ignore ID** value for a taxon, then that is explicitly rejecting
the parent taxon that we would naturally use from GBIF and so you will need to provide a
replacement. The example in the table is explicitly saying that _Morus rubra_ is not a
species of mulberry (genus [_Morus_](https://www.gbif.org/species/2984545)) but is a
species of gannet that belongs in the _Animalia_ genus
[_Morus_](https://www.gbif.org/species/2480962).

## Deleted taxa

The GBIF backbone also includes a large number of **deleted** taxon ids: these are a mix
of duplicated names, typos in the database and other errors. The GBIF ID of these taxa
are preserved, along with _some_ information, but we do not allow deleted taxa to be
used in the GBIFTaxa worksheet.

## My data doesn't contain taxa

Fine. You can omit either or both of the GBIFTaxa and NCBITaxa worksheets!

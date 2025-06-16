# The Taxa worksheets

Many datasets will involve data taken from organisms, whether that is a count of the
number of individuals or measurement of a trait such as body length. In order to help us
keep track of taxa, all datasets using taxa **must** contain a Taxa spreadsheet,
providing taxonomic information. There are two kinds:

* A [GBIFTaxa worksheet](gbif_taxa.md), which is typically used for observational data
* [Sequenced taxa worksheets](sequenced_taxa.md), which are used for data that has been
  obtained from genetic sequencing, the taxonomies of which can be generated from a wide
  variety of reference taxonomy databases

Both worksheet types can be provided for a single data set - for example, to capture
the species of trapped mammals and the sequences associated with their stomach contents.

Only a single GBIFTaxa worksheet can be provided, but multiple sequenced taxa worksheets
can be included. Sequenced taxonomies should only be split over multiple sheets when
they have be obtained via different sequencing approaches (e.g. 16S vs 18S) that require
them to be compared against either different taxonomy databases or different versions of
the same taxonomy database. Information about the database and version used for each
sequenced taxa worksheet **must** be provided in [the relevant part of the
summary](summary.md#sequenced-taxa-sheets-block) .

Note that both GBIF and sequenced taxa worksheets are expected to provide details on all
the taxa referenced in the dataset and **only the taxa referenced in the dataset**. This
is to ensures that the taxonomic index for a dataset is accurate and also double checks
that the omission of a taxon from the data worksheets is not an error.

## Duplicate taxa

You must not include duplicate taxa in a single taxon worksheet but it is perfectly
possible for a taxa to be identified both through field observation and through
sequencing, potentially even multiple times through different sequencing datasets.

When this happens, we require that you use a **different taxon worksheet name**
in the different taxon worksheets. This is to ensure that every taxon referred to in the
Data worksheets can be unambiguously linked back to a taxonomic data source. As an
example, if your data include both field observations of *Vulpes vulpes* and
identification of fox faecal samples through sequencing, you could use `v_vulpes` for
the GBIFTaxa worksheet and `v_vulpes_faeces` in the sequenced taxon worksheet.

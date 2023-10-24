# The Taxa worksheets

Many datasets will involve data taken from organisms, whether that is a count of the
number of individuals or measurement of a trait such as body length. In order to help us
keep track of taxa, all datasets using taxa **must** contain a Taxa spreadsheet,
providing taxonomic information. There are two kinds:

* A [GBIFTaxa worksheet](gbif_taxa.md), which is typically used for observational data
* A [NCBITaxa worksheet](ncbi_taxa.md), which is used for data that has been obtained
  from genetic sequencing which very often cannot be easily mapped onto the GBIF taxonomic
  backbone.

Both worksheet types can be provided for a single data set - for example, to capture
the species of trapped mammals and the sequences associated with their stomach contents.

A taxon can be defined in both sheets (i.e. when a taxon has been both directly observed
and discovered in sequencing data). However, in this case the taxon must be given a
different name in each Taxa worksheet, e.g. `v_vulpes_ob` and `v_vulpes_seq` for
observed and sequenced instances of *Vulpes vulpes*, respectively. This is to ensure
that every taxon referred to in the Data worksheets can be unambiguously linked back to
a data source.

Note that you must only provide details for taxa actually used in the data worksheets.
This ensures that the taxonomic index for a dataset is accurate and also double checks
that the omission of a taxon from the data worksheets is not an error.

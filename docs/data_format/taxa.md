# The Taxa worksheets

Many datasets will involve data taken from organisms, whether that is a count of the number of individuals or measurement of a trait such as body length. In order to help us keep track of taxa, all datasets using taxa **must** contain a Taxa spreadsheet, providing taxonomic information. This can be a [GBIFTaxa worksheet](gbif_taxa.md) or a [NCBITaxa worksheet](ncbi_taxa.md). In general the GBIFTaxa worksheet is better for observational data, whereas the NCBITaxa worksheet is better for sequencing data. For this reason both sheets can be provided in a single data set, and specific taxa can be defined in both sheets (i.e. when a taxon has been both directly observed and discovered in sequencing data). However, in this case the taxon must be given a different name in each Taxa worksheet, e.g. `v_vulpes_ob` and `v_vulpes_seq` for observed and sequenced instances of *Vulpes vulpes*, respectively. This is to ensure that the every taxon referred to in the Data worksheets can be unambiguously linked back to a data source.

Note that you must only provide details for taxa actually used in the data worksheets. This
ensures that the taxonomic index for a dataset is accurate and also double checks that it the
omission of a taxon from the data worksheets is not an error.

# Frequently asked questions

The error messages from the data formatting checker can be a bit terse, so here are some
common questions.

## GBIFTaxa worksheet

### XXX found in GBIF backbone, but additional parent information is incompatible

You only have to provide parent information when GBIF can't find the taxon itself. If
you do provide parent information for a taxon that is in GBIF, then we check the two are
compatible. This isn't usually a problem but it can go wrong if the taxon is considered
a synonym.

For example, the fish species _Barbodes sealei_ is considered by GBIF to be a synonym
for _Puntius sealei_, and the taxonomic hierarchy we get from GBIF is that the species
is in the genus _Puntius_. So, if you put _Barbodes_ as the parent genus, you will get
an error because _Barbodes_ is not in that taxonomic hierarchy. We do not take a view on
whether GBIF is correct - we record both synonyms and canon names in the taxon index for
the dataset.

The solution is to remove the Parent information - you would then get a warning about
the synonymy but no error.

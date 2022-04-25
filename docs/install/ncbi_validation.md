# NCBI validation notes

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

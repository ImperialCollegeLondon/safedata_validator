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

Unlike NCBI does not provide status codes to taxa. As describes above it just
merges outdated taxa names and taxa IDs with up to date names and IDs. Our
validation proceed makes a note of where this happens. However, we do record
status codes so that GBIF and NCBI taxonomic information can be stored in the
same format. The three status codes we define are as follows:

    accepted, merged, user

The first of these (accepted) is for taxa found using a currently valid name or
ID, the second of these (merged) is for valid taxa found using either a previously
valid name or ID that has been merged with another name ID. Finally, 'user' is used
for taxa that are not found in the NCBI database but have valid parent information.
This is particularly useful for potentially novel taxa, i.e. ones defined by the
user.

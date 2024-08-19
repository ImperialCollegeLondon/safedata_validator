# The `safedata_validator` API

The `safedata_validator` package is divided into submodules handling different
components of the validation process.

* The [`summary` module](./summary.md) handles the parsing and validation of the
  [Summary datasheet](../../data_providers/data_format/summary.md)
* The [`locations` module](./locations.md) handles the loading of the locations
  gazetteer and providing a Locations class to be used in validation of the
  [Locations datasheet](../../data_providers/data_format/locations.md).
* Similarly, the [`taxa` module](./taxa.md) handles the creation of taxonomic validation
  classes, using local databases to check the contents of the
  [GBIF](../../data_providers/data_format/gbif_taxa.md) and
  [NCBI](../../data_providers/data_format/ncbi_taxa.md) data worksheets.
* The [`field` module](./field.md) contains all of the code for reading and validating
  data worksheets, along with the main Dataset class used to load and validate entire
  dataset files.
* The [`resources` module](./resources.md) is used to load and validate the
  [configuration file](../../data_managers/install/configuration.md).
* The [`extent` module](./extent.md) is used to configure and track the [temporal and
  spatial
  extents](../../data_managers/install/configuration.md#validation-configuration) of a
  dataset.
* The [`logger` module](./logger.md) is used to set up logging of the validation
  process. Logging is a core component of the system, because the validation process is
  written to work through the whole file, logging issues as it goes, rather than exiting
  at the first problem.
* The [`taxondb` module](./taxondb.md) is used to download copies of the GBIF and NCBI
  taxonomy databases and build usable local SQLite3 databases from them.
* The [`zenodo` module](./zenodo.md) is used to communicate with the Zenodo in order to
  create and upload Zenodo deposits.

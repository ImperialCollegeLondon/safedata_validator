# Local NCBI taxonomy database

The NCBI taxonomy database is composed of a number of tables. The
present state of these tables are provided online as DMP files. To use
this data as a local data resource for taxon validation, the `safedata_validator`
package requires it to be loaded into a singles SQLite3 database.

You may need to download and install SQLite3. You want the 'bundle' installer
for your operating system.

[https://www.sqlite.org/download.html](https://www.sqlite.org/download.html)

Next, go to this page:

[https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)

The DMP files are provided in three different zip formats (`taxdmp.zip`, `taxdump.tar.Z`,
`taxdump.tar.gz`), you can download which ever of these is most convient to unzip.
This zip file contains the DMP files for a [number of tables](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt). The ones that we require are `nodes.dmp`, `names.dmp`, and `merged.dmp`, which
describe the taxonomic nodes, the names assigned to specific nodes, and superseeded nodes that have
now been merged into other nodes, respectively.

There are a number of steps needed to convert this data into a SQLite3 database.
For convenience, this conversion has been automated into a Python script
(`additional_scripts/create_ncbi_database_sqlite.py`). This may need updating if file
names and structures change, but the basic process should work.

* The nodes table contains a lot of fields unrelated to anything we are interested in (i.e. `division_id`). We therefore drop all fields other than `tax_id`, `parent_tax_id`, `rank`, and `comments`.

* The names table includes a field (`unique_name`) that gives a unique version of every name. Given that `safedata_validator` has built in synonym handling and that these unique names extremely specific (e.g. `Bacteria <bacteria>`), we do not believe this is a useful field for our purpose. It is, therefore, dropped.

* A single database is then generated with three tables: `nodes`, `names`, `merged`.

* The speed of the package is much improved by building covering indices to
   speed up the four kinds of searches used by `safedata_validator`:

      1.  searches for specific taxon ids in the nodes table,
      2.  searches for taxon names in the names table,
      3.  searches for scientific names on the names tables and
      4.  searches for specific taxon ids in the merged table.

Once you have an SQLite3 database, you will then need to edit the `ncbi_database`
entry in your [configuration file](configuration.md) to provide the path to your
new SQLite file or provide the path as an argument.

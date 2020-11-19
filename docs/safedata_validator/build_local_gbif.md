# Local GBIF backbone database

GBIF provide the backbone taxonomy data as a SQL file defining the table structure and a tab delimited text file containing the backbone data. To use this data as a local data resource for taxon validation, the `safedata_validator` package requires it to be loaded into a SQLite3 database.

You will need to do the following:

* You may need to download and install SQLite3. You want the 'bundle' installer for your operating system.

      [https://www.sqlite.org/download.html](https://www.sqlite.org/download.html)

* Download the current zipped backbone from the GBIF backbone archive:

      [http://rs.gbif.org/datasets/backbone/readme.html](http://rs.gbif.org/datasets/backbone/readme.html)

      You want the `current-backbone-simple.txt.gz` file. The following commands will download and extract it, but you can equally do this manually.

        curl -O  https://hosted-datasets.gbif.org/datasets/backbone/backbone-current-simple.txt.gz
        gunzip backbone-current-simple.txt.gz

* You also need to download the SQL table definition for the file from this link:

        curl -O https://raw.githubusercontent.com/gbif/checklistbank/master/checklistbank-mybatis-service/src/main/resources/backbone-ddl.sql

      At present, that file is missing a final semi-colon after the right bracket at the bottom, so open it and change `)` to `);`.

* The file `backbone-current-simple.txt` contains the field `name_published_in`. The values in this field include a lot of quotes, which are tricky to parse into an SQLite database. The field isn't used in this application, so the simplest thing to do is to delete it.

        cut -f 1-28,30 backbone-current-simple.txt > backbone-current-simple-truncate.txt

      You also now have to remove the line in `backbone-ddl.sql` that defines the field: `name_published_in text,`.

* Now, you can create a new database with that table:

        sqlite3 backbone-current-simple.sqlite < backbone-ddl.sql

* Now open up the new SQLite database and import the data. This will take a while as the file is large.

        sqlite3 backbone-current-simple.sqlite
        .mode tab
        .import backbone-current-simple-truncate.txt backbone

*  The file contains a lot of `\N` values, which is a PostgreSQL symbol for a null field. SQLite 3 treats these as strings, so they need to be reset to `null` for each field that contains them. It might seem easier to simple delete all the `\N` values in the file to leave empty fields but SQLite3 then imports these as empty strings, not as null.

        update backbone set id = null where id = '\N';
        update backbone set parent_key = null where parent_key = '\N';
        update backbone set basionym_key = null where basionym_key = '\N';
        update backbone set is_synonym = null where is_synonym = '\N';
        update backbone set status = null where status = '\N';
        update backbone set rank = null where rank = '\N';
        update backbone set nom_status = null where nom_status = '\N';
        update backbone set constituent_key = null where constituent_key = '\N';
        update backbone set origin = null where origin = '\N';
        update backbone set source_taxon_key = null where source_taxon_key = '\N';
        update backbone set kingdom_key = null where kingdom_key = '\N';
        update backbone set phylum_key = null where phylum_key = '\N';
        update backbone set class_key = null where class_key = '\N';
        update backbone set order_key = null where order_key = '\N';
        update backbone set family_key = null where family_key = '\N';
        update backbone set genus_key = null where genus_key = '\N';
        update backbone set species_key = null where species_key = '\N';
        update backbone set name_id = null where name_id = '\N';
        update backbone set scientific_name = null where scientific_name = '\N';
        update backbone set canonical_name = null where canonical_name = '\N';
        update backbone set genus_or_above = null where genus_or_above = '\N';
        update backbone set specific_epithet = null where specific_epithet = '\N';
        update backbone set infra_specific_epithet = null where infra_specific_epithet = '\N';
        update backbone set notho_type = null where notho_type = '\N';
        update backbone set authorship = null where authorship = '\N';
        update backbone set year = null where year = '\N';
        update backbone set bracket_authorship = null where bracket_authorship = '\N';
        update backbone set bracket_year = null where bracket_year = '\N';
        update backbone set issues = null where issues = '\N';

*  Now enter these two commands to build two indices to create covering indices to speed up the two kinds of searches used by `safedata_validator`: i) searches on the canonical name  and rank of a taxon and ii) searches on the taxon id. That is the last step so then quit from SQLite3.

        create index backbone_name_rank on backbone (canonical_name, rank);
        create index backbone_id on backbone (id);
        .quit

* You can now delete the `backbone-ddl.sql`, `backbone-current-simple.txt` and `backbone-current-simple-truncate.txt` file.

* Edit the `gbif_database` entry in your configuration file to provide the path to your new SQLite file or provide the path as an argument.


# Installing the SAFE Dataset Checker


The following steps should allow you to run the checks yourself before submitting a dataset:

1. If you're using a Mac or Linux, then it is almost certain that you already have Python installed. You may also do if you are using Windows. To find out, open a command line window (run Terminal on a Mac, or `cmd` on Windows) and then type `python` (on Mac or Unix) or `python.exe` on Windows).

2. If some text and a prompt (`>>>`) appears then you have Python. Check the version number in the first line: if it doesn't start Python 2.7 then you currently need to install Python 2.7 to run the checker. You can have multiple versions of Python installed, but it is going to be more complicated than is covered here.

    If you've got Python 2.7 then type `quit()` and skip to step 4.
    
    If you get a line that says the command is not found then you need to install Python. Download a copy from here:

    [https://www.python.org/downloads/](https://www.python.org/downloads/)

    The code is currently written to use Python 2.7, so make sure you download an installer for Python 2.7.14 and not the more recent Python 3.6 versions. For Windows, choose one of the MSI installer options.

3. Repeat the command line check from the first step: if this still doesn't work then you probably just need to tell the computer where to find Python: search online for instructions to add python to the `PATH` environment variable. On Windows, you will want to add (using the typical install location) `C:\Python27` and `C:\Python27\scripts`.

4. The `safe_dataset_checker` program mostly uses commands from the Python Standard Library - a set of code packages that are installed with Python - but does use four extra packages that can be installed using the `pip` package installer. At the command line, type:

        pip install openpyxl requests simplejson shapely

    Those packages allow Python to: read Excel files, get validation data over the internet, handle JSON formatted data and validate WKT formatted  GIS vector data.

5. Now create a folder to keep your data checking code in and download the `safe_file_checker.py` file into it from [here](https://raw.githubusercontent.com/ImperialCollegeLondon/safe_dataset_checker/master/safe_dataset_checker.py)

6. Now open a command line terminal, change to your data checking directory and run the following:

        python ./safe_dataset_checker.py -h

    In Windows, you will need to change `python` to `python.exe` here and in commands below that start `python`. You should see the usage  instructions.

Note that you can also install `safe_dataset_checker` as a package, which will allow you to use the functions and classes in the package within your own code.

## Local GBIF database

If you want to run the checker offline then you will need to download a copy of the backbone taxonomy and build a SQLite3 database from it. Using a local database is also much faster than using the GBIF API online. This isn't particularly hard, but the resulting database is around 1.6GB, so you'll need file space! The steps are:

* You may need to download and install SQLite3. You want the 'bundle' installer for your operating system.

      [https://www.sqlite.org/download.html](https://www.sqlite.org/download.html)

* Download the current zipped backbone from the GBIF backbone archive:

      [http://rs.gbif.org/datasets/backbone/readme.html](http://rs.gbif.org/datasets/backbone/readme.html)

      You want the `current-backbone-simple.txt.gz` file. The following commands will download and extract it, but you can equally do this manually.

        curl -O  http://rs.gbif.org/datasets/backbone/backbone-current-simple.txt.gz
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

*  Now enter these two commands to build two indices to create covering indices to speed up the two kinds of searches used by `safe_dataset_checker.py`: i) searches on the canonical name  and rank of a taxon and ii) searches on the taxon id. That is the last step so then quit from SQLite3.

        create index backbone_name_rank on backbone (canonical_name, rank);
        create index backbone_id on backbone (id);
        .quit

* You can now delete the `backbone-current-simple.txt` and `backbone-current-simple-truncate.txt` file. You should now be able to use local taxonomy validation by telling `safe_dataset_checker.py` where to find the SQLite3 database:

        python safe_dataset_checker.py -g backbone-current-simple.sqlite My_Excel_File.xlsx


## Local gazetteer database

The current list of valid location names and bounding boxes is automatically downloaded directly from the SAFE Gazetteer. To work offline, get a copy of this data from the following link:
	
[https://www.safeproject.net/call/json/get_locations_bbox](https://www.safeproject.net/call/json/get_locations_bbox)

Save the output as a file in your data checker folder (e.g. `SAFE_locations.json`). You will then be able to run the program using the following:

    python safe_dataset_checker.py -l SAFE_locations.json My_Excel_File.xlsx


## Fully offline use

 If you've done both the above steps then the following example would validate a file using both local datasets, and won't need the internet at all.

    python safe_dataset_checker.py -g backbone-current.sqlite -l SAFE_locations.json My_Excel_File.xlsx

You might want to create a shortcut to this command. On Windows, create a file called `local_dataset_check.cmd` with the following contents:

     python.exe safe_dataset_checker.py -g backbone-current.sqlite -l SAFE_locations.json %1

You can now run the program using the local checking from a terminal by typing:

     local_data_check.cmd My_Excel_File.xlsx

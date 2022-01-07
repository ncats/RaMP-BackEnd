# RaMP Creation Package #

This repository contains the source code that was used to create RaMP - a Relation database of Metabolic Pathways (https://www.ncbi.nlm.nih.gov/pubmed/29470400).  RaMP is a conglomerate of 4 different databases (KEGG, Reactome, HMDB, and WikiPathways), which includes pathways and annotations for genes and metabolites. The database also includes chemical properties and chemical class annotations from HMDB, ChEBI and LIPIDMAPS.  

#### Note that the files here are NOT RaMP but merely the scripts used to create RaMP. ####

These scripts are used to build the RaMP database from scratch. 
If you are looking for the RaMP mysql dump, you will find it here: https://github.com/ncats/RaMP-DB/tree/master/inst/extdata
If you are looking to access RaMP through our shiny app, look here: https://github.com/ncats/RaMP-DB/

### The Basics ###
Here is the overall workflow for getting the mySQL database up and running:

  1. git clone this repository
  2. Run main.py successfully (purpose is to create intermediate data files in the misc/output folder for each data type.) Reading resource files from the various data sources is memory intensive. 
  3. Run EntityBuilder.py to merge data to build entity relationships. This produces a file for each table in the RaMP DB Schema in the /misc/sql/loader.
  4. Create an empty database in mySQL (e.g. called ramp) 
  5. Import the sql files into the database using src/util/rampDBBulkLoader. See code section for details: [DB Loading Code](https://github.com/ncats/RaMP-BackEnd/blob/9e0ab9c719f3a690272fc7a0ae669b6f11d74b7a/src/util/rampDBBulkLoader.py#L393)
  6. You can now query the database!

Keep scrolling down for the details...

### Environment set-up ###

#### Python ####
RaMP was built using using python 3.8, which must be downloaded here: https://www.python.org/downloads/
Make certain that you have python 3.8 or later.

#### Integrated Development Environment ####
Using an IDE makes writing code much simpler. 

One option is to use pyDev for Eclipse for python:
  - download Eclipse: https://www.eclipse.org/downloads/
  - install the pyDev software into Eclipse: http://www.wikihow.com/Install-Pydev-on-Eclipse

#### Documentation ####
The package is being documented using Sphinx: https://samnicholls.net/2016/06/15/how-to-sphinx-readthedocs/ and http://autoapi.readthedocs.io/

Set up requires adding the Sphinx to your PYTHONPATH. Traditionally, this would mean adding a line to your bash_rc. However if you are using eclipse to run code it doesn't register the PYTHONPATH set by the bash_rc. You will need to add the PYTHONPATH to the project using the eclipse interface. Here are the steps:
  1. Import Sphinx into eclipse.
  2. Right click the package in the left-hand pane and select properties.
  3. Select "Pydev - PYTHONPATH" from the options.
  4. Click "add source" button and select that *src* directory of the project
  5. Click apply and ok buttons.

#### Building the database ####
*main.py*

The main script to build the database is 'main.py'. This script calls the other classes and necessary code to build the database. 

The first time 'main.py' is run, one must swith the "getDatabaseFiles" parameter to "True".  Setting this parameter to true will download the content of the databases (HMDB, Reactome, WikiPathways) on your computer.  By default, we have set "getDatabasefiles" to FALSE because you may need to gather the information from each database again if you add new functionalities to the scripts (in which case there is no need to download the original files again). Note that the Reactome resource updates monthly and changes the name of their files to add the date.  The src/parse/reactomeData.py file's getDatabaseFile method contains the variable to set.

The output of running 'main.py' are intermediate files that are created in data-source specific files in misc/output. 

The next step is building entities and final files for database loading. The files in misc/output are input into an entity harmonization process that enforces entity curation patches, aggregates duplicate gene and compound entities from different data soruces and builds entity relationships. This step is kicked off by running /src/util/EntityBuilder.py. This step will refer to a text file that contains manual curation results that captures known problems in entity mappings held in certain data sources. The compound curation list prevents incorrect mappings between compounds from being introduced to the RaMP database. 

*importing sql tables*

Import the sql files into the database using src/util/rampDBBulkLoader. See code section for details: [DB Loading Code](https://github.com/ncats/RaMP-BackEnd/blob/9e0ab9c719f3a690272fc7a0ae669b6f11d74b7a/src/util/rampDBBulkLoader.py#L393)

[Additional information coming soon on bringing up mysql and initial creation of the database schema]

You may also need to download mySQL: https://www.mysql.com/downloads/. 



#### Overview of folders in this repo ####
The repo contains the following folders

    **src**: contains classes used by main.py (most of the code)

    **main**: contains main.py -- uses the classes in src (code)

    **docs**: documentation/user manual

    **tests**: unit testing

    **misc/data**: folders for files from hmdb, wikipathways, kegg, and reactome stored here

    **misc/sql**: .sql files output from running main.py. These .sql files make up the mySQL RaMP database

    **misc/output**: any output file that is not a .sql file -- some functions have other outputs, such as venndiagrams or lists of converted genes (can help with debugging)

    **misc/queryMySQL**: place to keep files that query the mySQL database once it already exists 

*Note: when running main.py, the data, output, and sql folders contain lots of data which is not included with this repo due to size restrictions.*

#### Updates ####
The first update to the RaMP (01/24/2019):

1. The major change is to adapt to the recent HMDB updates which changed the hierarchy of Ontology.
2. Exclusion of KEGG database.
3. Handled few compound mismatch that had occurred in the old version.

### Contact ###
* John Braisted (john.braisted@nih.gov)
* Ewy Mathé (ewy.mathe@nih.gov)

Elizabeth Baskin, Senyang Hu, and Bofei Zhang were also involved in developing this code.

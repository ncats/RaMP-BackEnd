# RaMP Creation Package #

This repository contains the source code that was used to create RaMP - a Relation database of Metabolic Pathways (https://www.ncbi.nlm.nih.gov/pubmed/29470400).  RaMP is a conglomerate of 4 different databases (Kegg, Reactome, HMDB, and WikiPathways), which includes pathways and annotations for genes and metabolites.  

#### Note that the files here are NOT RaMP but merely the scripts used to create RaMP. ####
If you are looking for the RaMP mysql dump, you will find it here: https://github.com/ncats/RaMP-DB/tree/master/inst/extdata
If you are looking to access RaMP through our shiny app, look here: https://github.com/ncats/RaMP-DB/

### The Basics ###
Here is the overall workflow for getting the mySQL database up and running:

  1. git clone this repository
  2. Run main.py successfully (purpose is to create sql files in misc/sql folder)
  3. Create an empty database in mySQL (e.g. called ramp) 
  4. Import the sql files into the database using Program.cs found in mySQLBulkCsharp folder
  5. You can now query the database!

Keep scrolling down for the details...

### Environment set-up ###

#### Python ####
RaMP was built using using python 3.5, which must be downloaded here: https://www.python.org/downloads/

If are using MacOS, python is already installed but it is python 2.7. You will need to install python 3.5 in addition to the 2.7 version already installed. You need to run both simultaneously because MacOS relies on python 2.7 for essential functions.

#### Integrated Development Environment ####
Using an IDE makes writing code much simpler. If you are familiar with R, Rstudio is the IDE for R. 

Another option is to use pyDev for Eclipse for python:
  - download Eclipse: https://www.eclipse.org/downloads/
  - install the pyDev software into Eclipse: http://www.wikihow.com/Install-Pydev-on-Eclipse

Yet another option is to use PyCharm. 

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

The first time 'main.py' is run, one must swith the "getDatabaseFiles" parameter to "True".  Setting this parameter to true will download the content of the databases (HMDB, KEGG, Reactome, WikiPathways) on your computer.  By default, we have set "getDatabasefiles" to FALSE because you may need to gather the information from each database again if you add new functionalities to the scripts (in which case there is no need to download the original files again).  

The output of running 'main.py' are sql files that then need to be imported into one database.

*importing sql tables*

A C# script, mySQLBulkCsharp, imports the sql files in bulk.  The "Program.cs" file is the main script. 

You may need to download the Visual Studio C# IDE: https://www.visualstudio.com/vs/
(may also need to download .net/mono...)

You may also need to download mySQL: https://www.mysql.com/downloads/. 

*plots*

Most of the code is written in Python but two scripts are written in R (src/threeVenn.R and fourVenn.R). These scripts produce plots that show the overlapping genes and metabolites between all databases in RaMP. You must have R installed on your computer, as well as these R packages: VennDiagram and grDevices.


#### Overview of folders in this repo ####
The repo contains the following folders

    **src**: contains classes used by main.py (most of the code)

    **main**: contains main.py -- uses the classes in src (code)

    **docs**: documentation/user manual

    **tests**: unit testing

    **mySQLBulkCsharp**: contains the C# script that will import the .sql files found in misc/sql into a mySQL database

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

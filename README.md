# RaMP Creation Package #


### What is this repository for? ###

This repository contains a package whose classes and functions work together to create RaMP (Relation database of Metabolic Pathways). RaMP is a conglomerate of 4 different databases (Kegg, Reactome, HMDB, and wikipathways) -- it pulls together the information scattered across them -- mainly compounds, genes, and pathways. Note that the files here are NOT RaMP but merely the scripts used to create RaMP. 

So...

The repository "MatheLabRaMP" = RaMP creation package: a python package for creating .sql files 

RaMP = the final mySQL databases made up of .sql files 

### How do I get set up? ###

If you are planning to help contribute to RaMP here is some information on how to get set up.

This RaMP creation package is using python 3.5, which must be downloaded here: https://www.python.org/downloads/

Note that if are using MacOS python is already installed but it is python 2.7. You will need to install python 3.5 in addition to the 2.7 version already installed. You need to run both simultaneously because MacOS relies on python 2.7 for essential functions.

Using an IDE makes writing code much simpler. If you are familiar with R, Rstudio is the IDE for R. I currently use pyDev for Eclipse for python.  

You must first download Eclipse: https://www.eclipse.org/downloads/

Then install the pyDev software into Eclipse: http://www.wikihow.com/Install-Pydev-on-Eclipse

The IDE PyCharm may be more popular/prettier/more stream-lined than pyDev. I just went with Eclipse because it is what is used in OSU computer science classes and it was more familiar to me. It doesn't really matter what IDE you use. 

The package is being documented using sphinx: https://samnicholls.net/2016/06/15/how-to-sphinx-readthedocs/ and http://autoapi.readthedocs.io/

Set up requires adding the package to your PYTHONPATH. Traditionally, this would mean adding a line to your bash_rc. However if you are using eclipse to run code it doesn't register the PYTHONPATH set by the bash_rc. You will need to add the PYTHONPATH to the project using the eclipse interface. Here are the steps:

1. Import the package into eclipse.
2. Right click the package in the left-hand pane and select properties.
3. Select "Pydev - PYTHONPATH" from the options.
4. Click "add source" button and select that *src* directory of the project
5. Click apply and ok buttons.

The script to build the database is main.py. When run, it calls the other classes and all the classes work together to build the database. The script takes ~10 minutes to run. If you want to run the script on the supercomputer how I do it is below: https://www.osc.edu/resources/technical_support/supercomputers/oakley

Here is the script I use to build the database on oakley:

```#PBS -l walltime=3:00:00

#PBS -l nodes=2:ppn=12

module load python/3.4.2

python /users/PAS1143/osu0084/mathelabramp/src/main.py
```

Save this in a file called "RaMP.sh".

Then to run the script:

```qsub RaMP.sh```

If you are running the script on eclipse, open the main.py file and click the "play" triangle icon. 

This package is mostly created using the python programming but there are two scipts in the src (threeVenn.R and fourVenn.R) that are written in the R programming language. R is known as the superior language for plotting graphs, although there are plotting modules available in python as well. The reason I moved to R is because there are no packages in python available to plot venn diagrams with 4 circles/entities (in our case, the four databases). 

Therefore, you must have R installed on your computer, as well as these R packages: VennDiagram and grDevices

In order to gather information from the databases (aka run the main.py script), you must download the content of the databases onto your computer. This process has been automated via a getDatabaseFiles function in each of the four database classes and running these four functions is the first step in the main.py script. However, you must toggle the getDatabaseFiles functions to "on" by switching the "getDatabaseFiles" parameter to "True." It is default false because although you only need to get the databases files once, you may need to run the information gathering portions of the package many times as you add new functionalities to the script and debug and it is a waste of time to re-download the files. 

After the the sql files are generated they need to be imported into the databases. We have a large number of lines in our files so a bulk import is more efficient. In order to do a bulk import you need a C# script.

Download the Visual Studio C# IDE: https://www.visualstudio.com/vs/
(may also need to download .net/mono...)

The C# script can be found in the mySQLBulkCsharp of this package. The "Program.cs" file is the main script. 

However, before running the script you will need to download mySQL database software: https://www.mysql.com/downloads/ and create the database "mathelabramp". 

Here is an overview of the folders in RaMP and what they contain:

**src**: contains classes used by main.py (most of the code)

**main**: contains main.py -- uses the classes in src (code)

**docs**: documentation/user manual

**tests**: unit testing

**mySQLBulkCsharp**: contains the C# script that will import the .sql files found in misc/sql into a mySQL database

**misc/data**: folders for files from hmdb, wikipathways, kegg, and reactome stored here

**misc/sql**: .sql files output from this package -- .sql files make up the mySQL RaMP database

**misc/output**: any output file that is not a .sql file -- some functions have other outputs, such as venndiagrams
or lists of converted genes (can help with debugging)

**misc/queryMySQL**: place to keep files that query the mySQL database once it already exists 

Note: since *data*,*output*,*sql* folder has large files, so they are not coming with this repo. Please expect a longer time for downloading data from beginning.

Here is the overall workflow for getting the mySQL database up and running:

1. git clone this repository
2. Run main.py successfully (purpose is to create sql files in misc/sql folder)
3. Create an empty database in mySQL called mathelabramp 
4. Import the sql files into the database using Program.cs found in mySQLBulkCsharp folder
5. You can now query the database!



### Who do I talk to? ###
* Repo owner: Bofei Zhang (zhang.5675@osu.edu)
* Repo owner: Elizabeth Baskin (Elizabeth.Baskin@osumc.edu)
* PI: Ewy Mathe (Ewy.Mathe@osumc.edu)

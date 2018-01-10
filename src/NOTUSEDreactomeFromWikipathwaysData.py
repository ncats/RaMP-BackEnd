import urllib.request
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree
import os

class NOTUSEDreactomeDataFromWikipathways():
    
    '''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!This is class no longer used since many of the reactome pathways could not be obtained via wikipathways!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    The  reactomeData class actually gets information from reactome via wikipathways.
        
    Although reactome IS a separate database, the wikipathways database has worked to curate the reactome database and incorporate it 
    into the wikipathways database. The reason it is preferable to get the reactome database via wikipathways is that it is better 
    formatted and easier to acquire the information we want from there.
    
    All the information  is acquired by downloading a number of xml files called gpml files. Each file represents 
    a biochemical pathway composed of metabolite and gene interactions. As a result, metabolites and genes can  be linked to pathways.
    
    The reactome files can be downloaded here: http://www.wikipathways.org/index.php/Download_Pathways​
    
    The files have already been downloaded and placed in this package. The path to the files has been hardcoded into the function calls.
    
    This class only has two functions:
    
        - 1) getEverything
        - 2) writeToFiles
        
    Due to the nature of the structure of this database (a series of xml files) it makes more sense to get all the information at one time 
    from one function call and then, as in the other classes, write everything to sqls for the RaMP database. 
    
    
    
    '''
    
    def __init__(self):
    
        #key: ID for pathway, Value: pathway name
        #Right now this ID is generated in order of processing the file because a wikipathways ID hasn't been found
        self.pathwayDictionary = dict()
        
        #hsaID for pathway, value: pathway category (all will be "NA" for reactome)
        self.pathwayCategory = dict()
        
        #key: metabolite, value: list of pathways
        self.metabolitesWithPathwaysDictionary = dict()
        
        #empty for reactome
        self.metabolitesWithSynonymsDictionary = dict()
        
        #key: metabolite, value: metabolite mapping
        self.metaboliteIDDictionary = dict()
        
        #key: pathwayID, value: list of genes
        self.pathwaysWithGenesDictionary = dict()
        
        #key: gene, value: gene mapping
        self.geneInfoDictionary = dict()

        #key: pathwayId, value: list of metabolites
        self.pathwayWithMetabolitesDictionary = dict()
        
        #stays empty for this class
        self.pathwayOntology = dict()
        
        #stays empty for this class
        self.biofluidLocation = dict()
        
        #stays empty for this class
        self.biofluid = dict()
        
        #stays empty for this class
        self.cellularLocation = dict()
        
        #stays empty for this class
        self.cellular = dict()
        
       

        
    def getEverything(self):
        
      ''' This function gets all the necessary information from the gpml files and places it into dictionaries. 
          Unlike in other classes, only one function is required to get all the necessary information.
        
      '''
        
        
      pathreactome = "./mathelabramp/data/wikipathways/wikipathways_Homo_sapiens_Curation-Reactome_Approved__gpml/"
      
      
      
      dictionaryOfFiles = dict()
      
      for filename in os.listdir(pathreactome):
          fullpath = pathreactome + filename
          dictionaryOfFiles[fullpath] = "reactome"
      
      
      for filename in dictionaryOfFiles:
          
          
          print("current file: " + str(filename))
          
          tree = ET.parse(filename)
          
              
          root = tree.getroot()
          
          metaboliteMapping = {
                    "chebi_id": "NA", 
                    "drugbank_id": "NA", 
                    "drugbank_metabolite_id": "NA", 
                    "phenol_explorer_compound_id": "NA", 
                    "phenol_explorer_metabolite_id": "NA", 
                    "foodb_id": "NA", 
                    "knapsack_id": "NA", 
                    "chemspider_id": "NA",
                    "kegg_id": "NA",
                    "biocyc_id": "NA",
                    "bigg_id": "NA",
                    "wikipidia": "NA",
                    "nugowiki": "NA",
                    "metagene": "NA",
                    "metlin_id": "NA",
                    "pubchem_compound_id": "NA",
                    "het_id": "NA",
                    "hmdb_id": "NA",
                    "CAS": "NA"}
          
          
          geneMapping = {
                         "kegg": "NA",
                         "common_name": "NA",
                         "Ensembl": "NA", 
                         "HGNC": "NA", 
                         "HPRD": "NA", 
                         "NCBI-GeneID": "NA", 
                         "NCBI-ProteinID": "NA", 
                         "OMIM": "NA", 
                         "UniProt": "NA", 
                         "Vega": "NA", 
                         "miRBase": "NA", 
                         "HMDB_protien_accession": "NA",
                         "Entrez" : "NA",
                         "Enzyme Nomenclature": "NA"}
          
          
          
          
          
       
          
         
          

          
          
          listOfGenes = []
          listOfMetabolites = []
          
          pathwayname = root.get("Name")
          
          
          for child in root:
              childtag = child.tag.replace("{http://pathvisio.org/GPML/2013a}", "")
              
              
              if childtag == "Attribute":
                  reactomeID = child.get("Value")
                  pathwayID = reactomeID
                  
                  #this keeps track of the current pathway. 
                  #this is important later for placing pathways into the metaboliteswithpathways dictionary 
                  currentpathway = pathwayID
                  self.pathwayDictionary[pathwayID] = pathwayname
                  self.pathwayCategory[pathwayID] = "NA"
              
              if childtag == "DataNode":
                  metaboliteorgene = child.get("TextLabel")
                  type = child.get("Type")
                  for child2 in child:
                      database = child2.get("Database")
                      databaseID = child2.get("ID")
                      childtag = child2.tag
                      childtag = child2.tag.replace("{http://pathvisio.org/GPML/2013a}", "")
                      if childtag == "Xref":
                          
                          #"proteins" are the only gene-like things in reactome files
                          if type == "Protein":
                              geneMapping["common_name"] = metaboliteorgene
                              
                              if databaseID is not "" and database == "Uniprot-TrEMBL": 
                                  geneMapping["UniProt"] = databaseID
            
                                  if databaseID not in listOfGenes:
                                      listOfGenes.append(databaseID)
                             
                                  
                                  if databaseID not in listOfGenes:
                                      listOfGenes.append(databaseID)
                              
                              self.geneInfoDictionary[databaseID] = geneMapping
                          
                          if type == "Metabolite":

                              if database == "ChEBI":
                                  #remove prefix
                                  databaseID = databaseID.replace("CHEBI:", "")
                                  #kegg makes a list of chebi ids since there are more than one. For consistency, and so both databases
                                  #can use the ID conversion class, we will also make a "list" of chebi ids here, although I have not seen
                                  #an instance where there is more than one chebi in this database
                                  metaboliteMapping["chebi_id"] = [databaseID]
                              
                              
                              if databaseID not in listOfMetabolites:
                                  listOfMetabolites.append(databaseID)
                            
                              self.metaboliteIDDictionary[databaseID] = metaboliteMapping
                              
                              
                              
                              if databaseID not in self.metabolitesWithPathwaysDictionary:
                                  listOfPathways = []
                                  listOfPathways.append(pathwayID)
                                  self.metabolitesWithPathwaysDictionary[databaseID] = listOfPathways
                                  currentpathway = pathwayID
                                  
                              #this is what should happen if the metabolite has already been seen but we are on a NEW pathway file. If it is the same 
                              #pathway file it will not be recorded again.     
                              elif currentpathway != pathwayID:
                                  
                                  value = self.metabolitesWithPathwaysDictionary[databaseID]
                                  print(databaseID)
                                  value = value.append(pathwayID)
                                  self.metabolitesWithPathwaysDictionary[databaseID] = value
                                  
            
          self.pathwaysWithGenesDictionary[pathwayID] = listOfGenes
          self.pathwayWithMetabolitesDictionary[pathwayID] = listOfMetabolites
          
          
                            
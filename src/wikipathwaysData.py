import urllib.request as RE
import libchebipy
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree
import os
import zipfile
import time
from MetabolomicsData import MetabolomicsData



class wikipathwaysData(MetabolomicsData):
    
    '''
    The wikipathwaysData class actually gets information from wikipathways via number of xml files called gpml files. Each file represents 
    a biochemical pathway composed of metabolite and gene interactions. As a result, metabolites and genes can  be linked to pathways.
    
    The wikipathways files can be downloaded here: http://www.wikipathways.org/index.php/Download_Pathways​
    The reactome files can be downloaded here: http://www.wikipathways.org/index.php/Download_Pathways​
    
    The files have already been downloaded and placed in this package. The path to the files has been hardcoded into the function calls.
    
    To fill the data of this class, just call:
    self.getEverything()
        This function will check if the file are downloaded.
        Then, parse all files to fill dict(). Finally query chebi API
        to get all common names for metabolites.
    
        
    Due to the nature of the structure of this database (a series of xml files) it makes more sense to get all the information at one time 
    from one function call and then, as in the other classes, write everything to sqls for the RaMP database. 
    
    Update @ 5/9/2018
    
    This class is deprecated since we are not using GPML file to import data.
    
    '''
    
    def __init__(self):
        
        super().__init__()
        # key: ID for metabolites Value: Common Name (the only name in this database)
        self.metaboliteCommonName = dict()
        #key: ID for pathway, Value: pathway name
        self.pathwayDictionary = dict()
        
        #key: ID for pathway, category: various categories such as cellular process, metabolic process 
        self.pathwayCategory = dict()
        
        #key: gene, value: gene mapping
        self.geneInfoDictionary = dict()
        
        #key: metabolite, value: metabolite mapping
        self.metaboliteIDDictionary = dict()
        
        #key: pathwayID, value: list of genes
        self.pathwaysWithGenesDictionary = dict()
        
        #key: pathwayId, value: list of metabolites
        self.pathwayWithMetabolitesDictionary = dict()
        
        #empty for reactome
        self.metabolitesWithSynonymsDictionary = dict()
        
        #only not empty when a catalyzed class exists 
        #empty
        self.metabolitesLinkedToGenes = dict()
        
        #key: metabolite, value: list of pathways
        self.metabolitesWithPathwaysDictionary = dict()
        
        #key: current pathway, value: list of pathways linked to this pathway
        self.pathwayOntology = dict()
        
        #stays empty for this class
        self.biofluidLocation = dict()
        
        #stays empty for this class
        self.biofluid = dict()
        
        #stays empty for this class
        self.cellularLocation = dict()
        
        #stays empty for this class
        self.cellular = dict()
        
        #stays empty
        self.exoEndoDictionary = dict()
        self.exoEndo =dict()
        #tissue location stay empty
        self.tissue = dict()
        self.tissueLocation = dict()
        # show how id distributes
        # key database value: frequency
        '''
        self.metaboliteFromWhichDB = {
                                      "kegg_id":0,
                                      "hmdb_id":0,
                                      "CAS":0,
                                      "chebi_id":0,
                                      "pubchem_compound_id":0,
                                      "chemspider_id":0
                                      }

        self.setOfType = set()
        self.idSetFromGeneProducts = dict()
        self.idSetFromProtein = dict()
        '''
    
    def IDconversion(self,this_id):
        '''
        This functions replace some special character for the ID text
        Also add two zeroes to HMDB ids
        '''
        to_convert = {'\n':'',
                      '\t':'',
                      ' ':'',
                      'HMDB':'HMDB00',
                      'CHEBI:':'chebi:'}  
        for key,value in to_convert.items():
            this_id = this_id.replace(key,value)
        return this_id 
    
    def prepend(self,database,id):
        '''
        Prepend some prefix to avoid id mapping problem when creating ramp id
        e.g. pure number id are overlapped but from different source 
        '''
        prefix ={
                 'pubchem_compound_id':'pubchem:',
                 'Entrez':'entrez:',
                 'chemspider_id':'chemsipder:'
                 }
        if database in prefix:
            pre = prefix[database]
            id = pre + id
        return id
    def getDatabaseFiles(self):
        
        '''
        This function gets the files that make up wikipathways and place them to wikipathways folder.
        
        The data file is stored in the url: data.wikipathways.org/ in the order of published date
        
        '''
        
        path = "../misc/data/wikipathways/"
        url = 'http://data.wikipathways.org/20180210/gpml/wikipathways-20180210-gpml-Homo_sapiens.zip'
        file = "wikipathways-20180210-gpml-Homo_sapiens.zip"
        
        self.check_path(path)
        existed = os.listdir(path)
        print(existed)
        if file not in existed:
            self.download_files(url, 
                                path+file)
            
            with zipfile.ZipFile(path + file,"r") as zip_ref:
                zip_ref.extractall(path)
        else:
            print("Already Downloaded ....")
        
        
    def getEverything(self):
        
        ''' This function gets all the necessary information from the gpml files and places it into dictionaries. 
        Unlike in other classes, only one function is required to get all the necessary information.
        The update date needs to be hardcoded to the url.
        Please make sure going to the website to see which one is the latest version for 
        wikipathways data.
        '''
        # make sure the database is gotten.
        self.getDatabaseFiles()
        print("Call wikipathways getEverything......")  
        pathwikipathways = "../misc/data/wikipathways/"
      
        dictionaryOfFiles = dict()
        possible_attrType = set()
        databaseType = dict()
        '''
        All possible collection of analytes type with their possible database sources are listed here
        {'Protein': {'', 'Entrez Gene', 'ChEMBL compound', 'PubChem-compound', 'miRBase mature sequence', 'Uniprot-SwissProt', 'Ensembl', 'NCBI Protein', 'Uniprot-TrEMBL', 'Wikipedia', 'Ensembl Human', 'Enzyme Nomenclature', 'Wikidata', 'InterPro'}, 
        None: {'', 'Ensembl', 'Wikipedia', 'Uniprot-TrEMBL', 'OMIM', 'Wikidata', 'ChEBI', 'NanoParticle Ontology'}, 
        'Rna': {'', 'Entrez Gene', 'HGNC', 'miRBase Sequence', 'miRBase mature sequence', 'Ensembl'}, 
        'Pathway': {'', 'KEGG Pathway', 'Reactome', 'WikiPathways', 'Wikipedia'}, 
        'Metabolite': {'', 'KEGG Drug', 'Chemspider', 'HMDB', 'ChEMBL compound', 'PubChem-compound', 'LIPID MAPS', 'DrugBank', 'PubChem-substance', 'Wikidata', 'ChEBI', 'CAS', 'KEGG Compound'}, 
        'GeneProduct': {'', 'Entrez Gene', 'HGNC', 'RefSeq', 'miRBase Sequence', 'miRBase mature sequence', 'Kegg ortholog', 'GeneOntology', 'Enzyme Nomenclature', 'Wikidata', 'Pfam', 'undefined', 'Ensembl Human', 'pato', 'Uniprot-SwissProt', 'Wikipedia', 'Uniprot-TrEMBL', 'EcoGene', 'KEGG Genes', 'Ensembl', 'NCBI Protein', 'Other'}, 
        'Complex': {'', 'Wikipedia', 'Reactome', 'Uniprot-SwissProt'}}
        '''
        # all data types are collected are shown in this dict
        target_analytes = {
            'Protein': {'Entrez Gene','Uniprot-SwissProt', 
                        'Ensembl','Uniprot-TrEMBL','Ensembl Human', 
                        'Enzyme Nomenclature', 
                        'Wikidata'},
            'GeneProduct': {'Entrez Gene','Enzyme Nomenclature', 'Wikidata', 
                            'Ensembl Human', 'Uniprot-SwissProt', 'Uniprot-TrEMBL',
                            'KEGG Genes','Ensembl'},
            'Metabolite': {'Chemspider','HMDB','PubChem-compound','LIPID MAPS', 
                            'Wikidata', 'ChEBI', 'CAS', 'KEGG Compound'}
            
            }
        for filename in os.listdir(pathwikipathways):
            if ".zip" in filename:
                continue  # the original downloaded file is also in same path
            # and not parsed to the later process
            fullpath = pathwikipathways + filename
            dictionaryOfFiles[fullpath] = "wikipathways"
          
        print("done ... for input files name.")  
        for filename in dictionaryOfFiles:
          
            tree = ET.parse(filename)
          
              
            root = tree.getroot()

           
            #how to get the pathwayID from filename
            last = filename.rfind('_')
            first = filename[:last].rfind('_')
            pathwayID = filename[first+1:last]
          
          
            #this keeps track of the current pathway. 
            #this is important later for placing pathways into the metaboliteswithpathways dictionary 
            currentpathway = pathwayID
          

            
          
            listOfGenes = []
            listOfMetabolites = []
            listOfPathwaysForOntology = []
          
            pathwayname = root.get("Name")
            # link all pathways dictionary
            self.pathwayDictionary[pathwayID] = pathwayname
            self.pathwayCategory[pathwayID] = 'NA'
            self.pathwaysWithGenesDictionary[pathwayID] = []
            print('Parsing pathway {}:{}'.format(pathwayID,pathwayname))
            prefix = "{http://pathvisio.org/GPML/2013a}"
            dataNode = root.findall(prefix + 'DataNode')
            for dn in dataNode:
                textLabel = dn.get('TextLabel')
                textLabel = self.IDconversion(textLabel)
                type = dn.get('Type')
                xref = dn.find(prefix + 'Xref')
                database = xref.get('Database')
                databaseID = xref.get('ID')
                databaseID = self.IDconversion(databaseID)
                if type in ['Protein','GeneProduct']:
                    if database in target_analytes[type]:
                        if databaseID is None or databaseID is '':
                            print('No id ?????')
                            print(databaseID)
                            continue
                            time.sleep(1)
    
                        geneMapping = {"kegg": "NA",
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
                                        "HMDB_protein_accession": "NA",
                                        "Entrez" : "NA",
                                        "Enzyme Nomenclature": "NA",
                                        'WikiData':"NA"}
                        convert = {'Entrez Gene':'Entrez',
                                   'Uniprot-SwissProt':'UniProt', 
                                   'Ensembl':'Ensembl',
                                   'Uniprot-TrEMBL':'UniProt',
                                   'Ensembl Human':'Ensembl', 
                                   'Enzyme Nomenclature':'Enzyme Nomenclature', 
                                   'Wikidata':'WikiData',
                                   'KEGG Genes':'kegg'}
                        key_for_mapping = convert[database]
                        # wiki store kegg pathway as a gene product # might be wrong for ramp
                        if 'hsa' in databaseID:
                            print('Analyte {} is {} from database {} with id {}'.\
                                  format(textLabel,type,database,databaseID))
                            continue
                        if databaseID not in self.geneInfoDictionary and databaseID is not '':
                            databaseID = self.prepend(key_for_mapping, databaseID)
                            geneMapping[key_for_mapping] = databaseID
                            geneMapping['commonName'] = textLabel
                            self.geneInfoDictionary[databaseID] = geneMapping
                            self.pathwaysWithGenesDictionary[pathwayID].append(databaseID)
                        
                elif type in ['Metabolite']:
                    if database in target_analytes[type]:
                        if databaseID is None or databaseID is '':
                            print('No id ?????')
                            print(databaseID)
                            time.sleep(1)
                            continue
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
                                  "wikipedia": "NA",
                                  "nugowiki": "NA",
                                  "metagene": "NA",
                                  "metlin_id": "NA",
                                  "pubchem_compound_id": "NA",
                                  "het_id": "NA",
                                  "hmdb_id": "NA",
                                  "CAS": "NA",
                                  'LIPIDMAPS':'NA',
                                  'WikiData':'NA'}
                        convert = {'Chemspider':'chemspider_id'
                                   ,'HMDB':'hmdb_id',
                                   'PubChem-compound':'pubchem_compound_id',
                                   'LIPID MAPS':'LIPIDMAPS', 
                                   'Wikidata':'WikiData', 
                                   'ChEBI':'chebi_id', 
                                   'CAS':'CAS', 
                                   'KEGG Compound':'kegg_id'}
                        key_for_mapping = convert[database]
                        databaseID = self.prepend(key_for_mapping, databaseID)
                        metaboliteMapping[key_for_mapping] = [databaseID]
                        if 'hsa' in databaseID:
                            print('Analyte {} is {} from database {} with id {}'.\
                                  format(textLabel,type,database,databaseID))
                            continue
                        if databaseID not in self.metaboliteIDDictionary:
                            self.metaboliteIDDictionary[databaseID] = metaboliteMapping
                            self.metaboliteCommonName[databaseID] = textLabel
                            self.metabolitesWithSynonymsDictionary[databaseID] = [textLabel]
                            '''
                        else:
                            print('#### Repeated metabolites {}:{}'.format(textLabel,databaseID))
                            '''
                        if databaseID not in self.metabolitesWithPathwaysDictionary:
                            self.metabolitesWithPathwaysDictionary[databaseID] = [pathwayID]
                        else:
                            if pathwayID not in self.metabolitesWithPathwaysDictionary[databaseID]:
                                self.metabolitesWithPathwaysDictionary[databaseID].append(pathwayID)
                            #time.sleep(0.1)
                        
                        '''
                        else:
                            print('Analyte {} is {} from database {} with id {}'.\
                                  format(textLabel,type,database,databaseID))
                        '''
        self.getCommonNameForChebi()
        print('Total genes {};Total metabolites {}'\
              .format(len(self.geneInfoDictionary),len(self.metaboliteIDDictionary)))
        print('Relations path-gene/path-meta is {} and {}'\
              .format(len(self.pathwaysWithGenesDictionary),
                      len(self.metabolitesWithPathwaysDictionary)))
        print('Common name for metabolites: {}'.format(len(self.metaboliteCommonName)))
        
                            if Attributetype == "Protein" or Attributetype == 'GeneProduct':
                                geneMapping = {"kegg": "NA",
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
                                             "HMDB_protein_accession": "NA",
                                             "Entrez" : "NA",
                                             "Enzyme Nomenclature": "NA"}
                                
                                geneMapping["common_name"] = metaboliteorgene
                              
                                if databaseID is not "" and database == "Entrez Gene":
                                    databaseID = 'entrez:' + databaseID.replace(' ','') 
                                    geneMapping["Entrez"] = [databaseID]
                                    geneMapping['kegg'] = [databaseID.replace('entrez:','')]
                                if databaseID not in listOfGenes:
                                    listOfGenes.append(databaseID)
                                  
                                  
                                if databaseID is not "" and database == "Enzyme Nomenclature": 
                                    geneMapping["Enzyme Nomenclature"] = databaseID
                                    if databaseID not in listOfGenes:
                                        listOfGenes.append(databaseID)
                                  
                                  
                                if databaseID is not "" and database == "Ensembl": 
                                    geneMapping["Ensembl"] = [databaseID]
                                    if databaseID not in listOfGenes:
                                        listOfGenes.append(databaseID)
                                  
                                  
                                if databaseID is not "" and database == "Uniprot-TrEMBL": 
                                    geneMapping["UniProt"] = [databaseID]
                                    if databaseID not in listOfGenes:
                                        listOfGenes.append(databaseID)
                                
                                self.geneInfoDictionary[databaseID] = geneMapping
                                
                            if Attributetype == "GeneProduct":
                                
                                geneMapping = {"kegg": "NA",
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
                                             "HMDB_protein_accession": "NA",
                                             "Entrez" : "NA",
                                             "Enzyme Nomenclature": "NA"}
                                
                                geneMapping["common_name"] = [metaboliteorgene]
                              
                                if databaseID is not "" and database == "Entrez Gene": 
                                    databaseID = 'entrez:' + databaseID.replace(' ','')
                                    geneMapping["Entrez"] = [databaseID]
                                    
                                if databaseID not in listOfGenes:
                                    listOfGenes.append(databaseID)
                                  
                                  
                                if databaseID is not "" and database == "Enzyme Nomenclature": 
                                    geneMapping["Enzyme Nomenclature"] = databaseID
                                    if databaseID not in listOfGenes:
                                        listOfGenes.append(databaseID)
                                  
                                  
                                if databaseID is not "" and database == "Ensembl": 
                                    geneMapping["Ensembl"] = [databaseID]
                                    if databaseID not in listOfGenes:
                                        listOfGenes.append(databaseID)
                                  
                                  
                                if databaseID is not "" and database == "Uniprot-TrEMBL": 
                                    geneMapping["UniProt"] = [databaseID]
                                    if databaseID not in listOfGenes:
                                        listOfGenes.append(databaseID)
                                self.geneInfoDictionary[databaseID] = geneMapping
                                      
                              
                                  
                              
                             
                              
                              
                          
                            if Attributetype == "Metabolite":
                                metaboliteorgene = 'NA' 
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
                                  "CAS": "NA",
                                  'LIPIDMAPS':'NA'}
                                # some of the common name has some special character follows '&'
                                
                                  
                                if database == "HMDB":
                                    
                                    # databaseID is the source ID 
                                    databaseID = databaseID.replace("HMDB","HMDB00")
                                    metaboliteMapping["hmdb_id"] = [databaseID]
                                    self.metaboliteIDDictionary[databaseID] = metaboliteMapping
                                    self.metaboliteCommonName[databaseID] = metaboliteorgene
                                    if databaseID not in listOfMetabolites:
                                        listOfMetabolites.append(databaseID)
                                        
                                
                                    self.metabolitesWithSynonymsDictionary[databaseID] = [metaboliteorgene]
                                  
                                  
                                if database == "CAS":
                                    metaboliteMapping["CAS"] = databaseID
                                    self.metaboliteIDDictionary[databaseID] = metaboliteMapping
                                  
                                    if databaseID not in listOfMetabolites:
                                        listOfMetabolites.append(databaseID)
                                        #self.metaboliteFromWhichDB["CAS"] +=1
                                    self.metabolitesWithSynonymsDictionary[databaseID] = [metaboliteorgene]
                                    self.metaboliteCommonName[databaseID] = metaboliteorgene

                                  
                                if database == "ChEBI":
                                    #remove prefix
                                    if 'CHEBI' in databaseID:
                                        databaseID = databaseID.replace("CHEBI:", "chebi:")
                                    else:
                                        databaseID = 'chebi:' +databaseID
            
                                    #kegg makes a list of chebi ids since there are more than one. For consistency, and so both databases
                                    #can use the ID conversion class, we will also make a list of chebi ids here
                                    metaboliteMapping["chebi_id"] = [databaseID]
                                    self.metaboliteIDDictionary[databaseID] = metaboliteMapping
                                    self.metaboliteCommonName[databaseID] = metaboliteorgene

                                    if databaseID not in listOfMetabolites:
                                        listOfMetabolites.append(databaseID)
                                        #self.metaboliteFromWhichDB["chebi_id"] += 1
                                      
                                    self.metabolitesWithSynonymsDictionary[databaseID] = [metaboliteorgene]
                                  
                                if database == "KEGG Compound":
                                    databaseID = databaseID.replace(' ','')
                                    metaboliteMapping["kegg_id"] = databaseID
                                    self.metaboliteIDDictionary[databaseID] = metaboliteMapping
                                    self.metaboliteCommonName[databaseID] = metaboliteorgene
                                    if databaseID not in listOfMetabolites:
                                        listOfMetabolites.append(databaseID)
                                        #self.metaboliteFromWhichDB["kegg_id"] += 1
                                    self.metabolitesWithSynonymsDictionary[databaseID] = [metaboliteorgene]
                                  
                                if database == "PubChem-compound":
                                    databaseID = 'pubchem:'+databaseID
                                    metaboliteMapping["pubchem_compound_id"] = databaseID
                                    self.metaboliteIDDictionary[databaseID] = metaboliteMapping
                                    self.metaboliteCommonName[databaseID] = metaboliteorgene

                                    if databaseID not in listOfMetabolites:
                                        listOfMetabolites.append(databaseID)
                                        #self.metaboliteFromWhichDB["pubchem_compound_id"] += 1
                                    self.metabolitesWithSynonymsDictionary[databaseID] = [metaboliteorgene]
                                  
                                if databaseID is not "" and database == "Chemspider": 
                                    databaseID = 'chemspider:' + databaseID
                                    metaboliteMapping["chemspider_id"] = databaseID
                                    self.metaboliteIDDictionary[databaseID] = metaboliteMapping
                                    self.metaboliteCommonName[databaseID] = metaboliteorgene

                                    if databaseID not in listOfMetabolites:
                                        listOfMetabolites.append(databaseID)
                                        #self.metaboliteFromWhichDB["chemspider_id"] += 1
                                      
                                    self.metabolitesWithSynonymsDictionary[databaseID] = [metaboliteorgene]
                                # Not collecting pubchem-substance id due to inability to differentiate it from
                                # pubchem compound id.
                                if databaseID not in self.metabolitesWithPathwaysDictionary:
                                    listOfPathways = []
                                    listOfPathways.append(pathwayID)
                                    if database in ["HMDB", "CAS", "ChEBI", "KEGG Compound", "PubChem-compound", "Chemspider", "PubChem-substance"]:
                                        if database == "HMDB" and len(databaseID) < 11:
                                            databaseID = databaseID.replace("HMDB","HMDB00")
                                        self.metabolitesWithPathwaysDictionary[databaseID] = listOfPathways
                                        currentpathway = pathwayID
                                  
                                #this is what should happen if the metabolite has already been seen but we are on a NEW pathway file. If it is the same 
                                #pathway file it will not be recorded again.     
                                elif currentpathway != pathwayID:
                                   
                                    value = self.metabolitesWithPathwaysDictionary[databaseID]
                                    value = value.append(pathwayID)
                                    if database in ["HMDB", "CAS", "ChEBI", "KEGG Compound", "PubChem-compound", "Chemspider", "PubChem-substance"]:
                                        self.metabolitesWithPathwaysDictionary[databaseID] = listOfPathways
                                    self.metabolitesWithPathwaysDictionary[databaseID] = value
                                  
                            if Attributetype == "Pathway":
                                if database == "WikiPathways":
                                    listOfPathwaysForOntology.append(databaseID) 
          
            self.pathwayOntology[pathwayID] = listOfPathwaysForOntology        
            self.pathwaysWithGenesDictionary[pathwayID] = listOfGenes
            self.pathwayWithMetabolitesDictionary[pathwayID] = listOfMetabolites
        self.getCommonNameForChebi()
     
    
    def getCommonNameForChebi(self):
        '''
        Get common name from chebi API,
        This function is called inside self.getEverything
        '''
        for key in self.metaboliteIDDictionary:
            value = self.metaboliteIDDictionary[key]
            chebi = value["chebi_id"]
            name = "NA"
            for each in chebi:
                try:
                    chebiToSearch = str(each)
                    chebiToSearch2 = libchebipy.ChebiEntity(chebiToSearch)
                    name = chebiToSearch2.get_name()
                    #print(chebiToSearch)
                    #print(name)
                    #time.sleep(3)
                except:
                    pass
            if key not in self.metabolitesWithSynonymsDictionary:
                self.metabolitesWithSynonymsDictionary[key] = [name]
            else:
                self.metabolitesWithSynonymsDictionary[key].append(name)     
          

                
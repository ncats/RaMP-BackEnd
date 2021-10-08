import urllib.request 
import urllib.error as ER
import libchebipy
import time
import os
from os import path
import xml.etree.ElementTree as ET
from parse.MetabolomicsData import MetabolomicsData
from multiprocessing.dummy import Pool

class reactomeData(MetabolomicsData):
    '''
    This class contains all dict stored and distributed from Reactome Official websites.
    
    '''
    def __init__(self, resConfig):
        super().__init__()
        
        self.resourceConfig = resConfig
        
        # key: Chebi ID from reactome Value: common Name by querying Chebi Database
        self.metaboliteCommonName = dict()
        #key: ID for pathway, Value: pathway name
        
        self.pathwayDictionary = dict()
        
        #hsaID for pathway, value: pathway category (all will be "NA" for reactome)
        self.pathwayCategory = dict()
        
        #key: metabolite, value: list of pathways
        self.metabolitesWithPathwaysDictionary = dict()
        
        #key: chebi, value, synonyms
        self.metabolitesWithSynonymsDictionary = dict()
        
        #key: metabolite, value: metabolite mapping
        self.metaboliteIDDictionary = dict()
        
        
        #only not empty when a catalyzed class exists 
        #empty
        self.metabolitesLinkedToGenes = dict()
        
        #key: pathwayID, value: list of genes
        self.pathwaysWithGenesDictionary = dict()
        
        #key: gene, value: gene mapping
        self.geneInfoDictionary = dict()
        
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
        
        #stays empty
        self.exoEndoDictionary = dict()
        self.exoEndo = dict()
        self.tissueLocation = dict()
        self.tissue = dict()
    
    def getEverything(self,writeToFile = False):
        '''
        This function runs all functions below to fill all dictionaries
        
        '''
        self.getDatabaseFiles()
        # print("Getting genes...")
        # self.getGenes()
        print("Getting metabolites...")
        self.getMetabolites()
        
        print("Getting common names...")
        self.getCommonNameForChebi()
        
        
        print("Getting genes ...")        
        self.getGenes("proteins")
        
        # separate query on genes is not needed.
        # get gene ids (NCBI GeneID) when getting gene symbol from Uniprot files.
        # self.getGenes("genes")
        
        print("Getting common names for genes1 ...")
        self.downloadCommonNameFromUniprot()
        print("Getting common names for genes 2...")
        self.getCommonNameFromUniprot()
        print("Getting common names for genes 3...")
        if writeToFile:
            self.write_myself_files('reactome')

        print("Done with ReactomeT ...")
    def getDatabaseFiles(self):
        
        '''This function gets the files that make up reactome and places them into the reactome folder. 

        '''
        proteinConfig = self.resourceConfig.getConfig("reactome_gene")
        metConfig = self.resourceConfig.getConfig("reactome_met")
        
        
        url_proteins,dir_proteins,file_proteins = (proteinConfig.sourceURL, proteinConfig.localDir, proteinConfig.extractFileName)
        
#         url_genes,dir_proteins,file_genes = ("http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt",
#                                                             "../misc/data/reactome/",
#                                                             "NCBI2Reactome_All_Levels.txt")
        
        url_metabolites, dir_metabolites, file_metabolites = (metConfig.sourceURL, metConfig.localDir, metConfig.extractFileName)
#        existed = os.listdir(dir_proteins)

#        if self.check_path(dir_proteins) :  
            
#             if file_metabolites not in existed or file_proteins not in existed or file_genes not in existed:                                                
#                 self.download_files(url_proteins,dir_proteins+file_proteins)
#                 self.download_files(url_genes,dir_proteins+file_genes)
#                 self.download_files(url_metabolites, dir_metabolites+file_metabolites)
#             else:
#                 print("Already downloaded ...")

#        else:
        existed = os.listdir(dir_proteins) 
        if file_metabolites not in existed or file_proteins not in existed:                                                
            self.download_files(url_proteins,dir_proteins+file_proteins)
            self.download_files(url_metabolites, dir_metabolites+file_metabolites)
        else:
            print("Already downloaded ...")
            
        time.sleep(1)
        
        
    def getGenes(self, proteinsOrGenes): 
            
        if proteinsOrGenes == "proteins":    
            reactomeFile = open("../misc/data/reactome/UniProt2Reactome_All_Levels.txt", encoding="utf-8")
        else:
            reactomeFile = open("../misc/data/reactome/NCBI2Reactome_All_Levels.txt", encoding="utf-8")
            
        
        for line in reactomeFile:
            
            splitline = line.split("\t")       
            
            if len(splitline) > 2:
                if "Homo sapiens" in splitline[5]:
                    gene = splitline[0]
                    # gene = 'uniprot:'+gene
                    #print("gene:", gene)
                    pathwayID = splitline[1]
                    #print("pathwayID:", pathwayID)
                    pathwayName = splitline[3]
                    #print("pathwayName:", pathwayName)
                    
                    mapping = { 'kegg': 'NA',
                                'common_name': 'NA',
                                'Ensembl': 'NA', 
                                'HGNC': 'NA', 
                                'HPRD': 'NA', 
                                'NCBI-GeneID': 'NA', 
                                'NCBI-ProteinID': 'NA', 
                                'OMIM': 'NA', 
                                'UniProt': ['uniprot:'+gene], 
                                'Vega': 'NA', 
                                'miRBase': 'NA', 
                                'HMDB_protein_accession': 'NA',
                                'Entrez': 'NA',
                                'Enzyme Nomenclature': 'NA'}

                    if pathwayID not in self.pathwaysWithGenesDictionary:
                        if proteinsOrGenes == "proteins":
                            self.pathwaysWithGenesDictionary[pathwayID] = ['uniprot:'+gene]
                        else:
                            self.pathwaysWithGenesDictionary[pathwayID] = ['entrez:'+gene]                                
                           
                           
                            
                        self.pathwayDictionary[pathwayID] = pathwayName
                        self.pathwayCategory[pathwayID] = "NA"
                        self.geneInfoDictionary['uniprot:'+gene] = mapping
                        #'uniprot:'+
                    else: 
                        listOfGenes = self.pathwaysWithGenesDictionary[pathwayID]
                        if proteinsOrGenes == "proteins":
                            listOfGenes.append('uniprot:'+gene)
                            self.geneInfoDictionary['uniprot:'+gene] = mapping
                        else:
                            listOfGenes.append('entrez:'+gene)
                            self.geneInfoDictionary['entrez:'+gene] = mapping
                            
                        self.pathwaysWithGenesDictionary[pathwayID] = listOfGenes
                 
    
    
    def getMetabolites(self):
           
        reactomeFile = open("../misc/data/reactome/ChEBI2Reactome_All_Levels.txt", encoding="utf-8")
        
        for line in reactomeFile:
            splitline = line.split("\t")
            if len(splitline) > 2:
                if "Homo sapiens" in splitline[5]:
                    metabolite = splitline[0]
                    pathwayID = splitline[1]
                    pathwayName = splitline[3]
                    metabolite = 'chebi:'+metabolite
                    if metabolite not in self.metabolitesWithPathwaysDictionary:
                        self.metabolitesWithPathwaysDictionary[metabolite] = [pathwayID]

                        mapping = {"chebi_id": [metabolite], 
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
            
                         
                         
                         
                         
                         
                        self.metaboliteIDDictionary[metabolite] = mapping
                    else: 
                        listOfPathways = self.metabolitesWithPathwaysDictionary[metabolite]
                        listOfPathways.append(pathwayID)
                        self.metabolitesWithPathwaysDictionary[metabolite] = listOfPathways
                    if pathwayID not in self.pathwayDictionary:
                        self.pathwayDictionary[pathwayID] = pathwayName
                        self.pathwayCategory[pathwayID] = "NA"
                         
    
    def getCommonNameForChebi(self):
        
        for key in self.metaboliteIDDictionary:
            
            value = self.metaboliteIDDictionary[key]
            chebi = value["chebi_id"]
            name = "NA"
            commonName = None
            for each in chebi:
                try:
                    chebiToSearch = str(each)
                    # do we have 'chebi:' prefix?
                    chebi = chebiToSearch.split(":")
                    if(len(chebi) > 1):
                        chebiToSearch = chebi[1]
                    else:
                        chebiToSearch = chebi[0]
                    
                    chebiToSearch2 = libchebipy.ChebiEntity(chebiToSearch)
                    name = chebiToSearch2.get_name()
                    commonName = name
                    
                    
                    #time.sleep(1)
                    if commonName is not None:
                        self.metaboliteCommonName[each] = commonName
                    else:
                        self.metaboliteCommonName[each] = "NA" 
                except:
                    pass
            self.metabolitesWithSynonymsDictionary[key] = [name]  
            if commonName is not None:
                self.metaboliteCommonName[key] = commonName
            else:
                self.metaboliteCommonName[key] = "NA" 
                
                
    '''
    Query Uniprot database to download all files to fill common name of 
    genes by using their RESTFUL APT
    - hmdbdict gene information from HMDB database
    - keggdict gene information from KEGG database
    - wikidict gene information from wikipathways database
    If attributes != None, it will find common name from existed dictionary 
    to fill the Reactome.
    
    By default, all the synonyms are found by querying Uniprot API.
    
    The function should be called:
     1) self.downloadCommonNameFromUniprot()
     2) self.getCommonNameFromUniprot()
    '''
    def downloadCommonNameFromUniprot(self,hmdbdict=None,keggdict=None,wikidict = None):
        reactGeneIds = list()
        otherdatabase = {
            "hmdb": hmdbdict,
            "kegg": keggdict,
            "wiki": wikidict,
            }
        for key in self.geneInfoDictionary:
            splitline = key.split(":")
            reactGeneIds.append(splitline[1])
        #print("Total " + str(len(reactGeneIds)) +" without a common name")
        uniprot_commonName = dict()
        entrez_commonName = dict()
        for key in otherdatabase:
            geneInfo = otherdatabase[key]
            if geneInfo is not None:
                for key2 in geneInfo:
                    #print("gene uniprot takeawy:", key)
                    mapping = geneInfo[key2]
                    uniprot = mapping["UniProt"]
                    commonName = mapping["common_name"]
                    if uniprot != "NA":
                        for id in uniprot:
                            if commonName != "NA":
                                uniprot_commonName[id] = commonName

                    entrez = mapping["entrez"]
                    commonName = mapping["common_name"]
                    if entrez != "NA":
                        for id in entrez:
                            if commonName != "NA":
                                entrez_commonName[id] = commonName
                                  
                                
            #print("Found Uniprot in " + key +": "+ str(len(uniprot_commonName)))
        for id in reactGeneIds:
            #print("********id", id)
            if id in uniprot_commonName:
                name = uniprot_commonName[id]
                mapping = self.geneInfoDictionary['uniprot:'+id]
                mapping["common_name"] = name
                self.geneInfoDictionary['uniprot:'+id] = mapping
                reactGeneIds.remove(id)
                
            if id in entrez_commonName:
                name = entrez_commonName[id]
                mapping = self.geneInfoDictionary['entrez:'+id]
                mapping["common_name"] = name
                self.geneInfoDictionary['entrez:'+id] = mapping
                reactGeneIds.remove(id)    
        
        #print("Unfound genes for name are " + str(len(reactGeneIds)))
        
        Ids = reactGeneIds
        url ="http://www.uniprot.org/uniprot/"
        dir = "../misc/data/Uniprot/"
        self.check_path(dir)
        files = os.listdir(dir)
        num = len(files)
        query = []
        file_dir = []
        files_name = []
        for id in Ids:
           
            if id + ".xml" not in files:
                #print("Downloading ..." + id)
                query.append(url + id +".xml")
                file_dir.append(dir)
                files_name.append(id+".xml")
                #self.download_files(query, dir + id + ".xml")
                #print( str(num) +"/" +str(len(Ids)))
                num = num + 1
        with Pool(50) as p:
            p.starmap(self.download_files,
                  zip(query,file_dir,files_name))
            
            
    '''
    Parse XML files of Uniprot to get common name of genes.
    Fill self.geneInfoDict['common_name']
    '''
    def getCommonNameFromUniprot(self):
        files = os.listdir("../misc/data/Uniprot/")
        path = "../misc/data/Uniprot/"   
        i = 0
        #print('Parsing UniProt files ...')
        for f in files:
            i = i + 1
            #if i % 1000 == 0:
                #print('Processing {} files'.format(i))
            try:
                tree = ET.parse(path + f)
                geneid = f.replace(".xml","")
                root = tree.getroot()
                childs = root.iter("{http://uniprot.org/uniprot}entry")
                
                for child in childs:
                    for child2 in child:
                        childtag = child2.tag.replace("{http://uniprot.org/uniprot}","")
                        if childtag == "gene":
                            for name in child2:
                                type = name.get("type")
                                if type == "primary":
                                    #print(geneid+":"+name.text)
                                    try:
                                        mapping = self.geneInfoDictionary['uniprot:'+geneid]
                                        mapping["common_name"] = "gene_symbol:"+name.text
                                    except KeyError:
                                        print("Raw data does not have this ID ...")
                                        print(geneid)
                                        
                        # we now have uniprot to 'common_name', really gene id.
                        # now we want to grab the NCBI/Entrez 'GeneID'                 
                        if childtag == "dbReference":
                            if child2.get("type") == "GeneID":
                                geneId = child2.get("id")
                                geneId = 'entrez:'+geneId
                                # protein to gene can be 1:n, so they have to be stored as a list
                                # lets check for a value
                                idList = mapping.get("small_e_entrez", None)
                                if(idList == None):
                                    idList = list()
                                    mapping["small_e_entrez"] = idList

                                idList.append(geneId)  
                                
                        
                                        
            except ET.ParseError:
                print("Skip {} ...".format(f))
                pass
             
#     def checkFiles(self):
#         print("In reactome land")
#         files = os.listdir("../misc/data/Uniprot/")
#         print(path.exists("../misc/data/Uniprot/"))
#         print(len(files))
             
# rd = reactomeData()
# rd.getDatabaseFiles()
# rd.getEverything(True)          

#print("In reactome land")
#files = os.listdir("../misc/data/Uniprot/")
#print(path.exists("../misc/data/Uniprot/"))
#print(files)

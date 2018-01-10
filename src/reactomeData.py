import urllib.request 
import urllib.error as ER
import libchebipy
import time
from os import listdir
import xml.etree.ElementTree as ET
class reactomeData():
        
    def __init__(self):
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
        
    def getDatabaseFiles(self):
        
        '''This function gets the files that make up reactome and places them into the reactome folder. 

        '''
        
        urllib.request.urlretrieve("http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt", "../misc/data/reactome/UniProt2Reactome_All_Levels.txt")
        urllib.request.urlretrieve("http://www.reactome.org/download/current/ChEBI2Reactome_All_Levels.txt", "../misc/data/reactome/ChEBI2Reactome_All_Levels.txt")

    def getGenes(self): 
            
        reactomeFile = open("../misc/data/reactome/UniProt2Reactome_All_Levels.txt", encoding="utf-8")
        
        for line in reactomeFile:
            splitline = line.split("\t")       
            if len(splitline) > 2:
                 if "Homo sapiens" in splitline[5]:
                     gene = splitline[0]
                     pathwayID = splitline[1]
                     pathwayName = splitline[3]
                    
                     mapping = { 'kegg': 'NA',
                                'common_name': 'NA',
                                'Ensembl': 'NA', 
                                'HGNC': 'NA', 
                                'HPRD': 'NA', 
                                'NCBI-GeneID': 'NA', 
                                'NCBI-ProteinID': 'NA', 
                                'OMIM': 'NA', 
                                'UniProt': [gene], 
                                'Vega': 'NA', 
                                'miRBase': 'NA', 
                                'HMDB_protein_accession': 'NA',
                                'Entrez': 'NA',
                                'Enzyme Nomenclature': 'NA'}

                     if pathwayID not in self.pathwaysWithGenesDictionary:
                         self.pathwaysWithGenesDictionary[pathwayID] = [gene]
                         self.pathwayDictionary[pathwayID] = pathwayName
                         self.pathwayCategory[pathwayID] = "NA"
                         self.geneInfoDictionary[gene] = mapping
                         
                     else: 
                         listOfGenes = self.pathwaysWithGenesDictionary[pathwayID]
                         listOfGenes.append(gene)
                         self.geneInfoDictionary[gene] = mapping
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
                                     "CAS": "NA"}
            
                         
                         
                         
                         
                         
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
                    chebiToSearch2 = libchebipy.ChebiEntity(chebiToSearch)
                    name = chebiToSearch2.get_name()
                    commonName = name
                    print("Getting ..." + each)
                    print(name)
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
    Find Common Name from all other datasbase geneInforDict
    And return the Uniprot IDs that does not have common name 
    
    '''
    def getCommonNameForGenes1(self,hmdbdict = None,keggdict = None,wikidict = None):
        reactGeneIds = list()
        otherdatabase = {
            "hmdb": hmdbdict,
            "kegg": keggdict,
            "wiki": wikidict,
            }
        for key in self.geneInfoDictionary:
            reactGeneIds.append(key)
        print("Total " + str(len(reactGeneIds)) +" without a common name")
        uniprot_commonName = dict()
        for key in otherdatabase:
            geneInfo = otherdatabase[key]
            if geneInfo is not None:
                for key2 in geneInfo:
                    mapping = geneInfo[key2]
                    uniprot = mapping["UniProt"]
                    commonName = mapping["common_name"]
                    if uniprot != "NA":
                        for id in uniprot:
                            if commonName != "NA":
                                uniprot_commonName[id] = commonName
            print("Found Uniprot in " + key +": "+ str(len(uniprot_commonName)))
        for id in reactGeneIds:
            if id in uniprot_commonName:
                name = uniprot_commonName[id]
                mapping = self.geneInfoDictionary[id]
                mapping["common_name"] = name
                self.geneInfoDictionary[id] = mapping
                reactGeneIds.remove(id)
        
        print("Unfound genes for name are " + str(len(reactGeneIds)))        
        return(reactGeneIds)
    '''
    Query Uniprot database to download all files to fill common name of 
    genes
    
    '''
    def downloadCommonNameFromUniprot(self,Ids):
        url ="http://www.uniprot.org/uniprot/"
        files = listdir("../misc/data/Uniprot/")
        num = len(files)
        for id in Ids:
            print("Downloading ..." + id)
            if id + ".xml" not in files:
                query = url + id +".xml"
                try:
                    urllib.request.urlretrieve(query,"../misc/data/Uniprot/" +id +".xml")
                except ER.HTTPError:
                    print(id + " Not Found ...")
                    continue
                print( str(num) +"/" +str(len(Ids)))
                num = num + 1
    '''
    Parse XML files of Uniprot to get common name of genes.
    Fill self.geneInfoDict['common_name']
    '''
    def getCommonNameForGenes2(self):
        files = listdir("../misc/data/Uniprot/")
        path = "../misc/data/Uniprot/"   
        for f in files:
            print(path+f)
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
                                    mapping = self.geneInfoDictionary[geneid]
                                    mapping["common_name"] = name.text
                                except KeyError:
                                    print("Raw data does not have this ID ...")
                                    print(geneid)
                                #time.sleep(0.1)
             
          
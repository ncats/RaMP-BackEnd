import urllib.request
import time
from multiprocessing import Pool
import os
from MetabolomicsData import MetabolomicsData
'''
This module works to import a variety of data about metabolites from the kegg database.

The kegg database is available here: http://www.genome.jp/kegg/kegg1.html and information
about the api (used in this module) is available here: http://www.kegg.jp/kegg/rest/keggapi.html

'''


class KeggData(MetabolomicsData):

    '''
    KeggData's functions work together to get all required information from the kegg database. 
    
    
    This class contains six functions:
    
        - 1) getPathways()
        - 2) getMetabolites()
        - 3) getSynonyms()
        - 4) getGenes()
        - 5) getGeneInfo()
        - 6) writeToFiles()
        
    These functions must be called in the order listed to work correctly as each consecutive function relies on the previous function(s). 
        
    '''
    
    def __init__(self):
        # constructor
        super().__init__()
        # key: kegg metabolite ID CXXXXX Value: Common Name Assume the first appeared name is the common name
        self.metaboliteCommonName = dict() 
        #key: hsaID for pathway, value: pathway name
        self.pathwayDictionary = dict()
        
        #hsaID for pathway, value: pathway category (only THREE things will be kept: metabolism, cellular process, human disease)
        self.pathwayCategory = dict()
        
        #key: metabolite id , value: list of pathway ids
        self.metabolitesWithPathwaysDictionary = dict()
        
        #key: metabollite id, value: list of synonyms 
        self.metabolitesWithSynonymsDictionary = dict()
        
        #key: pathway id, value: list of gene ids
        self.pathwaysWithGenesDictionary = dict()
        
        #key: gene id, value: gene common name
        self.genesDictionary = dict()
        
        #key: metabolite id, value: chebi
        self.metaboliteIDDictionary = dict()
        
        #only not empty when a catalyzed class exists 
        #empty
        self.metabolitesLinkedToGenes = dict()
        
        #key: gene id, value: list of gene identifiers 
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
        
        '''This function gets the files that make up kegg and places them into the kegg folder. 

        '''
        
        #retrieve the kegg pathway file
        '''
        Download file hsa.txt that has 321 pathways in it
        '''
        pathway_dir = "../misc/data/kegg/"
        pathway_file = "hsa.txt"
        download_url = "http://rest.kegg.jp/list/pathway/hsa"
        
        self.check_path(pathway_dir)
    
        if os.path.exists(pathway_dir) and pathway_file not in os.listdir(pathway_dir):
            self.download_files(download_url,pathway_dir+pathway_file)
            #urllib.request.urlretrieve(download_url, pathway_dir+pathway_file)
            
        tempPathwayDictionary = dict()
        pathwayDictionary = dict()
        compoundDictionary = dict()
        geneDictionary = dict()
        
        print("Downloading pathways...")
        #GET PATHWAYS
        pathwayFile = open("../misc/data/kegg/hsa.txt")
        
        for line in pathwayFile:
            splitline = line.split("\t")       
            #remove the "path:" prefix, it's pointless 
            splitline[0] = splitline[0].replace("path:", "")
            #if the length of the line is greater than one
            if len(splitline) > 1:
                #get information (requires more splitting on "-" since all pathways end with "- [organism]" ) 
                splitagain = splitline[1].split(" - ")
                hsaID = splitline[0]           
                #sometimes the ID needs to have the prefix "hsa" and sometimes it needs to
                #have the prefix "map". My reducing it to just the number it will be easy to 
                #append the correct prefix when necessary. 
                hsaIDnumberonly = hsaID.replace("hsa", "")
                pathway = splitagain[0]
                organism = splitagain[1]
                if hsaIDnumberonly not in tempPathwayDictionary:
                    #place it in the dictionary 
                    tempPathwayDictionary[hsaIDnumberonly] = pathway
        pathwayFile.close()
        #FILTER PATHWAYS BETTER
        each_pathway_dir = "../misc/data/kegg/pathways/"
        self.check_path(each_pathway_dir) # check and create directory
        
        files = os.listdir(each_pathway_dir)
        
        '''
        Get individual pathway file from hsa.txt
        '''
        
        for key in tempPathwayDictionary:    
            url = 'http://rest.kegg.jp/get/' + "path:hsa" + key
            print("url=" +url)
            name = "pathwayhsa"+key+".txt"
            #time.sleep(1)
            if name not in files:
                pathToSavedFile = "../misc/data/kegg/pathways/" + "pathwayhsa" + key + ".txt"
                self.download_files(url, pathToSavedFile)
                onePathwayFile = open(pathToSavedFile)
                # assign the name not all pathways have CLASS
                for line in onePathwayFile:
                    if "CLASS" in line:
                        line = line.replace("CLASS       ", "")
                        splitline = line.split("; ")
                        category = splitline[0]
                        #if category == "Metabolism" or category == "Human Diseases" or category == "Cellular Processes":
                    pathwayDictionary[key] = tempPathwayDictionary[key]
                onePathwayFile.close()
            else:
                pathToSavedFile = "../misc/data/kegg/pathways/" + "pathwayhsa" + key + ".txt"
                # dont need to download anything
                onePathwayFile = open(pathToSavedFile)
                for line in onePathwayFile:
                    if "CLASS" in line:
                        line = line.replace("CLASS       ", "")
                        splitline = line.split("; ")
                        category = splitline[0]
                        #if category == "Metabolism" or category == "Human Diseases" or category == "Cellular Processes":
                    pathwayDictionary[key] = tempPathwayDictionary[key]
                onePathwayFile.close()
            
        print("Downloading compound names...")
        #GET COMPOUNDS IN PATHWAY  
        '''
        Download file to pathwaysWithCompound folder
        Store pathway compound relations ...
        
        '''   
        
        pathway_cpd_dir = "../misc/data/kegg/pathwaysWithCompounds/"
        self.check_path(pathway_cpd_dir)
        files = os.listdir(pathway_cpd_dir)           
        for key in pathwayDictionary:
            url = 'http://rest.kegg.jp/link/cpd/' + "map" + key
            pathToSavedFile = "../misc/data/kegg/pathwaysWithCompounds/" + "cpdmap" + key + ".txt"
            name = "cpdmap"+key+".txt"
            
            print("download files ..." + name)
            compoundFile = None
            if name not in files:
                self.download_files(url, pathToSavedFile)
                compoundFile = open(pathToSavedFile)
            else:
                compoundFile = open(pathToSavedFile)

           
            
            for line in compoundFile:
                line = line.rstrip('\n')
                splitline = line.split("\t")
                if len(splitline) > 1:
                        #remove the "path" prefix and "cpd:" prefix
                        splitline[0] = splitline[0].replace("path:", "")
                        splitline[1] = splitline[1].replace("cpd:", "")
                        compound = splitline[1]
                        compoundDictionary[compound] = "compound"
                
            compoundFile.close()
            
        print("Downloading additional information about compound...")    
        #GET INDIVIDUAL COMPOUND FILES
        cpd_dir = '../misc/data/kegg/compounds/'
        self.check_path(cpd_dir)
        cfiles = os.listdir(cpd_dir)
        for metabolite in compoundDictionary:
            url = 'http://rest.kegg.jp/get/' + metabolite
            file = metabolite +'.txt'
            pathToSavedFile = "../misc/data/kegg/compounds/" + metabolite + ".txt"
            if file not in cfiles:
                self.download_files(url, pathToSavedFile)
            print(url)
       
        
        #GET DICTIONARY OF GENES (no new download)
        found = False
        for key in pathwayDictionary:
            pathToSavedFile = "../misc/data/kegg/pathways/" + "pathwayhsa" + key + ".txt"
            onePathwayFile = open(pathToSavedFile)

                 #There are two cases of what we are looking for:
                 #1) It is a gene that is also on the line marked gene
                 #2) It is a gene that is not on the line marked gene
            for line in onePathwayFile:
                 #A list of things that will signal the end of the genes 
                if "COMPOUND" in line or "REFERENCE" in line or "KO_PATHWAY" in line:
                    found = False
                    break
                 #on line marked gene
                if "GENE" in line:
                    found = True
                    splitline = line.split("  ")
                    geneid = splitline[4]
                    if geneid not in geneDictionary:
                        geneDictionary[geneid] = "gene"
                    continue
                 #not on line marked gene (but gene HAS BEEN found)
                if found:
                    splitline = line.split("  ")
                    geneid = splitline[6]
                    if geneid not in geneDictionary:
                        geneDictionary[geneid] =  "gene"
                     #add geneid to list for this pathway
                    continue
            onePathwayFile.close()
             
        print("Downloading genes...")
        totalgenes = len(geneDictionary)
        genes_dir = '../misc/data/kegg/genes/'
        self.check_path(genes_dir)
        files = os.listdir(genes_dir)
        currentgene = len(files)
        for gene in geneDictionary:
            
            url = "http://rest.kegg.jp/get/" +"hsa:" + gene
            file = gene + '.txt'
            
            pathToSavedFile = "../misc/data/kegg/genes/" + gene + ".txt"
            print("current Kegg gene download:"+file+"---" + str(currentgene) + "/" + str(totalgenes))
            currentgene = currentgene + 1
            if file not in files:
                self.download_files(url, pathToSavedFile)
                
    
    def getPathways_with_genes(self):
        '''
        getDatabaseFiles2(self)
        Following getDatabaseFiles1, this function download links between pathways and genes 
        to new Folder pathwayWithGenes...
        
        1. call self.getPathways() first
        2. Call self.getPathways_with_genes() to get all pathway gene mapping file...
        
        '''            
        path = "../misc/data/kegg/pathwaysWithGenes/"
        if self.check_path(path):
            dir = os.listdir(path)
            for key in self.pathwayDictionary:
                filename = "hsa"+key+".txt"
                if filename not in dir:
                    url ="http://rest.kegg.jp/link/hsa/hsa" + key
                    print("download ... " + url)
                    self.download_files(url,path + "hsa" +key+".txt")
                    
            # Update genes dict
            genedir = os.listdir("../misc/data/kegg/genes/")
            for file in dir:
                filepath = path + file
                pathwayFile = open(filepath)
                for line in pathwayFile:
                    splitline = line.split("\t")
                    if(len(splitline) > 1):
                        pathway = splitline[0].replace("path:hsa","")
                        geneid = splitline[1].replace("hsa:","")
                        geneid = geneid.replace("\n","")
                        geneFile = geneid +".txt"
                        if geneFile not in genedir:
                            url = "http://rest.kegg.jp/get/" +"hsa:" + geneFile.replace(".txt","")
                            
                
                            pathToSavedFile = "../misc/data/kegg/genes/" + geneFile
                         
                        
                            self.download_files(url, pathToSavedFile)
                        
                pathwayFile.close()   
                
                   
        
    
    def getPathways(self):
        '''
        The function getPathways gets a list of all HUMAN pathways in kegg (there are other, non-human pathways in kegg as well)
        
        Overall: the function gets a list of the human pathways from the api in a tab-delimited format. The pathwayID and
        the common name of the pathway are placed in a dictionary (self.pathwayDictionary) with the pathwayID as the key. 
        
        There is no way to get compounds from the kegg api that are listed SPECIFICALLY as human DIRECTLY. The indirect way is
        to find all the human pathways and then find all the compounds linked to the human pathway. This is the approach. 
        
        If we want to find the metabolites, we find pathway first ...
        '''
        
        #create a temporary pathway dictionary that has all human pathways regardless of category
        tempPathwayDictionary = dict()
    
        hsaFile = open("../misc/data/kegg/hsa.txt")
        for line in hsaFile:
                #split into columns via the tab character
                line = line.rstrip('\n')
                splitline = line.split("\t")       
                #remove the "path:" prefix, it's pointless 
                splitline[0] = splitline[0].replace("path:", "")
                #if the length of the line is greater than one
                if len(splitline) > 1:
                    #get information (requires more splitting on "-" since all pathways end with "- [organism]" ) 
                    splitagain = splitline[1].split(" - ")
                    hsaID = splitline[0]           
                    #sometimes the ID needs to have the prefix "hsa" and sometimes it needs to
                    #have the prefix "map". My reducing it to just the number it will be easy to 
                    #append the correct prefix when necessary. 
                    hsaIDnumberonly = hsaID.replace("hsa", "")
                    pathway = splitagain[0]
                    #print("PathwayName: "+pathway)
                    #time.sleep(3)
                    organism = splitagain[1]
                    if hsaIDnumberonly not in tempPathwayDictionary:
                        #place it in the dictionary 
                        tempPathwayDictionary[hsaIDnumberonly] = pathway
        
        #find the pathway category and only keep certain pathways: metabolism, cellular process, human disease
        # Not limited by categoires from 11/17/2017 Update
        files = os.listdir("../misc/data/kegg/pathways")
        for key in tempPathwayDictionary:    
            pathToSavedFile = "../misc/data/kegg/pathways/" + "pathwayhsa" + key + ".txt"
            onePathwayFile = open(pathToSavedFile)
            name = "pathwayhsa"+key+".txt"
            #if name not in files:
            # not all pathways have a categories ...
            foundClass = False
            for line in onePathwayFile:
                line = line.rstrip('\n')
                
                if "CLASS" in line:
                    line = line.replace("CLASS       ", "")
                    splitline = line.split("; ")
                    category = splitline[0]
                    # Initially filter this BUT later is considerred not necessary
                    #if category == "Metabolism" or category == "Human Diseases" or category == "Cellular Processes":
                    self.pathwayDictionary[key] = tempPathwayDictionary[key]
                    self.pathwayCategory[key] = category
                    foundClass = True
            if not foundClass:
                self.pathwayDictionary[key] = tempPathwayDictionary[key]
                self.pathwayCategory[key] = "NA"
            
            
            onePathwayFile.close()
        list1 = self.pathwayDictionary.keys()
        list2 = tempPathwayDictionary.keys()
        print(set(list2) - set(list1))
                         
                        
    def getMetabolites(self):
        '''
        The function getMetabolites creates a list of metabolites that are linked to the human pathways found in the getPathways function. 
        
        Overall: the function takes the keys from in self.pathwayDictionary, queries the kegg api with them, and reports back the 
        metabolites that are linked to those pathways and places them in a list of Metabolites.  
        '''
        
        ##keeping track of current pathway

        
        for key in self.pathwayDictionary:
            pathToSavedFile = "../misc/data/kegg/pathwaysWithCompounds/" + "cpdmap" + key + ".txt"
            onePathwayFile = open(pathToSavedFile) 
            for line in onePathwayFile:
                line = line.rstrip('\n')
                #split into columns via the tab character
                splitline = line.split("\t")
                print(splitline)
                #time.sleep(3)
                if len(splitline) > 1:
                    #remove the "path" prefix and "cpd:" prefix
                    splitline[0] = splitline[0].replace("path:", "")
                    splitline[1] = splitline[1].replace("cpd:", "")
                    mapID = splitline[0]
                    compound = splitline[1]
                    pathway = "hsa" + key
                    #put the compound in the dictionary if it has never been seen before 
                    
                    
                    #If the compound has never been seen before
                    if compound not in self.metabolitesWithPathwaysDictionary:
                        #make a new list and this is the first item
                        listOfPathways = []
                        listOfPathways.append(key)
                        self.metabolitesWithPathwaysDictionary[compound] = listOfPathways
                    
                    #If the compound has been seen before:
                    else:
                        #bring up that compound in a map
                        currentPathwaysForCompound = self.metabolitesWithPathwaysDictionary[compound]
                        #add it to the list of pathways
                        if key not in currentPathwaysForCompound:
                            currentPathwaysForCompound.append(key)
                            self.metabolitesWithPathwaysDictionary[compound] = currentPathwaysForCompound
            onePathwayFile.close()
        
                
                    
                     
                        
    def getSynonymsAndCHEBI(self):
        '''The function getSynonymsAndCHEBI gets all of the "common names" and chebi ids of metabolites found via the getMetabolites function.
        
        Overall: the function takes all the items in the list self.metabolitesWithPathwaysDictionary and queries the api to find the "common names" and chebi id. 
        Since there are multiple common names per metaboliteID the synonyms are stored in a dictionary of lists with the metaboliteID 
        as the key and the synonyms in the list. 
        
        
        In order to get synonyms and chebi the same url must be called, and it is the url that is unique to each compound.
        This is the most time-consuming step of kegg because there are so many compounds. That's why chebi and synonyms 
        are done together. 
        
        '''
        #We will iterate through all the compounds and and find their "human-interpretable"
#        names. Many/(All?) of the compounds have multiple names so we will need to store them all.

#        Here the for loop is looking through the text file for certain key words. The keyword
#        "NAME" will signal that the lines coming up will be synonyms for the compound of 
#        interest. There are multiple lines that signal that the compound lines have ended.
#        The script checks for those and stops reading in at those lines.
        
        ##keeping track of current metabolite

        
     
        
        
        
        found = False
        for metabolite in self.metabolitesWithPathwaysDictionary: 
            
            
            
            
            
            
            pathToSavedFile = "../misc/data/kegg/compounds/" + metabolite + ".txt"
            compound = open(pathToSavedFile)
            compoundList = []
            #for every line
            for line in compound:
                line = line.rstrip('\n')          
                #check if the loop should end for this file because it onto the next section
                commonName = None
                if "FORMULA" in line or "COMMENT" in line or "REMARK" in line or "REACTION" in line or "SEQUENCE" in line or "DBLINKS" in line:
                    found = False
                    break           
                #if not, and if what we are looking for is found
                if "NAME" in line:
                    found = True
                    line = line.replace(" ", "")
                    line = line.replace("NAME", "")
                    #place it in the list
                    compoundList.append(line)
                    commonName = line
                    if commonName is not None:
                        self.metaboliteCommonName[metabolite] = commonName
                    else:
                        self.metaboliteCommonName[metabolite] = "NA"
                    
                    #get next line
                    continue            
                #continue placing the lines into the list as long as found is true 
                if found:
                    line = line.replace(" ", "")
                    compoundList.append(line)
                    continue
            
            self.metabolitesWithSynonymsDictionary[metabolite] = compoundList
           
            ####CHEBI####
            
            listOfChebi = []
            
            metaboliteMapping = {"chebi_id": "NA", 
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
            
            for line in compound:
                line = line.rstrip('\n')
                if "ChEBI" in line:
                    splitline = line.split(" ")
                    length = len(splitline)

                    #first chebi id
                    #(There may be more than one chebi id)
                    start = 13
                    
                    while start < length:
                        chebiid = splitline[start]
                        start = start +    1
                        listOfChebi.append(chebiid)
                        
                if "CAS" in line:
                    casid = line[line.find(":")+1:]
                    casid = casid.replace(" ","")
                    metaboliteMapping["CAS"] = casid
                if "PubChem" in line:
                    pubchem = line[line.find(":")+1:]
                    pubchem = pubchem.replace(" ","")
                    metaboliteMapping["pubchem_compound_id"] = pubchem    
                metaboliteMapping["chebi_id"] = listOfChebi
                
            metaboliteMapping["kegg_id"] = metabolite
            self.metaboliteIDDictionary[metabolite] = metaboliteMapping
            compound.close()
                         
                        
                         
           

    def getGenes(self):
        
        '''The getGenes function finds genes for the pathways queried.'''
        
        
        ##keeping track of current pathway
   
        kegg_id = []
        found = False
        for key in self.pathwayDictionary:
             
            geneList = []
             
            #keeping track of number of genes 
     
            pathToSavedFile = "../misc/data/kegg/pathways/" + "pathwayhsa" + key + ".txt"
             
             
             
             
            onePathwayFile = open(pathToSavedFile)
                 
                 
            geneList = []     
                #There are two cases of what we are looking for:
                #1) It is a gene that is also on the line marked gene
                #2) It is a gene that is not on the line marked gene
            for line in onePathwayFile:
                line = line.rstrip('\n')
                 
                #A list of things that will signal the end of the genes 
                if "COMPOUND" in line or "REFERENCE" in line or "KO_PATHWAY" in line:
                    found = False
                    break
                #on line marked gene
                if "GENE" in line:
                     
                    mapping = {'kegg': 'NA',
                         'common_name': 'NA',
                         'Ensembl': 'NA', 
                         'HGNC': 'NA', 
                         'HPRD': 'NA', 
                         'NCBI-GeneID': 'NA', 
                         'NCBI-ProteinID': 'NA', 
                         'OMIM': 'NA', 
                         'UniProt': 'NA', 
                         'Vega': 'NA', 
                         'miRBase': 'NA', 
                         'HMDB_protein_accession': 'NA',
                         'Entrez' : 'NA',
                         'Enzyme Nomenclature': 'NA'}
                     
                     
                    found = True
                     
                    splitline = line.split("  ")
                     
                    geneid = splitline[4]
                    genefullname = splitline[5]
                    #remove [KO:K05757]-like attachment to the full name -- unsure what the point of it is 
                    genefullname = genefullname.split(" [")
                    genefullname = genefullname[0]


                    #add to dictionary that keeps track of gene ids: key: geneid, value: genefullname 
                    mapping['common_name'] = [genefullname]
                    mapping['kegg'] = geneid   
                    kegg_id.append(geneid)    
                    if geneid not in self.geneInfoDictionary:
                        print(geneid)
                        #time.sleep(3)
                        self.geneInfoDictionary[geneid] = mapping
                     
                    #add geneid to list for this pathway
                    if geneid not in geneList:
                        geneList.append(geneid)
                                 
                    continue
                #not on line marked gene (but gene HAS BEEN found)
                if found:
                    print(splitline)
                    #time.sleep(3)
                    mapping = {'kegg': 'NA',
                         'common_name': 'NA',
                         'Ensembl': 'NA', 
                         'HGNC': 'NA', 
                         'HPRD': 'NA', 
                         'NCBI-GeneID': 'NA', 
                         'NCBI-ProteinID': 'NA', 
                         'OMIM': 'NA', 
                         'UniProt': 'NA', 
                         'Vega': 'NA', 
                         'miRBase': 'NA', 
                         'HMDB_protein_accession': 'NA',
                         'Entrez' : 'NA',
                         'Enzyme Nomenclature': 'NA'}
                     
                    splitline = line.split("  ")
                    geneid = splitline[6]
                    #print(splitline)
                    print(geneid)
                    #time.sleep(3)
                    genefullname = splitline[7]
                    #remove [KO:K05757]-like attachment to the full name -- unsure what the point of it is 
                    genefullname = genefullname.split(" [")
                    genefullname = genefullname[0]
                    
                    #add to dictionary that keeps track of gene ids: key: geneid, value: genefullname 

                    #add to dictionary that keeps track of gene ids: key: geneid, value: genefullname 
                    mapping['common_name'] = [genefullname]       
                    mapping['kegg'] = geneid
                    if geneid not in self.geneInfoDictionary:
                        self.geneInfoDictionary[geneid] = mapping
                    #add geneid to list for this pathway
                    if geneid not in geneList:
                        geneList.append(geneid)
                     
                    continue

            if key not in self.pathwaysWithGenesDictionary:
                self.pathwaysWithGenesDictionary[key] = geneList
                
        
        print(kegg_id)
        print(len(kegg_id))
        
    '''
    self.getGenes + self.getGeneInfo are not enough to map all genes to pathways
    Here adds third functions following the two above
    This function parses all files of pathwaysLinkTogenes to build up pathway-genes relations
    '''
    def getPathwayLinkedToGene(self):
        files = os.listdir("../misc/data/kegg/pathwaysWithGenes/")
        for file in files:
            pathwayFile = open("../misc/data/kegg/pathwaysWithGenes/" + file)
            
            for line in pathwayFile:
                splitline = line.split("\t")
                if(len(splitline) >1):
                    print(splitline)
                    pathway = splitline[0].replace("path:hsa","")
                    geneid = splitline[1].replace("hsa:","")
                    geneid = geneid.replace("\n","")
                    if pathway in self.pathwaysWithGenesDictionary:
                        geneList = self.pathwaysWithGenesDictionary[pathway]
                        if geneid not in geneList:
                            print("Add extra gene to list ..." + geneid)
                            geneList.append(geneid)
                            self.pathwaysWithGenesDictionary[pathway] = geneList
                    else:
                        print("Find new pathway " + pathway)
                        print("Add gene " + geneid)
                        print(self.pathwaysWithGenesDictionary)
                        self.pathwaysWithGenesDictionary[pathway] = [geneid]
            
            pathwayFile.close()        
                #time.sleep(3)
            
    def getGeneInfo(self):
        
        '''
        Get information for the genes. 
        
        '''

        
        for gene in self.geneInfoDictionary:
            pathToSavedFile = "../misc/data/kegg/genes/" + gene + ".txt"
            found = False
            geneFile = open(pathToSavedFile)
            keywords = ['Ensembl', 'HGNC', 'HPRD', 'NCBI-GeneID', 'NCBI-ProteinID', 'OMIM', 'UniProt', 'Vega', 'miRBase']
            mapping = self.geneInfoDictionary[gene]
                           
                
            for line in geneFile:
                line = line.rstrip('\n')
                if "NAME" in line:
                    print(line)
                    splitline = line.split("   ")
                    names = splitline[len(splitline) - 1]
                    names = names.split(",")
                    print(names)
                    mapping = self.geneInfoDictionary[gene]
                    for name in names:
                        name = name.replace(" ","")
                        print(name)
                        mapping["common_name"].append(name)
                        self.geneInfoDictionary[gene] = mapping 
                        #print(splitline)
                        #time.sleep(3)
                if "AASEQ" in line or "NTSEQ" in line or "STRUCTURE" in line:
                    found = False
                    break
                if "DBLINKS" in line:
                    found = True
                    splitline = line.split(" ")
                    key = splitline[5]
                    #chop off last character (a colon)
                    key = key[:-1]
                    value = splitline[6]
                    print(key+":"+value)
                    #time.sleep(3)
                    if key in keywords:
                        if key in ['Ensembl', 'UniProt']:
                            mapping[key] = [value]
                        else:
                            mapping[key] = value
                    continue
                if found:
                    splitline = line.split(" ")
                    key = splitline[12]
                    #chop off last character (a colon)
                    key = key[:-1]
                    value = splitline[13]
                    print(key+":"+value)
                    #time.sleep(3)
                     
                     
                    if key in keywords:
                        if key in ['Ensembl', 'UniProt']:
                            mapping[key] = [value]
                        else:
                            mapping[key] = value
                    continue
            
            
            self.geneInfoDictionary[gene] = mapping
            geneFile.close()
             

            
#This line exists for sphinx to work properly (documentation): http://autoapi.readthedocs.io/                
__all__ = ['KeggData']
                    
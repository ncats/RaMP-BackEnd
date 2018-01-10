import xml.etree.ElementTree as ET
import urllib.request
from xml.etree.ElementTree import ElementTree
import codecs
from test.test_iterlen import NoneLengthHint
import zipfile
import time
import os
from MetabolomicsData import MetabolomicsData

class hmdbData(MetabolomicsData):
    
    '''
    HMDBData's functions work together to get all required information from the hmdb database.
    r
    The hmdb database can be queried by parsing an xml file which contains all the information for the database.
    The xml file for metabolites can be obtained here: http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
	The xml file for proteins can be obtained here: http://www.hmdb.ca/system/downloads/current/hmdb_proteins.zip
    
    This file has already been downloaded and can be found in the data folder of this package. Its location in the package 
    is hardcoded into function calls. 
    
    Parsing xml files is not straight-forward if you do not have experience. 
    Here is a tutorial on xml file parsing: https://pymotw.com/2/xml/etree/ElementTree/parse.html
    
    This class contains five functions:
    
        - 1) getMetaboliteOtherIDs()
        - 2) getPathwaysandSynonyms()
        - 3) getGenes()
        - 4) getBiofluidCellularLocationDisease()
        - 5) getPathwayLinkedToGene()
        - 6) WriteToFiles()
    
    In summary, the first four functions parse the xml file, while the last function takes the information acquired by parsing and writes it 
    to sql files for the RAMP database. 
    
    The functions are mostly independent of each other. They do not rely on one another due to the nature of the hmdb database (xml tree).
    However, the final function *IS* dependent on the previous four functions. 
    
    Due to the structure of the data in hmdb's database (xml tree) it is often easier and quicker to get some information together (for example, pathways and synonyms) 
    and therefore a variety of information gathering is often grouped into one function. However, there has been some effort made to separate the different information gathering
    into separate functions for easier readability. The time-intensive step for hmdb information gathering is simply opening/parsing through the xml file
    so it is better to limit the number of this this occurs -- this is the benefit of grouping information together in one function. The drawback is readability (the code can get messy
    and hard to follow). 
    
    '''
    
    def __init__(self):
        
        super().__init__()
        #self.tree = ET.parse('../misc/data/hmdb/hmdb_metabolites.xml')
        #self.tree2 = ET.parse('../misc/data/hmdb/hmdb_proteins.xml')
        ####DICTIONARIES IN COMMON WITH OTHER CLASSES######################################
        # common name dictionary Key: HMDB ID Value: Common Name
        self.metaboliteCommonName = dict()
        #pathway dictionary. Key: hsaID for pathway, Value: pathway name
        self.pathwayDictionary = dict()
        # pathway id mapping Key SMP ID Value: Kegg id
        self.SMPToKegg = dict()
        #hsaID for pathway, value: pathway category (all will be "NA" for hmdb)
        self.pathwayCategory = dict()
        
        #key: metabolite id , value: list of pathway id
        self.metabolitesWithPathwaysDictionary = dict()
        
        #key: metabollite id, value: list of synonyms 
        self.metabolitesWithSynonymsDictionary = dict()
        
        #key: metabolite id, value: mapping of other ids
        self.metaboliteIDDictionary = dict()
        
        # key: pathway SMP id, value: list of gene HMDBP id
        self.pathwaysWithGenesDictionary = dict()      
        

        #key: gene id, value: list of FOUR gene identifiers 
        self.geneInfoDictionary = dict()
        
        #only not empty when a catalyzed class exists 
        #key: matabole, value: list of genes
        self.metabolitesLinkedToGenes = dict()
        
        ###################################################################
        
        #stays empty for this class
        self.pathwayOntology = dict()
        
        #key: metabolite id, value: list of metabolite locations
        self.biofluidLocation = dict()
        
        #key: biofluid location, value: the string "placeholder"
        self.biofluid = dict()
        
        #key: metabolite id, value: list of cellular locations
        self.cellularLocation = dict()
        
        #key: cellular location, value: the string "placeholder"
        self.cellular = dict()
        
        #key: metaboliteID, value: exo/endo/Drug metabolite
        self.exoEndoDictionary = dict()
        
        #key: origin, value: exo/endo/ Drug ,etc.
        self.exoEndo = dict()
        #key: metaboliteID, value: tissue_locations
        self.tissueLocation = dict()
        
        #key: tissue location, value : "placeholder"
        self.tissue = dict()
        self.idDictForMetabolite = dict()

        
        
    def getDatabaseFiles(self):
        '''
		This function gets the files that make up hmdb and places them into the hmdb folder. 
		
        '''
        
        file_metabolites = "hmdb_metabolites.zip"
        file_proteins = "hmdb_proteins.zip"
        download_url = "http://www.hmdb.ca/system/downloads/current/"
        dir = "../misc/data/hmdb/"
        if not os.path.exists(dir):
            try:
                os.makedirs(dir) # check if the directory exists, create one if not
            except OSError as e: # Trap the OS error and show the embedded error code
                if e.errno != errno.EEXIST:
                    raise				
        else:
            if file_metabolites in os.listdir(dir) and file_proteins in os.listdir(dir):
                print("Files are already downloaded ...")
                return
 
        self.download_files(download_url+file_metabolites, dir+file_metabolites)
        self.download_files(download_url+file_proteins, dir+file_proteins)
        with zipfile.ZipFile(dir+file_metabolites,"r") as zip_ref:
            zip_ref.extractall(dir)
            
        with zipfile.ZipFile(dir+file_proteins,"r") as zip_ref:
            zip_ref.extractall(dir)
        
        os.remove(dir+file_metabolites)
        os.remove(dir+file_proteins)
                            
        
    def getMetaboliteOtherIDs(self):     
        '''
        This functions finds a number of alternative ids for the main metabolite identifier and places them into: 
            
            - self.metaboliteIDDictionary. 
            - self.metaboliteCommonName
            '''
        print("Start parsing ...")
        now = time.time()
        tree = ET.parse('../misc/data/hmdb/hmdb_metabolites.xml')
        root = tree.getroot()
        print("Finish parsing ..." + str(time.time() - now))   
        metabohmdbid = "not yet found"
               
        #we will need to iterate through the xml tree to find the information we are looking for
        for metabolite in root:
        #XML trees sometimes have namespace prefixes on the nodes. They need to be removed for parsing.
            #That's the point of the find and replace. We are removing the namespace string "{http://www.hmdb.ca}"

            metabolitetag = metabolite.tag.replace("{http://www.hmdb.ca}", "")
            
            
            #find the accession number (metabolite id)
            if metabolitetag == "metabolite":
                
                
            ##here are some of the things we will be looking for in the xml tree
                mapping = {"chebi_id": "NA", 
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
                commonName = None
                #find other ids for metabolite 
                for child in metabolite:
                    childtag = child.tag.replace("{http://www.hmdb.ca}", "")
                    if childtag == "name":
                        commonName = child.text
                        #print(metabohmdbid)
                        #print(commonName)
                        #time.sleep(0.5)
                    #THIS WILL BE THE KEY
                    if childtag == "accession":
                        metabohmdbid = child.text 
                        mapping["hmdb_id"] = [metabohmdbid]
                        
                    if childtag == "chebi_id":
                        chebiid = child.text  
                        if chebiid == None:
                            chebiid = "NA"   
                        mapping["chebi_id"] = [chebiid]
                       
                    if childtag == "drugbank_id":
                        drugbankid = child.text
                        if drugbankid == None:
                            drugbankid = "NA"
                        mapping["drugbank_id"] = drugbankid
                    
                    if childtag == "drugbank_metabolite_id":
                        drugbankmetaboliteid = child.text
                        if drugbankmetaboliteid == None:
                            drugbankmetaboliteid = "NA"
                        mapping["drugbank_metabolite_id"] = drugbankmetaboliteid
                     
                    if childtag == "phenol_explorer_compound_id":
                        drugbankmetaboliteid = child.text
                        if drugbankmetaboliteid == None:
                            drugbankmetaboliteid = "NA" 
                        mapping["phenol_explorer_compound_id"] = drugbankmetaboliteid                      
                        
                    if childtag == "phenol_explorer_metabolite_id":
                        phenolexplorermetaboliteid = child.text
                        if phenolexplorermetaboliteid == None:
                            phenolexplorermetaboliteid = "NA"
                        mapping["phenol_explorer_metabolite_id"] = phenolexplorermetaboliteid
                   
                    if childtag == "foodb_id":
                        fooddbid = child.text
                        if fooddbid == None:
                            fooddbid = "NA"
                        mapping["foodb_id"] = fooddbid
                    
                    if childtag == "knapsack_id":
                        knapsackid = child.text
                        if knapsackid == None:
                            knapsackid= "NA"
                        mapping["knapsack_id"] = knapsackid
                       
                    if childtag == "chemspider_id":
                        chemspiderid = child.text
                        if chemspiderid == None:
                            chemspiderid = "NA"
                        mapping["chemspider_id"] = chemspiderid
                        
                    if childtag == "kegg_id":
                        keggid = child.text
                        if keggid == None:
                            keggid = "NA"
                        mapping["kegg_id"] = keggid
                        
                    if childtag == "biocyc_id":
                        biocycid = child.text
                        if biocycid == None:
                            biocycid = "NA"
                        mapping["biocyc_id"] = biocycid
                    
                    if childtag == "bigg_id":
                        biggid = child.text
                        if biggid == None:
                            biggid = "NA"
                        mapping["bigg_id"] = biggid
                       
                    if childtag == "wikipidia":
                        wikipedia = child.text
                        if wikipedia == None:
                            wikipedia = "NA"
                        mapping["wikipidia"] = wikipedia
                 
                    if childtag == "nugowiki":
                        nugowiki = child.text
                        if nugowiki == None:
                            nugowiki = "NA"
                        mapping["nugowiki"] = nugowiki
                       
                    if childtag == "metagene":
                        metagene = child.text
                        if metagene == None:
                            metagene = "NA"
                        mapping["metagene"] = metagene                   
                        
                    if childtag == "metlin_id":
                        metlin_id = child.text
                        if metlin_id == None:
                            metlin_id = "NA"
                        mapping["metlin_id"] = metlin_id     
                        
                    if childtag == "pubchem_compound_id":
                        pubchem_compound_id = child.text
                        if pubchem_compound_id == None:
                            pubchem_compound_id = "NA"
                        mapping["pubchem_compound_id"] = pubchem_compound_id
                           
                    if childtag == "het_id":
                        het_id = child.text
                        if het_id == None:
                            het_id = "NA"
                        mapping["het_id"] = het_id
                    
                    if childtag == "cas_registry_number":
                        cas_id = child.text
                        if cas_id == None:
                            cas_id = "NA"
                        mapping["CAS"] = cas_id
                        #time.sleep(3)
                    

                 
                                    
             
                        
                    
                #place all id information in this mapping        
                if metabohmdbid not in self.metaboliteIDDictionary:
                    self.metaboliteIDDictionary[metabohmdbid] = mapping
                if commonName is not None:
                    self.metaboliteCommonName[metabohmdbid] = commonName
                else:
                    self.metaboliteCommonName[metabohmdbid] = "NA"

    
    def getPathwaysandSynonyms(self):
        '''
        This functions finds pathways and synonyms for the metabolites and places them in:
            
            - self.metabolitesWithPathwaysDictionary
            - self.metabolitesWithSynonymsDictionary
            
        Additionally it creates a mapping between the pathwayid and the pathway name and places it in:
        
            - self.pathwayDictionary
            
        
        '''
        
  
        tree = ET.parse('../misc/data/hmdb/hmdb_metabolites.xml')
        root = tree.getroot()
        
        
        
        
        
        metabohmdbid = "not yet found"
        
        
        #we will need to iterate through the xml tree to find the information we are looking for
        for metabolite in root:
          
            #XML trees sometimes have namespace prefixes on the nodes. They need to be removed for parsing.
            #That's the point of the find and replace. We are removing the namespace string "{http://www.hmdb.ca}"

            metabolitetag = metabolite.tag.replace("{http://www.hmdb.ca}", "")
            # Check if this metabolite has pathway
            haspathway = False
            #find the accession number (metabolite id)
            if metabolitetag == "metabolite":
                
                #find other ids for metabolite 
                 for child in metabolite:
                    childtag = child.tag.replace("{http://www.hmdb.ca}", "")
                    
                    #THIS WILL BE THE KEY
                    if childtag == "accession":
                        metabohmdbid = child.text 
                        print("Getting ..." + metabohmdbid)    
                   #find the pathways 
                    
                    listOfPathways = []
                    if childtag =="pathways":
                        
                        for pathway in child:
                            pathwayName = None
                            smpid = None
                            keggid = None
                            
                            for info in pathway:
                                infotag = info.tag.replace("{http://www.hmdb.ca}", "")
                                
                                if infotag == "name":
                                    if info.text is not None:
                                        pathwayName = info.text
                                if infotag == "smpdb_id":
                                    if info.text is not None:
                                        smpid = info.text
                                        listOfPathways.append(smpid)
                                     #place all pathways in a dictionary with the pathwayid as the key and the common name as the value
                                        if smpid not in self.pathwayDictionary:
                                            self.pathwayDictionary[smpid] = pathwayName
                                            self.pathwayCategory[smpid] = "NA"
                                '''       
                                if infotag == "kegg_map_id":
                                    keggid = info.text
                                    if keggid is not None:
                                        keggid = keggid.replace("map","")
                                        listOfPathways.append(keggid)
                                        # map kegg id with smpid
                                        if smpid not in self.SMPToKegg and smpid is not None:
                                            self.SMPToKegg[smpid] = keggid
                                '''            
                            
                            
                            #time.sleep(3)
                                        # map kegg id with name
                                        #if keggid not in self.pathwayDictionary:
                                        #    self.pathwayDictionary[keggid] = pathwayName 
                        
                                              
                                    
                        if len(listOfPathways) >0 :           
                            self.metabolitesWithPathwaysDictionary[metabohmdbid] = listOfPathways

                    #find the synonyms                    
                    listOfSynonyms = []
                    
                    if childtag  == "synonyms":
                        for synonym in child:
                            synonym = synonym.text
                            listOfSynonyms.append(synonym)
                    
                        #place all synonyms in a dictionary 
                        self.metabolitesWithSynonymsDictionary[metabohmdbid] = listOfSynonyms
        
                             
                           
             
    def getGenes(self):
        '''
        This function finds genes linked to metabolites and places them in:
        
            -self.MetabolitesLinkedToGenes
            
        Additionally, it links the uniprotid to the gene name and place it in:
            
            -self.geneInfoDictionary
        
        And, finally, it finds other ids for every gene and places this in:
        
            -self.geneInfoDictionary
        
        '''
        
        #key: metabolite id, value: list of gene ids
        self.metabolitesLinkedToGenes = dict()
        
        tree = ET.parse('../misc/data/hmdb/hmdb_metabolites.xml')
        root = tree.getroot()
  
        metabohmdbid = "not yet found"
        
        
        #we will need to iterate through the xml tree to find the information we are looking for
        for metabolite in root:
          
            #XML trees sometimes have namespace prefixes on the nodes. They need to be removed for parsing.
            #That's the point of the find and replace. We are removing the namespace string "{http://www.hmdb.ca}"

            metabolitetag = metabolite.tag.replace("{http://www.hmdb.ca}", "")

            
            listofgenes = []
            
            
            
             
            if metabolitetag == "metabolite":
                
  
                #find other ids for metabolite 
                for child in metabolite:
                    childtag = child.tag.replace("{http://www.hmdb.ca}", "")

                    #THIS WILL BE THE KEY
                    if childtag == "accession":
                        metabohmdbid = child.text 
                        
                    
                    if childtag == "protein_associations":
                        
                        for protein in child:
                            
                            mapping = { 'kegg': 'NA',
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
                                       'Entrez': 'NA',
                                       'Enzyme Nomenclature': 'NA'}
                            
                            
                            
                            
                            for proteininfo in protein:
                                proteininfotag = proteininfo.tag.replace("{http://www.hmdb.ca}", "")
                                #print(proteininfotag)
                                #print(proteininfo.text)
                                #time.sleep(1)
                                if proteininfotag == "uniprot_id":
                                    uniprotid = proteininfo.text
                                    mapping["UniProt"] = [uniprotid]
                                    
                                    
                                if proteininfotag == "protein_accession":
                                    proteinacc = proteininfo.text
                                    
                                    if proteinacc is None:
                                        proteinacc = "NA"
                                    mapping["HMDB_protein_accession"] = proteinacc
                                    #self.geneInfoDictionary[proteinacc] = mapping 
                                    
                                    listofgenes.append(proteinacc)
                                    
                                    
                                if proteininfotag == "gene_name":
                                    genename = proteininfo.text
                                    if genename is None:
                                        genename = "NA"
                                        
                                    mapping["common_name"] = genename
                            self.geneInfoDictionary[proteinacc] = mapping
                                
                                    
                                
                         
                        #using the uniprot id as the gene id 
            
             
            self.metabolitesLinkedToGenes[metabohmdbid] = listofgenes
   
                
        
        
        
        
        
    def getBiofluidCellularLocationDisease(self):
        
        '''This function finds biofluid and cellular location infromation for every metabolite and places them in:
        
            -self.cellularLocation
            -self.biofluidLocation
            
        Additionally, a running list of all biofluid and cellular locations are kept:
        
            -self.cellular
            -self.biofluid
        
        
        
        '''
        
        tree = ET.parse('../misc/data/hmdb/hmdb_metabolites.xml')
        root = tree.getroot()
  
        metabohmdbid = "not yet found"
        
        
        #we will need to iterate through the xml tree to find the information we are looking for
        for metabolite in root:
          
            #XML trees sometimes have namespace prefixes on the nodes. They need to be removed for parsing.
            #That's the point of the find and replace. We are removing the namespace string "{http://www.hmdb.ca}"

            metabolitetag = metabolite.tag.replace("{http://www.hmdb.ca}", "")

            #find the accession number (metabolite id)
            if metabolitetag == "metabolite":
               
                #find other ids for metabolite 
                 for child in metabolite:
                    childtag = child.tag.replace("{http://www.hmdb.ca}", "")

                    #THIS WILL BE THE KEY
                    if childtag == "accession":
                        metabohmdbid = child.text 
                    
                    if childtag == "ontology":
                        for cellularlocations in child:
                            cellularlocationstag = cellularlocations.tag.replace("{http://www.hmdb.ca}", "")
                            if cellularlocationstag == "cellular_locations":
                                listOfInfo = []
                                for cellularlocation in cellularlocations:
                                    cellularlocationtext = cellularlocation.text
                                    
                                    if cellularlocationtext not in self.cellular:
                                        self.cellular[cellularlocationtext] = "placeholder"
                                    
                                    listOfInfo.append(cellularlocationtext)
                                self.cellularLocation[metabohmdbid] = listOfInfo
                            
                            if cellularlocationstag == "origins":
                                listOfInfo = []
                                for origin in cellularlocations:
                                    origintext = origin.text
                                    if origintext not in self.exoEndo:
                                        self.exoEndo[origintext] = "placeholder"
                                    listOfInfo.append(origintext)
                                self.exoEndoDictionary[metabohmdbid] = listOfInfo           
                                
                                
                    biofluidList = []
                    if childtag == "biofluid_locations":
                         for biofluid in child: 
                             biofluidtext = biofluid.text
                             biofluidList.append(biofluidtext)   
                             if biofluidtext not in self.biofluid:
                                 self.biofluid[biofluidtext] = "placeholder"
                                                       
                         self.biofluidLocation[metabohmdbid] = biofluidList
                    tissueList = []
                    if childtag == "tissue_locations":
                        for tissuelocation in child:
                            tissuetext = tissuelocation.text
                            tissueList.append(tissuetext)
                            if tissuetext not in self.tissue:
                                self.tissue[tissuetext] = "placeholder"
                        self.tissueLocation[metabohmdbid] = tissueList
    def getAllId(self):
        tree = ET.parse('../misc/data/hmdb/hmdb_metabolites.xml')
        root = tree.getroot()
        for metabolite in root: 
            metatag = metabolite.tag.replace("{http://www.hmdb.ca}","")
            if metatag == "metabolite":
                for child in metabolite:
                    childtag = child.tag.replace("{http://www.hmdb.ca}", "")
                    if "_id" in childtag and child.text is not None and childtag not in self.idDictForMetabolite:
                        self.idDictForMetabolite[childtag] = []
                        self.idDictForMetabolite[childtag].append(child.text)
                    elif "_id" in childtag and child.text is not None and childtag in self.idDictForMetabolite:
                        self.idDictForMetabolite[childtag].append(child.text)                      
    '''
    Found protein file at 11/1/2017 From HMDB
    Initially leave this dict empty
    Try to link hmdb gene to pathways in this case
    Key: pathway id SMP_NUM Value: Geneid HMDBP_NUM
    This function fills dictionary:
        self.PathwayLinkedToGene key: pathway id Value: HMDBP ID
        self.pathwayDictionary key pathway id Value: Name
    '''
    def getPathwaysLinkedToGene(self):
        tree = ET.parse('../misc/data/hmdb/hmdb_proteins.xml')
        keggnum = 0
        root = tree.getroot() 
        for protein in root:
            accession = protein.find('{http://www.hmdb.ca}accession')
            proteintag = protein.tag.replace("{http://www.hmdb.ca}",'')
            for pathways in protein.iter('{http://www.hmdb.ca}pathways'):
                for pathway in pathways:
                    pathwaytag = pathway.tag.replace('{http://www.hmdb.ca}','')
                    print(pathwaytag)
                    pathwayName = pathway.find('{http://www.hmdb.ca}name').text
                    smpid = pathway.find('{http://www.hmdb.ca}smpdb_id').text
                    if smpid is not None:
                        if pathwayName is not None:
                            self.pathwayDictionary[smpid] = pathwayName
                            self.pathwayCategory[smpid] = 'NA'
                    keggid = None
                    #print(pathwayName)
                    #print(smpid)
                    #time.sleep(3)

                    for info in pathway:
                        print(pathwayName)
                        print(smpid)
                        #time.sleep(1)
                        infotag = info.tag.replace('{http://www.hmdb.ca}','')
                        if infotag == "name":
                            pathwayName = info.text
                        if infotag == 'smpdb_id':
                            if info.text is not None:
                                smpid = info.text
                                print(smpid)
                                print(pathwayName)
                                if smpid not in self.pathwayDictionary and smpid is not None:
                                    self.pathwayDictionary[smpid] = pathwayName
                                if smpid not in self.pathwaysWithGenesDictionary:
                                    self.pathwaysWithGenesDictionary[smpid] = []
                                    self.pathwaysWithGenesDictionary[smpid].append(accession.text)
                                else:
                                    if accession.text not in self.pathwaysWithGenesDictionary[smpid]:
                                        self.pathwaysWithGenesDictionary[smpid].append(accession.text)
                    
                    #time.sleep(3)                    

  
    def findMetabolitesWithPathways(self):
        tree = ET.parse('../misc/data/hmdb/hmdb_metabolites.xml')
        root = tree.getroot()

        for metabolite in root:
            metaboliteid = metabolite.find('{http://www.hmdb.ca}accession')
            print(metaboliteid.text)

            time.sleep(3)

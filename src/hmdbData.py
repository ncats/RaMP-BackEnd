#import xml.etree.ElementTree as ET
import urllib.request
#from xml.etree.ElementTree import ElementTree
from lxml import etree as ET
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
        # Define file name and url for downloading file
        file_metabolites = "hmdb_metabolites.zip"
        file_proteins = "hmdb_proteins.zip"
        download_url = "http://www.hmdb.ca/system/downloads/current/"
        dir = "../misc/data/hmdb/"
        # check if this path exists
        self.check_path(dir)
        if file_metabolites.replace('.zip', '.xml') not in os.listdir(dir) or file_proteins.replace('.zip', '.xml') not in os.listdir(dir):
            # Download files from given url and path
            print('####### Downloading HMDB source file #######')
            self.download_files(download_url+file_metabolites, dir+file_metabolites)
            self.download_files(download_url+file_proteins, dir+file_proteins)
            # Open zip file if file are downloaded.
            with zipfile.ZipFile(dir+file_metabolites,"r") as zip_ref:
                zip_ref.extractall(dir)
                
            with zipfile.ZipFile(dir+file_proteins,"r") as zip_ref:
                zip_ref.extractall(dir)
        else:
            print('HMDB source files are ready ...')
            
        
                            
        
    def getMetaboliteOtherIDs(self,tree = None,dir = 'hmdb_metabolites.xml'):     
        '''
        This functions finds a number of alternative ids for the main metabolite identifier and places them into: 
            
            - self.metaboliteIDDictionary. 
            - self.metaboliteCommonName
            
        param elementTree tree A parsed XML HMDB file
        param str dir a string that specifies which XML file this function parses.
            
        return:
            the tree object it parsed initially. Reuse large object to save memory
        '''
        print("Start parsing ...")
        now = time.time()
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/' + dir)
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
                           "wikipedia": "NA",
                           "nugowiki": "NA",
                           "metagene": "NA",
                           "metlin_id": "NA",
                           "pubchem_compound_id": "NA",
                           "het_id": "NA",
                           "hmdb_id": 'NA',
                           "CAS": "NA",
                           'LIPIDMAPS':'NA'}
                # Store target tag and key value for mapping dictionary
                idtag = {"chebi_id":'chebi_id',
                         "kegg_id":'kegg_id',
                         "accession":'hmdb_id',
                         "chemspider_id":'chemspider_id',
                         "biocyc_id":'biocyc_id',
                         "cas_registry_number":"CAS",
                         "metlin_id":"metlin_id",
                         "pubchem_compound_id":"pubchem_compound_id"
                         }
                
                commonName = metabolite.find('{http://www.hmdb.ca}name').text
                
                #print(metabohmdbid)
                #find other ids for metabolite 
                for child in metabolite:
                    childtag = child.tag.replace("{http://www.hmdb.ca}", "")
                    
                    #print(childtag)
                    if childtag == "accession":
                        metabohmdbid = child.text
                        mapping['hmdb_id']  = metabohmdbid
                    # if this tag is in the id we are looking for
                    elif childtag in idtag:
                        source = idtag[childtag]
                        # if has id in the tag, append it to the list 
                        if type(mapping[source]) is not list and child.text is not None:
                            if childtag == 'chemspider_id':
                                mapping[source] = ['chemspider:'+child.text]
                            elif childtag == 'pubchem_compound_id':
                                mapping[source] =['pubchem:' + child.text]
                            elif childtag == 'chebi_id':
                                mapping[source] =['chebi:'+child.text]
                            else:
                                mapping[source] =[child.text]
                        elif type(mapping[source]) is list and child.text is not None:
                            mapping[source].append(child.text)
                #place all id information in this mapping        
                if metabohmdbid not in self.metaboliteIDDictionary:
                    self.metaboliteIDDictionary[metabohmdbid] = mapping
                if commonName is not None:
                    self.metaboliteCommonName[metabohmdbid] = commonName
                else:
                    self.metaboliteCommonName[metabohmdbid] = "NA"        
                        
        return tree 
             
                        
                    
                

    
    def getPathwaysandSynonyms(self,tree = None,dir = 'hmdb_metabolites.xml'):
        '''
        This functions finds pathways and synonyms for the metabolites and places them in:
            
            - self.metabolitesWithPathwaysDictionary
            - self.metabolitesWithSynonymsDictionary
            
        Additionally it creates a mapping between the pathwayid and the pathway name and places it in:
        
            - self.pathwayDictionary
            
        
        '''
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/'+dir)
        root = tree.getroot()
        metabohmdbid = "NA"
        #we will need to iterate through the xml tree to find the information we are looking for
        print('####### Get pathway #######')
        smpdb2 = self.getSMPDB_Category()
        for metabolite in root:
          
            #XML trees sometimes have namespace prefixes on the nodes. They need to be removed for parsing.
            #That's the point of the find and replace. We are removing the namespace string "{http://www.hmdb.ca}"
            
            metabolitetag = metabolite.tag.replace("{http://www.hmdb.ca}", "")
            # Check if this metabolite has pathway
            haspathway = False
            
            #find the accession number (metabolite id)
            if metabolitetag == "metabolite":
                accessiontag = metabolite.find('{http://www.hmdb.ca}accession')
                pathways = metabolite.find('{http://www.hmdb.ca}pathways')
                synonyms = metabolite.find('{http://www.hmdb.ca}synonyms')
                # Find accession number
                if accessiontag is not None and accessiontag.text is not None:
                    metabohmdbid = accessiontag.text
        
                    if metabohmdbid not in self.metabolitesWithSynonymsDictionary:
                        self.metabolitesWithSynonymsDictionary[metabohmdbid] =[]
                    if metabohmdbid not in self.metabolitesWithPathwaysDictionary:
                        self.metabolitesWithPathwaysDictionary[metabohmdbid] =[]
                else:
                    raise ValueError('Accession number cannot be None')
                if metabohmdbid is 'NA':
                    raise ValueError('Metabolite ID cannot be None Type')
                # find pathways
                if pathways is not None:
                    listOfPathways = []
                    for pathway in pathways:
                        pathwayNametag = pathway.find('{http://www.hmdb.ca}name')
                        smpidtag = pathway.find('{http://www.hmdb.ca}smpdb_id')
                        if pathwayNametag is not None and smpidtag is not None:
                            smpid = smpidtag.text
                            pathwayName = pathwayNametag.text
                            if smpid is not None and pathwayName is not None and smpid not in self.pathwayDictionary:
                                self.pathwayDictionary[smpid] = pathwayName
                                self.pathwayCategory[smpid] = 'NA'
                                if smpid not in self.metabolitesWithPathwaysDictionary[metabohmdbid]:
                                    self.metabolitesWithPathwaysDictionary[metabohmdbid].append(smpid)
                                    
                    
                else:
                    raise ValueError('Each metabolites tag has a pathways children')
                # find synonyms 
                if synonyms is not None:
                    for synonym in synonyms:
                        if synonym is not None and synonym.text is not None:
                            self.metabolitesWithSynonymsDictionary[metabohmdbid].append(synonym.text)
                        
                # find accession
                #find other ids for metabolite 
                '''
                for child in metabolite:
                    childtag = child.tag.replace("{http://www.hmdb.ca}", "")
                    
                    #THIS WILL BE THE KEY
                    if childtag == "accession":
                        metabohmdbid = child.text 
                        #print("Getting ..." + metabohmdbid)    
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
                                            self.pathwayCategory[smpid] = 'NA'
                                            
                                                         
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
                    '''    
        for key in self.pathwayCategory:
            if key in smpdb2:
                self.pathwayCategory[key] = 'smpdb2'
            else:
                self.pathwayCategory[key] = 'smpdb3'
        print('######### Finished for metaboliteIDDict ###########')
        print('{} items in pathwayDictionary.'.format(len(self.pathwayDictionary)))
        print('{} items in metabolitesWithSynonyms dictionary'.format(len(self.metabolitesWithSynonymsDictionary)))
        return tree                     
                           
             
    def getGenes(self,tree = None,dir = 'hmdb_metabolites.xml'):
        '''
        This function finds genes linked to metabolites and places them in:
        
            -self.MetabolitesLinkedToGenes
            
        Additionally, it links the uniprotid to the gene name and place it in:
            
            -self.geneInfoDictionary
        
        And, finally, it finds other ids for every gene and places this in:
        
            -self.geneInfoDictionary
            
        param elementTree tree parsed XML file from HMDB source file
        '''
        #key: metabolite id, value: list of gene ids
        self.metabolitesLinkedToGenes = dict()
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/'+ dir)
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
                            
                            idtag = {"uniprot_id":'Uniprot',
                                     "protein_accession":'HMDB_protein_accession',
                                     "gene_name":'common_name',
                            }
                            for key in idtag:
                                # Find all target tag in idtag.keys()
                                sourceid = protein.find('{http://www.hmdb.ca}' + key)
                                
                                # replace the name space
                                id_tag_key = sourceid.tag.replace('{http://www.hmdb.ca}','')
                                mapping_key = idtag[id_tag_key]
                                if sourceid.text is not None:
                                    mapping[mapping_key] = sourceid.text
                                    
                            
                            proteinacc = mapping['HMDB_protein_accession']
                            listofgenes.append(proteinacc)
                            self.geneInfoDictionary[proteinacc] = mapping
                        #using the uniprot id as the gene id 
            self.metabolitesLinkedToGenes[metabohmdbid] = listofgenes  
        print("Length of geneInfoDict is {}".format(str(len(self.geneInfoDictionary))))
        print('Length of metabolite-gene is {}'.format(len(self.metabolitesLinkedToGenes)))
                            
              
    def getBiofluidCellularLocationDisease(self,tree = None,dir = 'hmdb_metabolites.xml'):
        
        '''This function finds biofluid and cellular location infromation for every metabolite and places them in:
            -self.cellularLocation
            -self.biofluidLocation
            
        Additionally, a running list of all biofluid and cellular locations are kept:
        
            -self.cellular
            -self.biofluid
        
        param elementTree tree parsed XML file from HMDB source file
        
        
        '''
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/' + dir)
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
    '''
    Try to link hmdb gene to pathways in this case. Note: This part uses different 
    XML file parsed compared to previous functions.
    Key: pathway id SMP_NUM Value: Geneid HMDBP_NUM
    This function fills dictionary:
        self.PathwayLinkedToGene key: pathway id Value: HMDBP ID
        self.pathwayDictionary key pathway id Value: Name
        self.metabolitesLinkedToGenes key: HMDB ID Value: List of HMDBPID
        
    param elementTree tree parsed XML file from HMDB source file
    return:
        elementTree tree reused large object to save memory
    '''
    def getPathwaysLinkedToGene(self,tree = None):
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/hmdb_proteins.xml')
        smpdb2 = self.getSMPDB_Category()
        root = tree.getroot() 
        for protein in root:
            accession = protein.find('{http://www.hmdb.ca}accession')
            
            uniprotidtag = protein.find('{http://www.hmdb.ca}uniprot_id')
            if uniprotidtag is not None and uniprotidtag.text is not None:
                accessionnum = accession.text
                if accessionnum not in self.geneInfoDictionary:
                    mapping = list(self.geneInfoDictionary.values())[0]
                    mapping['HMDB_protein_accession'] = accessionnum
                    mapping['UniProt'] = uniprotidtag.text
                    self.geneInfoDictionary[accessionnum] = mapping
                else:
                    self.geneInfoDictionary[accessionnum]['UniProt'] = uniprotidtag.text
            for pathways in protein.iter('{http://www.hmdb.ca}pathways'):
                for pathway in pathways:
                    pathwaytag = pathway.tag.replace('{http://www.hmdb.ca}','')
                    pathwayName = pathway.find('{http://www.hmdb.ca}name').text
                    smpid = pathway.find('{http://www.hmdb.ca}smpdb_id').text
                    if smpid is not None:
                        if pathwayName is not None and smpid not in self.pathwayDictionary:
                            self.pathwayDictionary[smpid] = pathwayName
                        if smpid not in self.pathwayCategory:
                            if smpid in smpdb2:
                                self.pathwayCategory[smpid] = 'smpdb2'
                            else:
                                self.pathwayCategory[smpid] = 'smpdb3'     
                    keggid = None
                    # Pathway ID now have the kegg id
                    # Considering incorporate KEGG pathway in the future
                    for info in pathway:
                        infotag = info.tag.replace('{http://www.hmdb.ca}','')
                        if infotag == "name":
                            pathwayName = info.text
                        if infotag == 'smpdb_id':
                            if info.text is not None:
                                smpid = info.text
                                if smpid not in self.pathwayDictionary and smpid is not None:
                                    self.pathwayDictionary[smpid] = pathwayName
                                if smpid not in self.pathwaysWithGenesDictionary:
                                    self.pathwaysWithGenesDictionary[smpid] = []
                                    self.pathwaysWithGenesDictionary[smpid].append(accession.text)
                                else:
                                    if accession.text not in self.pathwaysWithGenesDictionary[smpid]:
                                        self.pathwaysWithGenesDictionary[smpid].append(accession.text)
            for metabolite in protein.find('{http://www.hmdb.ca}metabolite_associations'):
                hmdb_id = metabolite.find('{http://www.hmdb.ca}accession').text # find the metabolite id under this node
                # Use protein file to add more information to metablite-gene relations
                if hmdb_id not in self.metabolitesLinkedToGenes:
                    self.metabolitesLinkedToGenes[hmdb_id] = [accession.text]
                else:
                    if accession.text not in self.metabolitesLinkedToGenes[hmdb_id]:
                        self.metabolitesLinkedToGenes[hmdb_id].append(accession.text)
                
        print('After parsing protein file, geneInfo has {} items'.format(len(self.geneInfoDictionary)))
        print('After parsing protein file, metabolites-gene has {} items'.format(len(self.metabolitesLinkedToGenes)))
        return tree            
    def getSMPDB_Category(self):
        SMPDB2 = []
        with open('../misc/data/hmdb/SMPDB.txt','r') as f:
            for line in f:
                SMPDB2.append(line.rstrip('\n'))
        return SMPDB2
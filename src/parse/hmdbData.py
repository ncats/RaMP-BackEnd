# import xml.etree.ElementTree as ET
import urllib.request
# from xml.etree.ElementTree import ElementTree
from lxml import etree as ET
import zipfile
import time
import os
from parse.MetabolomicsData import MetabolomicsData
import pandas as pd
from rampConfig.RampConfig import RampConfig



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
        - 4) getPathwayLinkedToGene()
        - 5) WriteToFiles()
    
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
    
    def __init__(self, resConfig):
        
        super().__init__()
        # self.tree = ET.parse('../misc/data/hmdb/hmdb_metabolites.xml')
        # self.tree2 = ET.parse('../misc/data/hmdb/hmdb_proteins.xml')
        ####DICTIONARIES IN COMMON WITH OTHER CLASSES######################################

        self.resourceConfig = resConfig
        
        # common name dictionary Key: HMDB ID Value: Common Name
        self.metaboliteCommonName = dict()
        # pathway dictionary. Key: hsaID for pathway, Value: pathway name
        self.pathwayDictionary = dict()
        # pathway id mapping Key SMP ID Value: Kegg id
        self.SMPToKegg = dict()
        # hsaID for pathway, value: pathway category (all will be "NA" for hmdb)
        self.pathwayCategory = dict()
        
        # key: metabolite id , value: list of pathway id
        self.metabolitesWithPathwaysDictionary = dict()
        
        # key: metabollite id, value: list of synonyms 
        self.metabolitesWithSynonymsDictionary = dict()
        
        # key: metabolite id, value: mapping of other ids
        self.metaboliteIDDictionary = dict()
        self.metaInchi = dict()
        # key: pathway SMP id, value: list of gene HMDBP id
        self.pathwaysWithGenesDictionary = dict()      
        # key: pathway id, value: list of metabolites id with this pathways.
        self.pathwaysWithMetabolitesDictionary = dict()
        # key: gene id, value: list of FOUR gene identifiers 
        self.geneInfoDictionary = dict()
        
        # only not empty when a catalyzed class exists 
        # key: matabole, value: list of genes
        self.metabolitesLinkedToGenes = dict()
        
        # holds protein ids to their type designation
        # primarily used to populate the catalyzed resource.
        self.protein2type = dict()
        
        # self.inchiDict = dict()
        ###################################################################
        
        # stays empty for this class
        self.pathwayOntology = dict()
        
        # key: metabolite id, value: list of metabolite locations
        self.biofluidLocation = dict()
        
        # key: biofluid location, value: the string "placeholder"
        self.biofluid = dict()
        
        # key: metabolite id, value: list of cellular locations
        self.cellularLocation = dict()
        
        # key: cellular location, value: the string "placeholder"
        self.cellular = dict()
        
        # key: metaboliteID, value: exo/endo/Drug metabolite
        self.exoEndoDictionary = dict()
        
        # key: origin, value: exo/endo/ Drug ,etc.
        self.exoEndo = dict()
        # key: metaboliteID, value: tissue_locations
        self.tissueLocation = dict()

        # key: metaboliteId, value: organ_location        
        self.organLocation = dict()
        
        # key: tissue location, value : "placeholder"
        self.tissue = dict()
        
        self.idDictForMetabolite = dict()
        # key: source ID e.g. HMDBID, value: dictionary that has key of sub,class,super class
        # value as the class name
        self.metaboliteClass = dict()
        
        # holds industrial application, like 'Drug', 'Self Care Product'
        self.metaboliteApplication = dict()
    
        # holds industrial application, like 'Drug', 'Self Care Product'
        self.healthCondition = dict()
        
        # holds HMDB status for each metabolite
        self.metStatus = dict()
        
    def getEverything(self, writeToFile=False):
        '''
        Run all the function to get everything from hmdb source
        - param bool writeToFile if true, write all dictionaries to misc/output/hmdb/
        '''
        self.getDatabaseFiles()
        tree = self.getMetaboliteOtherIDs()

        self.getPathwaysandSynonyms(tree)
        
        self.getGenes(tree)

        self.getOntology(tree)
        self.getPathwaysLinkedToGene()
        self.getMetabolitesClasses(tree)
        self.getStatus(tree)
        
        # wth just adding a type to met -> prot -> type, by adding <tab>protein_type. :)
        self.annealProteinTypeToMet2ProtDict()
        
        
        if writeToFile:
            self.write_myself_files('hmdb')
        
    def getDatabaseFiles(self):
        '''
		This function gets the files that make up hmdb and places them into the hmdb folder. 
		
        '''
        # Define file name and url for downloading file
#         file_metabolites = "hmdb_metabolites.zip"
#         file_proteins = "hmdb_proteins.zip"
#         download_url = "https://hmdb.ca/system/downloads/current/"
        
        metConfig = self.resourceConfig.getConfig("hmdb_met")
        proteinConfig = self.resourceConfig.getConfig("hmdb_gene")
                
        file_metabolites = metConfig.sourceFileName
        mets_url = metConfig.sourceURL
        
        file_proteins = proteinConfig.sourceFileName
        proteins_url = proteinConfig.sourceURL
        
        localDir = metConfig.localDir
        
        print("getting HMDB files, localDir = " + localDir)
        
        # check if this path exists
        self.check_path(localDir)
        if file_metabolites.replace('.zip', '.xml') not in os.listdir(localDir) or file_proteins.replace('.zip', '.xml') not in os.listdir(localDir):
            # Download files from given url and path
            print('####### Downloading HMDB source file #######')
            self.download_files(mets_url, localDir + file_metabolites)
            self.download_files(proteins_url, localDir + file_proteins)
            # Open zip file if file are downloaded.
            with zipfile.ZipFile(localDir + file_metabolites, "r") as zip_ref:
                zip_ref.extractall(localDir)
                
            with zipfile.ZipFile(localDir + file_proteins, "r") as zip_ref:
                zip_ref.extractall(localDir)
        else:
            print('HMDB source files are ready ...')
                                    
    def getMetaboliteOtherIDs(self, tree=None, dir='hmdb_metabolites.xml'):     
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
            
        # we will need to iterate through the xml tree to find the information we are looking for
        for metabolite in root:
        # XML trees sometimes have namespace prefixes on the nodes. They need to be removed for parsing.
            # That's the point of the find and replace. We are removing the namespace string "{http://www.hmdb.ca}"

            metabolitetag = metabolite.tag.replace("{http://www.hmdb.ca}", "")
            
            # find the accession number (metabolite id)
            if metabolitetag == "metabolite":
                # #here are some of the things we will be looking for in the xml tree
                mapping = {
                           "smiles":'NA',
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
                           "hmdb_id": 'NA',
                           "CAS": "NA",
                           'LIPIDMAPS':'NA',
                           'hmdb_secondary_id':'NA'}
                # Store target tag and key value for mapping dictionary
                idtag = {"chebi_id":'chebi_id',
                         "kegg_id":'kegg_id',
                         "accession":'hmdb_id',
                         "chemspider_id":'chemspider_id',
                         "biocyc_id":'biocyc_id',
                         "cas_registry_number":"CAS",
                         "metlin_id":"metlin_id",
                         "pubchem_compound_id":"pubchem_compound_id",
                         'smiles':'smiles',
                         'LIPIDMAPS':'LIPIDMPAS'
                         }
                
                commonName = metabolite.find('{http://www.hmdb.ca}name').text
                commonName = commonName.strip()

                # find other ids for metabolite
                # prefix is the id we collect and would like to store it in RamP 
                prefix = {'chemspider_id':'chemspider:',
                          'pubchem_compound_id':'pubchem:',
                          'chebi_id':'chebi:',
                          'CAS':'CAS:',
                          'kegg_id':'kegg:',
                          'smiles':'smiles:'
                          }
                for child in metabolite:
                    childtag = child.tag.replace("{http://www.hmdb.ca}", "")
                    # print(childtag)
                    if childtag == "accession":
                        metabohmdbid = 'hmdb:' + child.text.strip()
                        mapping['hmdb_id'] = metabohmdbid
                        
                    # capture secondary ids    
                    elif childtag == "secondary_accessions":
                        children = child.findall('{http://www.hmdb.ca}accession')
                        if children is not None:
                            for kid in children:
                                if type(mapping['hmdb_secondary_id']) is not list:
                                    mapping['hmdb_secondary_id'] = list()                               
                                mapping['hmdb_secondary_id'].append('hmdb:' + kid.text.strip())
                            
                    # if this tag is in the id we are looking for
                    elif childtag in idtag:
                        source = idtag[childtag]
                        # if has id in the tag, append it to the list 
                        if type(mapping[source]) is not list and child.text is not None:
                            if source in prefix:
                                mapping[source] = [prefix[source] + child.text.strip()]
                        elif type(mapping[source]) is list and child.text is not None:
                            if source in prefix:
                                mapping[source].append(prefix[source] + child.text.strip())
                # place all id information in this mapping        
                if metabohmdbid not in self.metaboliteIDDictionary:
                    self.metaboliteIDDictionary[metabohmdbid] = mapping
                if commonName is not None:
                    self.metaboliteCommonName[metabohmdbid] = commonName

                    # append the commonName as another synonym
                    if metabohmdbid not in self.metabolitesWithSynonymsDictionary:
                        self.metabolitesWithSynonymsDictionary[metabohmdbid] = []                    
                    if commonName not in self.metabolitesWithSynonymsDictionary[metabohmdbid]:
                        self.metabolitesWithSynonymsDictionary[metabohmdbid].append(commonName)
                else:
                    self.metaboliteCommonName[metabohmdbid] = "NA"        
                        
        return tree 
    
    def getPathwaysandSynonyms(self, tree=None, dir='hmdb_metabolites.xml'):
        '''
        This functions finds pathways and synonyms for the metabolites and places them in:
            
            - self.metabolitesWithPathwaysDictionary
            - self.metabolitesWithSynonymsDictionary
            
        Additionally it creates a mapping between the pathwayid and the pathway name and places it in:
        
            - self.pathwayDictionary
            
        
        '''
        countGlobal = 0
        countMicrobe = 0
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/' + dir)
        root = tree.getroot()
        metabohmdbid = "NA"
        # we will need to iterate through the xml tree to find the information we are looking for
        print('####### Get pathway #######')
        # hmdbinchiKeyFile = open("../misc/sql/" + "hmdbinchiKeyFile.sql", 'wb')
        lipidCount = 0
        for metabolite in root:
            countGlobal += 1
            # print(type(metabolite))
            # print('Metaboliteed', metabolite)
            # XML trees sometimes have namespace prefixes on the nodes. They need to be removed for parsing.
            # That's the point of the find and replace. We are removing the namespace string "{http://www.hmdb.ca}"
            # print(metabohmdbid)
            metabolitetag = metabolite.tag.replace("{http://www.hmdb.ca}", "")
            # Check if this metabolite has pathway
            haspathway = False
            # print('metatag', metabolite)
            # find the accession number (metabolite id)

            if metabolitetag == "metabolite":
                accessiontag = metabolite.find('{http://www.hmdb.ca}accession')
                pathways = metabolite.find('{http://www.hmdb.ca}pathways')
                synonyms = metabolite.find('{http://www.hmdb.ca}synonyms')
                biological_properties = metabolite.find('{http://www.hmdb.ca}biological_properties')
                inchikey = metabolite.find('{http://www.hmdb.ca}inchikey')
                inchi = metabolite.find('{http://www.hmdb.ca}inchi')
                smiles = metabolite.find('{http://www.hmdb.ca}smiles')

                tax = metabolite.find('{http://www.hmdb.ca}taxonomy')
                if tax is not None:
                    superClass = tax.find('{http://www.hmdb.ca}super_class')
                    if superClass is not None:
                        if superClass.text == "Lipids and lipid-like molecules":
                            lipidCount = lipidCount + 1
                        # else:
                        #    print("no true*********")

                # end of experiment

                # Find accession number
                if accessiontag is not None and accessiontag.text is not None:
                    metabohmdbid = 'hmdb:' + accessiontag.text
                    inchiID = inchi.text
                    if inchiID is not None:
                        self.metaInchi[metabohmdbid] = inchiID
                    
                    if metabohmdbid not in self.metabolitesWithSynonymsDictionary:
                        self.metabolitesWithSynonymsDictionary[metabohmdbid] = []
                    if metabohmdbid not in self.metabolitesWithPathwaysDictionary:
                        self.metabolitesWithPathwaysDictionary[metabohmdbid] = []
                else:
                    raise ValueError('Accession number cannot be None')
                if metabohmdbid is 'NA':
                    raise ValueError('Metabolite ID cannot be None Type')
                # find pathways
                # print('Pathwaycsd', pathways)

                if pathways is not None:
                    listOfPathways = []
                    for pathway in pathways:
                        # Find pathway name and smp id
                        # kegg id is not collected
                        pathwayNametag = pathway.find('{http://www.hmdb.ca}name')
                        smpidtag = pathway.find('{http://www.hmdb.ca}smpdb_id')
                        # print('++++++++SMPD', smpidtag)
                        if pathwayNametag is not None and smpidtag is not None:
                            smpid = smpidtag.text
                            pathwayName = pathwayNametag.text
                            if smpid is not None and pathwayName is not None and smpid not in self.pathwayDictionary:
                                self.pathwayDictionary[smpid] = pathwayName
                                self.pathwayCategory[smpid] = 'NA'
                                # add pathways to metabolites With pathway dictionary
                                
                            # JCB: The above IF adds pathway to the pathway dictionary if not in there...
                            # JCB: But this addition to the M2P dictionary won't happen if the smp1d isn't in the dict.?
                            # JCB: Consider moving this if statement down, this should be run if the pathway is in the dictionary. 
                            # JCB: Added 'if smpid is not None:' to bring down to add mapping for all smp-path ids
                            if smpid is not None:
                                if smpid not in self.metabolitesWithPathwaysDictionary[metabohmdbid]:
                                    self.metabolitesWithPathwaysDictionary[metabohmdbid].append(smpid)

# new HMDB taxonomy for pathways - edit: Manju

                elif biological_properties is not None:
                    pathways = biological_properties.find('{http://www.hmdb.ca}pathways')
                    # print("for biological_properties pathways", accessiontag.text)
                    listOfPathways = []
                    for pathway in pathways:
                        # Find pathway name and smp id
                        # kegg id is not collected
                        pathwayNametag = pathway.find('{http://www.hmdb.ca}name')
                        smpidtag = pathway.find('{http://www.hmdb.ca}smpdb_id')
                        keggidtag = pathway.find('{http://www.hmdb.ca}kegg_map_id')
                        # print('++++++++SMPD', smpidtag)
                        if pathwayNametag is not None and smpidtag is not None:
                            smpid = smpidtag.text
                            pathwayName = pathwayNametag.text
                            if smpid is not None and pathwayName is not None and smpid not in self.pathwayDictionary:
                                self.pathwayDictionary[smpid] = pathwayName
                                self.pathwayCategory[smpid] = 'NA'
                                # add pathways to metabolites With pathway dictionary
                                
                            # JCB: The above IF addes pathway to the pathway dictionary if not in there...
                            # JCB: But this addition to the M2P dictionary won't happen if the smp1d isn't in the dict.?
                            # JCB: Consider moving this if statement down, this should be run if the pathway is in the dictionary.     
                            # JCB: Added 'if smpid is not None:' to bring down to add mapping for all smp-path ids
                            if smpid is not None:
                                if smpid not in self.metabolitesWithPathwaysDictionary[metabohmdbid]:
                                    self.metabolitesWithPathwaysDictionary[metabohmdbid].append(smpid)

                        if pathwayNametag is not None and keggidtag is not None:
                            keggid = keggidtag.text
                            pathwayName = pathwayNametag.text
                            if keggid is not None and pathwayName is not None and keggid not in self.pathwayDictionary:
                                self.pathwayDictionary[keggid] = pathwayName
                                self.pathwayCategory[keggid] = 'kegg'
                                # add pathways to metabolites With pathway dictionary
                                
                            # JCB: The above IF addes pathway to the pathway dictionary if not in there...
                            # JCB: But this addition to the M2P dictionary won't happen if the smp1d isn't in the dict.?
                            # JCB: Consider moving this if statement down, this should be run if the pathway is in the dictionary.                                     
                            # JCB: Added 'if keggid is not None:' to bring down to add all keggmaps
                            if keggid is not None:
                                if keggid not in self.metabolitesWithPathwaysDictionary[metabohmdbid]:
                                    self.metabolitesWithPathwaysDictionary[metabohmdbid].append(keggid)

                else:
                    # print("Pathway None..", accessiontag.text)
                    raise ValueError('Each metabolites tag has a pathways children')
                # find synonyms 
                if synonyms is not None:
                    for synonym in synonyms:
                        if synonym is not None and synonym.text is not None:
                            if synonym.text not in self.metabolitesWithSynonymsDictionary[metabohmdbid]:
                                # add synonym if not already there.
                                self.metabolitesWithSynonymsDictionary[metabohmdbid].append(synonym.text)
                                
        print("count global: ", countGlobal, " count micro ", countMicrobe)
        print("lipid super class", lipidCount)
        for key in self.pathwayCategory:
            if self.pathwayCategory[key] is 'NA':
                self.pathwayCategory[key] = 'smpdb3'
        # hmdbinchiKeyFile.close()
        print('{} items in pathwayDictionary.'.format(len(self.pathwayDictionary)))
        print('{} items in metabolitesWithSynonyms dictionary'.format(len(self.metabolitesWithSynonymsDictionary)))
        return tree

    def getGenes(self, tree=None, dir='hmdb_metabolites.xml'):
        '''
        This function finds genes linked to metabolites and places them in:
        
            -self.MetabolitesLinkedToGenes
            
        Additionally, it links the uniprotid to the gene name and place it in:
            
            -self.geneInfoDictionary
        
        And, finally, it finds other ids for every gene and places this in:
        
            -self.geneInfoDictionary
            
        param elementTree tree parsed XML file from HMDB source file
        '''
        # key: metabolite id, value: list of gene ids
        self.metabolitesLinkedToGenes = dict()
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/' + dir)
        root = tree.getroot()
        metabohmdbid = "not yet found"
        # we will need to iterate through the xml tree to find the information we are looking for
        for metabolite in root:
          
            # XML trees sometimes have namespace prefixes on the nodes. They need to be removed for parsing.
            # That's the point of the find and replace. We are removing the namespace string "{http://www.hmdb.ca}"

            metabolitetag = metabolite.tag.replace("{http://www.hmdb.ca}", "")
            listofgenes = []
            prefix = {"UniProt":'uniprot:',
                      'HMDB_protein_accession':'hmdb:',
                      'gene_name':'gene_symbol'}
            if metabolitetag == "metabolite":
                # find other ids for metabolite 
                for child in metabolite:
                    childtag = child.tag.replace("{http://www.hmdb.ca}", "")

                    # THIS WILL BE THE KEY
                    if childtag == "accession":
                        metabohmdbid = 'hmdb:' + child.text 
                        
                    if childtag == "protein_associations":
                        
                        for protein in child:
                            # print("Protein*******")
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
                                       'Enzyme Nomenclature': 'NA',
                                       'gene_name': 'NA'}
                            
                            idtag = {"uniprot_id":'UniProt',
                                     "protein_accession":'HMDB_protein_accession',
                                     "gene_name":'common_name',
                            }
                            for key in idtag:
                                # print("idTag*******")
                                # Find all target tag in idtag.keys()
                                sourceid = protein.find('{http://www.hmdb.ca}' + key)
                                # print(sourceid)
                                # replace the name space
                                id_tag_key = sourceid.tag.replace('{http://www.hmdb.ca}', '')
                                mapping_key = idtag[id_tag_key]
                                # print("******",id_tag_key, mapping_key)
                                # print(metabohmdbid)

                                if sourceid.text is not None:
                                    if mapping_key is not 'common_name':
                                        # if metabohmdbid == "hmdb:HMDB0004970":
                                         #   print("for HMDB0004970:", [prefix[mapping_key] + sourceid.text], mapping_key)
                                        mapping[mapping_key] = [prefix[mapping_key] + sourceid.text]
                                    else:
                                        if mapping_key == 'gene_name' or mapping_key == 'common_name':
                                            mapping[mapping_key] = "gene_symbol:" + sourceid.text
                                            mapping['common_name'] = "gene_symbol:" + sourceid.text
                            
                            proteinacc = mapping['HMDB_protein_accession'][0]
                            # print(proteinacc)
                            listofgenes.append(proteinacc)
                            self.geneInfoDictionary[proteinacc] = mapping
                        # using the uniprot id as the gene id
            # print("+++", metabohmdbid)
            self.metabolitesLinkedToGenes[metabohmdbid] = listofgenes  
        print("Length of geneInfoDict is {}".format(str(len(self.geneInfoDictionary))))
        print('Length of metabolite-gene is {}'.format(len(self.metabolitesLinkedToGenes)))

    def getOntology(self, tree=None, dir='hmdb_metabolites.xml'):
        # get disposition
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/' + dir)
        
        root = tree.getroot()
        
        for metabolite in root:
            
            metId = self.getMetaboliteIdTag(metabolite)
            
            if metId is None:
                continue
                        
            ontology = metabolite.find('{http://www.hmdb.ca}ontology')
                    
            if ontology is not None:
                self.parseSource(ontology, metId)
                self.parseTissue(ontology, metId)
                self.parseBiofluid(ontology, metId)
                self.parseCellLocation(ontology, metId)
                self.parseApplication(ontology, metId)
                self.parsePhysiology(ontology, metId)


    def parseSource(self, ontology, metId):
        
        endoExoList = [
            'Endogenous',
            'Exogenous',
            'Microbial',
            'Drug or steroid metabolite',
            'Drug metabolite',
            'Toxin/Pollutant',
            'Food',
            'Cosmetic',
            'Plant',
            'Drug',
            'Environmental',
            'Saccharomyces cerevisiae',
            'Synthetic',
            'Animal',
            'Fungi',
            'Microbe']
        
        sourceNode = self.findSingleDecendentNodeTerm(ontology, "Source")
        
        if sourceNode is not None:      
            for dec in sourceNode.iter('{http://www.hmdb.ca}descendant'):
                for termNode in dec.iter('{http://www.hmdb.ca}term'):            
                    term = termNode.text.strip()                    
                    if(term in endoExoList):                       
                        if metId in self.exoEndoDictionary:
                            if(term not in self.exoEndoDictionary[metId]):
                                self.exoEndoDictionary[metId].append(term)
                        else:
                            self.exoEndoDictionary[metId] = [term]
        
    def parseTissue(self, ontology, metId):

        keyTerms = ["Tissue and substructures", "Organ and components"]
        
        tAndSNode = self.findSingleDecendentNodeTerm(ontology, "Tissue and substructures")
        oAndCNode = self.findSingleDecendentNodeTerm(ontology, "Organ and components")
        
        if tAndSNode is not None:
            for tissueDec in tAndSNode.iter('{http://www.hmdb.ca}descendant'):
                termNode = tissueDec.find('{http://www.hmdb.ca}term')
                if termNode is not None:
                    term = termNode.text.strip()

                    if term in keyTerms:
                        continue

                    if metId in self.tissueLocation:
                        if term not in self.tissueLocation[metId]:
                            self.tissueLocation[metId].append(term)
                    else:
                        self.tissueLocation[metId] = [term]
                        
        if oAndCNode is not None:
            for organComp in oAndCNode.iter('{http://www.hmdb.ca}descendant'):
                termNode = organComp.find('{http://www.hmdb.ca}term')
                if termNode is not None:
                    term = termNode.text.strip()
                    
                    if term in keyTerms:
                        continue
             
                    if metId in self.organLocation:
                        if term not in self.organLocation[metId]:
                            self.organLocation[metId].append(term)
                    else:
                        self.organLocation[metId] = [term]               
                           
    def parseBiofluid(self, ontology, metId):
 
        bfTerm = "Biofluid and excreta"
        
        biofluidsNode = self.findSingleDecendentNodeTerm(ontology, bfTerm)
              
        if biofluidsNode is not None:
            for bioNodeDec in biofluidsNode.iter('{http://www.hmdb.ca}descendant'):
                termNode = bioNodeDec.find('{http://www.hmdb.ca}term')
                if termNode is not None:
                    term = termNode.text.strip()
                    
                    if term == bfTerm:
                        continue
                    
                    if metId in self.biofluidLocation:
                        if term not in self.biofluidLocation[metId]:
                            self.biofluidLocation[metId].append(term)
                    else:
                        self.biofluidLocation[metId] = [term]
        
    def parseCellLocation(self, ontology, metId):
        
        cellTerm = "Subcellular"
        
        cellNode = self.findSingleDecendentNodeTerm(ontology, cellTerm)
              
        if cellNode is not None:
            for cellNodeDec in cellNode.iter('{http://www.hmdb.ca}descendant'):
                termNode = cellNodeDec.find('{http://www.hmdb.ca}term')
                if termNode is not None:
                    term = termNode.text.strip()
                    
                    if term == cellTerm:
                        continue
                    
                    if metId in self.cellularLocation:
                        if term not in self.cellularLocation[metId]:
                            self.cellularLocation[metId].append(term)
                    else:
                        self.cellularLocation[metId] = [term]
        
    def parseApplication(self, ontology, metId):
    
        appTerm = "Industrial application"
                
        appNode = self.findSingleDecendentNodeTerm(ontology, appTerm)
        
        if appNode is not None:
            for appNodeDec in appNode.iter('{http://www.hmdb.ca}descendant'):
                termNode = appNodeDec.find('{http://www.hmdb.ca}term')
                if termNode is not None:
                    term = termNode.text.strip()
                    
                    if term == appTerm:
                        continue
                    
                    if metId in self.metaboliteApplication:
                        if term not in self.metaboliteApplication[metId]:
                            self.metaboliteApplication[metId].append(term)
                    else:
                        self.metaboliteApplication[metId] = [term]
 
    def parsePhysiology(self, ontology, metId):
    
        physTerm = "Health condition"
                
        physNode = self.findSingleDecendentNodeTerm(ontology, physTerm)
                
        if physNode is not None:
            for physNodeDec in physNode.iter('{http://www.hmdb.ca}descendant'):
                termNode = physNodeDec.find('{http://www.hmdb.ca}term')
                if termNode is not None:
                    term = termNode.text.strip()
                    
                    if term == physTerm:
                        continue
                    
                    if metId in self.healthCondition:
                        if term not in self.healthCondition[metId]:
                            self.healthCondition[metId].append(term)
                    else:
                        self.healthCondition[metId] = [term]       
    
    def findSingleDecendentNodeTerm(self, parentNode, queryTerm):
        foundit = False
   
        decs = parentNode.iter('{http://www.hmdb.ca}descendant')
        for dec in decs:            
            termNode = dec.find('{http://www.hmdb.ca}term')
            if termNode.text == queryTerm:
                return dec
        return None  
            
    def getMetaboliteIdTag(self, parentNode):
        child = parentNode.find('{http://www.hmdb.ca}accession')
        metId = None         
        if child is not None:
            metId = 'hmdb:' + child.text
        return metId

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

    def getPathwaysLinkedToGene(self, tree=None):
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/hmdb_proteins.xml')
        root = tree.getroot()
         
        for protein in root:
            accession = protein.find('{http://www.hmdb.ca}accession')
            pType = protein.find('{http://www.hmdb.ca}protein_type')        
            uniprotidtag = protein.find('{http://www.hmdb.ca}uniprot_id')
            genenametag = protein.find('{http://www.hmdb.ca}gene_name')
            id_mapping = {
                       'kegg': 'NA',
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
                       'Enzyme Nomenclature': 'NA',
                       'gene_name':'NA'}
            if uniprotidtag is not None and uniprotidtag.text is not None:
                # prepend hmdb: to the HMDBP ID
                accessionnum = 'hmdb:' + accession.text
                uniprotid = 'uniprot:' + uniprotidtag.text
                proteinType = pType.text
                self.protein2type[accessionnum] = proteinType
                
                if genenametag is not None and genenametag.text is not None:
                    genename = 'gene_symbol:' + genenametag.text
                else:
                    genename = None
                if accessionnum not in self.geneInfoDictionary:
                    mapping = id_mapping.copy()
                    mapping['HMDB_protein_accession'] = [accessionnum]
                    mapping['UniProt'] = [uniprotid]
                    if genename is not None:
                        mapping['gene_name'] = [genename]
                    
                    self.geneInfoDictionary[accessionnum] = mapping
                else:
                    if self.geneInfoDictionary[accessionnum]['UniProt'] is not 'NA' and\
                     uniprotid not in self.geneInfoDictionary[accessionnum]['UniProt']:
                        self.geneInfoDictionary[accessionnum]['UniProt'].append(uniprotid)
                    else:
                        self.geneInfoDictionary[accessionnum]['UniProt'] = [uniprotid]
                    
                    #handle gene name / symbol
                    if self.geneInfoDictionary[accessionnum]['gene_name'] is not 'NA' and\
                     genename not in self.geneInfoDictionary[accessionnum]['gene_name']:
                        self.geneInfoDictionary[accessionnum]['gene_name'].append(genename)
                    else:
                        self.geneInfoDictionary[accessionnum]['gene_name'] = [genename]
                        
                        
            for pathways in protein.iter('{http://www.hmdb.ca}pathways'):
                for pathway in pathways:
                    pathwaytag = pathway.tag.replace('{http://www.hmdb.ca}', '')
                    pathwayName = pathway.find('{http://www.hmdb.ca}name').text
                    smpid = pathway.find('{http://www.hmdb.ca}smpdb_id').text
                    keggidtag = pathway.find('{http://www.hmdb.ca}kegg_map_id')
                    
                    if smpid is not None:
                        if pathwayName is not None and smpid not in self.pathwayDictionary:
                            self.pathwayDictionary[smpid] = pathwayName
                        if smpid not in self.pathwayCategory:
                            self.pathwayCategory[smpid] = 'smpdb3'
                    keggid = None

                    if keggidtag is not None:
                        keggid = keggidtag.text
                        if keggid is not None:
                            if pathwayName is not None and keggid not in self.pathwayDictionary:
                                self.pathwayDictionary[keggid] = pathwayName
                            if keggid not in self.pathwayCategory:
                                self.pathwayCategory[keggid] = 'kegg'
                    
                    # Pathway ID now have the kegg id
                    # Considering incorporate KEGG pathway in the future
                    for info in pathway:
                        infotag = info.tag.replace('{http://www.hmdb.ca}', '')
                        if infotag == "name":
                            pathwayName = info.text
                        if infotag == 'smpdb_id':
                            if info.text is not None:
                                smpid = info.text
                                if smpid not in self.pathwayDictionary and smpid is not None:
                                    self.pathwayDictionary[smpid] = pathwayName
                                if smpid not in self.pathwaysWithGenesDictionary:
                                    self.pathwaysWithGenesDictionary[smpid] = []
                                    self.pathwaysWithGenesDictionary[smpid].append('hmdb:' + accession.text)
                                else:
                                    if accession.text not in self.pathwaysWithGenesDictionary[smpid]:
                                        self.pathwaysWithGenesDictionary[smpid].append('hmdb:' + accession.text)

                        if infotag == 'kegg_map_id':
                            if info.text is not None:
                                keggid = info.text
                                if keggid not in self.pathwayDictionary and keggid is not None:
                                    self.pathwayDictionary[keggid] = pathwayName
                                if keggid not in self.pathwaysWithGenesDictionary:
                                    self.pathwaysWithGenesDictionary[keggid] = []
                                    self.pathwaysWithGenesDictionary[keggid].append('hmdb:' + accession.text)
                                else:
                                    if accession.text not in self.pathwaysWithGenesDictionary[keggid]:
                                        self.pathwaysWithGenesDictionary[keggid].append('hmdb:' + accession.text)
                                if keggid not in self.pathwayCategory:
                                    self.pathwayCategory[keggid] = 'kegg'        
                                        
            for metabolite in protein.find('{http://www.hmdb.ca}metabolite_associations'):
                hmdb_id = 'hmdb:' + metabolite.find('{http://www.hmdb.ca}accession').text  # find the metabolite id under this node
                # Use protein file to add more information to metablite-gene relations
                if hmdb_id not in self.metabolitesLinkedToGenes:
                    self.metabolitesLinkedToGenes[hmdb_id] = ['hmdb:' + accession.text]
                else:
                    if accession.text not in self.metabolitesLinkedToGenes[hmdb_id]:
                        self.metabolitesLinkedToGenes[hmdb_id].append('hmdb:' + accession.text)
                
        print('After parsing protein file, geneInfo has {} items'.format(len(self.geneInfoDictionary)))
        print('After parsing protein file, metabolites-gene has {} items'.format(len(self.metabolitesLinkedToGenes)))
        return tree  

    def getMetabolitesClasses(self, tree=None, file='hmdb_metabolites.xml'):
        # if the source file is not parsed yet
        cols = ['hmdb_id', 'super_class', 'class', 'sub_class']
        result = pd.DataFrame(columns=cols)
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/' + file)
        root = tree.getroot()
        prefix = '{http://www.hmdb.ca}'
        i = 0
        for metabolite in root.findall(prefix + 'metabolite'):
            # print(metabolite.tag)
            hmdbid = metabolite.find(prefix + 'accession')
            taxonomy = metabolite.find(prefix + 'taxonomy')
            metabolites_class = {'ClassyFire_super_class':'NA',
                                 'ClassyFire_class':'NA',
                                 'ClassyFire_sub_class':'NA'}
            # if i % 1000 == 0:
            #    print('{} metabolites parsed'.format(i))
            if taxonomy is not None:
                super_clas = taxonomy.find(prefix + 'super_class')
                clas = taxonomy.find(prefix + 'class')
                sub_clas = taxonomy.find(prefix + 'sub_class')
                if super_clas is not None and super_clas is not None:
                    metabolites_class['ClassyFire_super_class'] = super_clas.text
                if clas is not None and clas.text is not None:
                    metabolites_class['ClassyFire_class'] = clas.text
                if sub_clas is not None and sub_clas.text is not None:
                    metabolites_class['ClassyFire_sub_class'] = sub_clas.text
            i = i + 1
            self.metaboliteClass['hmdb:' + hmdbid.text] = metabolites_class
            # print('metabolite {} has super class {} class {} subclass {}'\
            #      .format(hmdbid.text,super_clas.text,clas.text,sub_clas.text))
            '''
            result.loc[i,['hmdb_id','super_class','class','sub_class']] = [hmdbid.text,
                                                                           metabolites_class['super_class'],
                                                                           metabolites_class['class'],
                                                                           metabolites_class['sub_class']]
            i = i + 1
        
        result.to_csv('../misc/output/metabolites_class.csv')
            
        '''

    def getStatus(self, tree=None, file='hmdb_metabolites.xml'):
        if tree is None:
            tree = ET.parse('../misc/data/hmdb/' + file)
        root = tree.getroot()
        prefix = '{http://www.hmdb.ca}'

        for metabolite in root.findall(prefix + 'metabolite'):
            hmdbid = metabolite.find(prefix + 'accession')
            status = metabolite.find(prefix + 'status')
            
            if hmdbid != None and status != None:
                self.metStatus['hmdb:' + hmdbid.text] = status.text

    def annealProteinTypeToMet2ProtDict(self):
        for met in self.metabolitesLinkedToGenes:
            prots = self.metabolitesLinkedToGenes[met]
            protDict = dict()
            for prot in prots:                
                ptype = self.protein2type.get(prot, "Unknown")
                protDict[prot] = ptype
            self.metabolitesLinkedToGenes[met] = protDict

   

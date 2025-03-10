'''
Created on Nov 16, 2020

@author: braistedjc
'''
import csv
import os
from os import path
from os.path import exists
from rampEntity.Gene import Gene
from rampEntity.GeneList import GeneList
from rampEntity.Metabolite import Metabolite
from rampEntity.MetaboliteList import MetaboliteList
from rampEntity.Pathway import Pathway
from rampEntity.PathwayList import PathwayList
from rampEntity.Ontology import Ontology
from rampEntity.OntologyList import OntologyList
from chemprop.ChemWrangler import ChemWrangler
from rampEntity.Molecule import Molecule
from rampEntity.RheaReaction import RheaReaction
from rampEntity.Protein import Protein
from rampEntity.RheaCompound import RheaCompound


from pathlib import Path


import pandas as pd
import numpy as np
from pandas.io.common import file_path_to_url
from pandas.api.types import is_string_dtype
from pandas.io.html import _remove_whitespace
#from collections import defaultdict

import pubchempy as pcp
import time
from unittest.mock import inplace

class EntityBuilder(object):
    '''
    The EntityBuilder class is used to populate an harmonize objects from the rampEntity
    module that include Gene, Metabolite, Pathway, Molecule (chem properties) and related lists.
    The input data is the file set parsed from data sources into intermediate files.
    The output from the class is a collection of files that are formatted to upload to the database.
    The 'fullBuild' method runs through the process if building, harmonizing and exporting the data associated.
    During the process, all ramp data is held in one data structure for harmonization and review.
    
    Methods fall into three broad categories:
        Construction methods that build entities and relationships
        Write methods that dump the structure data into output files for database loading
        Utility methods use to evaluate and report on cases of improper mappings.
    '''

    __rampMetStartId = 0
    __rampGeneStartId = 0
    __rampPathStartId = 0
    __rampOntStartId = 0
    __rampRxnStartId = 0
    

    def __init__(self, resourceConfig):
        '''
        Constructor
        '''
        
        # config for data sources
        self.resConfig = resourceConfig
        
        # The full metabolite list object
        self.metaboliteList = MetaboliteList()
        
        # The gene list object
        self.geneList = GeneList()
        
        # The pathway list object. Note that pathway entities are also linked to metabolites and genes
        self.pathList = PathwayList()
        
        # ontology list
        self.ontologyList = OntologyList()
        
        # reaction dictionary
        self.reactionDict = dict()
        
        # List of DataSource objects. These hold data source configuration.
        self.sourceList = []
        
        # Note: The following populates reactome and wikipathway sources and appends to the default hmdb source
        # This data source list will eventually be populated by config file
        self.source = DataSource()        
        self.sourceList.append(self.source)

        if not os.path.exists(self.source.exportPath):
            os.makedirs(self.source.exportPath)

        self.dataSource2 = DataSource()
        self.dataSource2.sourceName = 'reactome'
        self.dataSource2.filePrefix = 'reactome'
        self.dataSource2.haveChemClassInfo = False
        self.dataSource2.sourceLocPath = '../misc/output/reactome'
         
        self.sourceList.append(self.dataSource2)
         
        self.dataSource3 = DataSource()        
        self.dataSource3.sourceName = 'wiki'
        self.dataSource3.filePrefix = 'wikipathwayRDF'
        self.dataSource3.haveChemClassInfo = False
        self.dataSource3.sourceLocPath = '../misc/output/wikipathwayRDF'
         
        self.sourceList.append(self.dataSource3)
         
        self.dataSource4 = DataSource()        
        self.dataSource4.sourceName = 'lipidmaps'
        self.dataSource4.filePrefix = 'lipidmaps'
        self.dataSource4.haveChemClassInfo = True
        self.dataSource4.sourceLocPath = '../misc/output/lipidmaps'
 
        self.sourceList.append(self.dataSource4)
        # End DataSource code
        
        self.dataSource5 = DataSource()        
        self.dataSource5.sourceName = 'rhea'
        self.dataSource5.filePrefix = 'rhea'
        self.dataSource5.haveChemClassInfo = False
        self.dataSource5.sourceLocPath = '../misc/output/rhea_reactions'
   
        self.sourceList.append(self.dataSource5)

        self.dataSource6 = DataSource()
        self.dataSource6.sourceName = 'pfocr'
        self.dataSource6.filePrefix = 'pfocr'
        self.dataSource6.haveChemClassInfo = False
        self.dataSource6.sourceLocPath = '../misc/output/pfocr'

        self.sourceList.append(self.dataSource6)

        self.dataSource7 = DataSource()
        self.dataSource7.sourceName = 'refmet'
        self.dataSource7.filePrefix = 'refmet'
        self.dataSource7.haveChemClassInfo = False
        self.dataSource7.sourceLocPath = '../misc/output/refmet'

        self.sourceList.append(self.dataSource7)
           
        # dictionary that holds data statistics
        self.geneToPathAssocSourceTallies = dict()
        self.metToPathAssocSourceTallies = dict()
        
        # mapping exclusion list and population of the list
        # The population of the exclusion list should be delegated to a method
        self.mappingExclustionList = MappingExclusionList()
        self.mappingExclustionList.populateExclusionList("../config/curation_mapping_issues_list.txt")
    
        # Collection of Molecule objects holding chemical properties.
        self.chemSourceRecords = dict()
        
        # Maps all source ids to their list of alternate ids.
        # Ths captures all primary associations between source id 
        self.sourceIdToIDDict = dict()
        
        # this counts the number of manually curated errors that are still found in the input.
        self.curationAvoidanceCount = 0

        # a simple list of uniprot accessions that are secondary
        # use this for rhea export... to skip any uniprot in this list
        self.uniprotSecondaryAccessions = set()
        
        # suport for delivering the correct proteins for rhea
        self.rheaToSwissprotDict = dict()
        self.rheaToTremblDict = dict()

    def fullBuild(self):
        """
        This high level method performs the entire process of entity construction
        associations and writing. The stages of construction are:
        1.) Constructing pathways and pathway category
        2.) Constructing genes, gene common names, synonyms.
        3.) Building gene/pathway edges.
        4.) Constructing/harmonizing metabolites.
        5.) Adding metabolites common names and synonyms.
        6.) Building metabolite/pathway edges.
        7.) Load chemical properties - currently (1/2021), HMDB and ChEBI
        8.) Associate chemical properties (molecules) with with ramp metabolites
        9.) Write methods: pathways, analyte-source info, synonyms, analyte to pathways, analyte registry, chemical properties.
        """
        # load pathways
        self.loadPathways()
        self.addPathwayCategory()

        # collect uniprot secondary accessions
        self.getUniprotSecondaryAccessions()

        # load genes
        self.loadGeneList()
        self.addGeneCommonNameAndSynonyms()
        self.buildGeneToPathwayConnections()

        # load metabolite list, over all sources and hamonization during build
        self.loadMetaboList()
        self.addMetaboliteCommonName()
        self.addMetaboliteHMDBStatus()
        self.addMetaboliteSynonyms()
        self.buildMetaboliteToPathwayConnections()

        self.loadOntolgies()

        self.loadMetaboliteToGene()
        self.metaboliteClassConnections()

        # load chemistry based on sources, resolveChemistry will attach chem props to metabolites and rampids
        # 1/2021 - currently hmdb and chebi sources
        self.loadChemstry(["hmdb", "chebi", "lipidmaps"])
        self.resolveChemistry(["hmdb", "chebi", "lipidmaps"])      
        
        self.metaboliteList.collapseMetsOnInchiKeyPrefix()

        # Rhea reactions # we have to handle reactions after all mergers because reactions aren't part of the metabolite class, and mergers don't work
        self.processRheaReactions()

        self.metaboliteList.determineBestNames()
        self.geneList.determineBestNames()

        # loader file writes
        self.writePathways()
        self.writeAnalyteSource()
        self.writeAnalyteSynonyms()
        self.writeAnalyteToPathway()
        self.writeAnalyte()
        self.writeMetGeneAssociation()
        self.writeOntologies()
        self.writeOntologyAssociations()
        self.writeChemProps()
        self.writeMetaboliteClass()
        
        self.writeReactionEntities()
        
        print("Number of problem associations skipped (curationAvoidanceCount): " + str(self.curationAvoidanceCount))
        
        
    def loadMetaboList(self, eqMetric = 0):
        """
        Loads the metabolite list for all data sources
        Note that as the list is built, metabolite entities are being merged (subsumed) into RaMP metabolite objects
        """
        Metabolite.__equalityMetric = eqMetric
        
        # Build list for each data source
        for src in self.sourceList:
            
            source = src.sourceName
            
            # Note that the input format is metabolite id dictionaries.
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteIDDictionary.txt"
            
            if not(path.exists(file)):
                continue
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for row in df.itertuples(index=False):
                row = [str(element) for element in row]
                currSourceId, type, altId = row

                if type == "smiles":
                    continue

                
                # add the sourceId and altId to the support dictionary
                if currSourceId not in self.sourceIdToIDDict:
                    self.sourceIdToIDDict[currSourceId] = list()                
                
                excludeMappingConnection = False
                
                # check exclusion mapping list
                excludeMappingConnection = self.mappingExclustionList.isMappingProblem(currSourceId, altId)

                # if it should be excluded, then it's ok to keep this connection
                if not excludeMappingConnection:
                    self.sourceIdToIDDict[currSourceId].append(altId)
                else:
                    self.curationAvoidanceCount = self.curationAvoidanceCount + 1
                
                metabolite = self.metaboliteList.getMetaboliteBySourceId(currSourceId)
                altMetabolite = self.metaboliteList.getMetaboliteBySourceId(altId)
                
                # if the source id isn't already captured, make a metabolite
                if(metabolite is None):
                    
                    if(altMetabolite is not None and not excludeMappingConnection):
                        metabolite = altMetabolite
                    else:
                        metabolite = Metabolite()
                        metabolite.rampId = self.generateRampId("C")
                        metabolite.sourceId = currSourceId
                        if not excludeMappingConnection:   
                            # we have the metabolite and we add it's altId         
                            metabolite.addId(altId, source)                        
                            self.metaboliteList.addMetaboliteByAltId(altId, metabolite)
                        
                    metabolite.addSource(source)
                    metabolite.addId(currSourceId, source)                    
                    self.metaboliteList.addMetaboliteByAltId(currSourceId, metabolite)

                    # this is a sourceId that already exists
                else:
                    # need to check if the alt id already exists as a key id
                    met2 = self.metaboliteList.getMetaboliteBySourceId(altId)
                    
                    if met2 is not None and not excludeMappingConnection:
                        met2.addId(altId, source)
                        met2.addSource(source)
                        met2.addId(currSourceId, source)
                       
                        if not excludeMappingConnection:
                            # this reasigns the primary source id and strands the 'metabolite' record
                            self.metaboliteList.addMetaboliteByAltId(currSourceId, met2)
                            
                            # we have two metabolites, with the same altID
                            # met2 already exists for the altId, 
                            # metabolite exits for the source id
                            # 
                            # if met2 != metabolite - we have two rampIds
                            # we don't want two records
                            # we need to consolidate metabolites... I think
                            if(met2 is not metabolite):
                                # keep the original metabolite (met2) and transfer info
                                met2.subsumeMetabolite(metabolite)
                                
                                # now we need to point references to metabolite to met2
                                #metabolite = met2
                                for id in metabolite.idList:
                                    self.metaboliteList.addMetaboliteByAltId(id, met2)
                                
                        
                    else:
                        # now we have a case where the altId is not already a metabolite
                        # we do have the source id metabolite
                        # we check if the altId should be added or not based on exclusion list 
                        
                        if not excludeMappingConnection:
                            metabolite.addId(altId, source)
                            # lets add the metabolite back in based on the id
                            self.metaboliteList.addMetaboliteByAltId(altId, metabolite)
                            # safe add, adds unique source to metabolite
                            metabolite.addSource(source)
                    

            
    def isaPrimaryIdMapping(self, sourceId, altId):
        """
        Returns True if a source and an alternate id are primary mapping id pairs from the data source.
        """
        return sourceId in self.sourceIdToIDDict and altId in self.sourceIdToIDDict[sourceId]
         
        
    def addMetaboliteCommonName(self):
        """
        Appends metabolite common names to metabolite records for all data sources
        """
        for src in self.sourceList:
            # print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteCommonName.txt"

            if not os.path.exists(file) or os.path.getsize(file) < 1:
                continue

            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
                
            for i,row in df.iterrows():
                met = self.metaboliteList.getMetaboliteBySourceId(row[0])
                if met is not None:
                    met.addCommonName(row[0], row[1], source)
    
        # resolve common name for ids without corresponding common names
        mets = self.metaboliteList.getAllMetabolites()
        for met in mets:
            met.resolveCommonNames()
    
    def addMetaboliteHMDBStatus(self):
        hmdbSrc = None
        for src in self.sourceList:
            if src.filePrefix == 'hmdb':
                hmdbSrc = src
        
        if hmdbSrc is not None:
            # capture id to status dictionary
            hmdbStatus = dict()
            file = hmdbSrc.sourceLocPath + "/" + hmdbSrc.filePrefix + "metStatus.txt"
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')

            for i,row in data.iterrows():
                hmdbStatus[row[0]] = row[1]
            
            # traverse metabolite list
            mets = self.metaboliteList.getAllMetabolites()
            for met in mets:
                idlist = met.idDict.get("hmdb", None)
                if idlist is not None:
                    for id in idlist:
                        status = hmdbStatus.get(id,None)
                        if status is not None:
                            met.setPriorityHMDBStatus(status)

            
           
    def addMetaboliteSynonyms(self):
        """
        Adds all metabolite synonyms for all data sources
        """
        for src in self.sourceList:
            print(src.sourceName)
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metabolitesWithSynonymsDictionary.txt"
            
            if not os.path.exists(file) or os.path.getsize(file) < 1:
                return
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
                
            for i,row in df.iterrows():
                met = self.metaboliteList.getMetaboliteBySourceId(row[0])
                if met is not None:
                    met.addSynonym(row[1], source)                
        
        
                    
    def generateRampId(self, type):
        """
        Generates a unique reamp id for the requested type.
        type is a character C = compound (metabolite), G = gene, P = pathway
        """
        if(type == "C"):
            self.__rampMetStartId = self.__rampMetStartId + 1
            return "RAMP_C_" + (str(self.__rampMetStartId)).zfill(9)
        elif(type == "G"):
            self.__rampGeneStartId = self.__rampGeneStartId + 1
            return "RAMP_G_" + (str(self.__rampGeneStartId)).zfill(9)
        elif(type == "P"):
            self.__rampPathStartId = self.__rampPathStartId + 1
            return "RAMP_P_" + (str(self.__rampPathStartId)).zfill(9)
        elif(type == "OL"):
            self.__rampOntStartId = self.__rampOntStartId + 1
            return "RAMP_OL_" + (str(self.__rampOntStartId)).zfill(9)
        elif(type == "R"):
            self.__rampRxnStartId = self.__rampRxnStartId + 1
            return "RAMP_R_" + (str(self.__rampRxnStartId)).zfill(9)

    def loadPathways(self):
        """
        Loads pathways from all data sources
        """
        for src in self.sourceList:

            #Load Pathway Dictionary First            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "pathwayDictionary.txt"
            if not(path.exists(file)):
                continue
        
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
                        
            for i,row in df.iterrows():
                pathway = Pathway()
                pathway.pathRampId = self.generateRampId("P")
                pathway.pathSource = source                
                pathway.pathSourceId = row[0]
                if(pathway.pathSourceId.startswith("map")):
                    pathway.pathSource = 'kegg'
                pathway.pathName = row[1]
                self.pathList.addPathway(row[0], pathway)
          
    
    def addPathwayCategory(self):
        """
        Loads pathway category for each data source
        """
        for src in self.sourceList:
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "pathwayCategory.txt"
            
            if not(path.exists(file)):
                continue
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, engine='python')

            for i,row in data.iterrows():
                pathway = self.pathList.getPathwayBySourceId(row[0])
                if pathway is not None:
                    if row[1] != np.nan:
                        pathway.pathCategory = row[1]
                    else:
                        pathway.pathCategory = "NA"
        
    
    
    def buildMetaboliteToPathwayConnections(self):
        """
        Loads metabolite/pathway connections.
        Pathway objects from the pathway list are associated with metbolites for each data source.
        """
        strandedMetSourceIds = list()
        strandedPathSourceIds = list()

        for src in self.sourceList:

            #Load Pathway Dictionary First
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metabolitesWithPathwaysDictionary.txt"
            
            if not(path.exists(file)):
                continue
        
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
            
            # print("Number of pathway associations = " + str(data.shape[0]))
        
            map = dict()
            
            self.metToPathAssocSourceTallies[source] = 0
            
            for i,row in df.iterrows():
                
                if row[0] in map:                   
                    map[row[0]].append(row[1])
                else:
                    newlist = list()
                    newlist.append(row[1])
                    map[row[0]] = newlist

            
            i = 0
            assCount = 0
            for metId in map.keys():
                met = self.metaboliteList.getMetaboliteBySourceId(metId)
                if met is not None and len(map[metId]) < 25000:
                    for pathId in map[metId]:
                        pathway = self.pathList.getPathwayBySourceId(pathId)
                        if pathway is not None:
                            met.addPathway(pathway, source)
                            assCount = assCount + 1
                            self.metToPathAssocSourceTallies[source] = self.metToPathAssocSourceTallies[source] + 1
#                            if(assCount % 10000 == 0):
#                                print("associations processsed = " + str(assCount), flush = True)
                                
                
                else:
                    print("high pathway metab: " + metId + " pathway count = " + str(len(map[metId])))


                i = i + 1     
#                if(i % 1000 == 0):
#                    print("metabolites processed = " + str(i), flush=True)

        # print("Finished met to path mapping stranded counts (mets and paths)")
        # print(str(len(strandedMetSourceIds)))
        # print(str(len(strandedPathSourceIds)))
        
    def getUniprotSecondaryAccessions(self):
        
        self.uniprotSecondaryAccessions = set()

        file = "../misc/output/uniprot_human/uniprot_acc_mapping.txt"
        
        data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
        
        for idx, row in data.iterrows():
            altId = row[0]
            id = row[1]
            
            if altId != id:
                self.uniprotSecondaryAccessions.add(altId)
        

    def loadGeneList(self, eqMetric = 0):
        """
        Populates the gene list from all data sources using the <source>geneInfoDictionary files.
        This builds gene entities and merges based on common ids.
        """
        f = open("geneList.log", 'w')

        Metabolite.__equalityMetric = eqMetric
        
        for src in self.sourceList:
#            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "geneInfoDictionary.txt"
            
            if not(path.exists(file)):
                print("in add gene list... geneInfoDictionary not found for :" + file)
                continue
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
            df.drop_duplicates(inplace = True)
                
            for i,row in df.iterrows():
                
                # common names (gene symbols) and secondary ids are ok, but proper names are synonyms
                # rhea has gene names which are synonyms, not proper ids
                if row[1] == 'protein_name' or row[1] == 'gene_name':
                    continue
                
                currSourceId = row[0]
                altId = row[2]
                gene = self.geneList.getGeneById(currSourceId)
                altGene = self.geneList.getGeneById(altId)

                if(gene is None):
                    if altGene is not None:
                        gene = altGene
                    else:    
                        gene = Gene()
                        gene.rampId = self.generateRampId("G")
                        gene.sourceId = currSourceId
                        gene.addId(altId, source)
                        self.geneList.addGene(altId, gene)

                        if gene.rampId == 'RAMP_G_000008086':
                            f.write("Creating our gene...RAMP_G_000008086\n")
                            f.write(gene.rampId+"\n")
                            f.write("\t".join(gene.idList)+"\n")
                            #f.write(gene.idDict)
                
                    gene.addId(currSourceId, source)                        
                    gene.addSource(source)                    
                    self.geneList.addGene(currSourceId, gene)

                    if gene.rampId == 'RAMP_G_000008086':
                       f.write("Adding IDs to our gene...RAMP_G_000008086\n")
                       f.write(gene.rampId + "\n")
                       f.write(currSourceId + "\n")
                       f.write(source + "\n")
                       f.write("\t".join(gene.idList) + "\n")


                    # this is a sourceId lets add 
                else:
                    # need to check if the alt id already exists as a key id
                    gene2 = self.geneList.getGeneById(altId)
                    if gene2 is not None:
                        gene2.addId(altId, source)
                        gene2.addSource(source)
                        gene2.addId(currSourceId, source)


                        if gene2.rampId == 'RAMP_G_000008086':
                            f.write("Linked and adding to our  gene...RAMP_G_000008086\n")
                            f.write(altId + "\n")
                            f.write(source + "\n")
                            f.write(currSourceId + "\n")

                        #metaboliteList.addMataboliteByAltId(altId, met2)
                        # this reasigns the primary source id and strands the 'metabolite' record
                        self.geneList.addGene(currSourceId, gene2)
                        
                        # we have two metabolites, with the same altID
                        # met2 already exists for the altId, 
                        # metabolite exits for the source id
                        # 
                        # if met2 != metabolite - we have two rampIds
                        # we don't want two records
                        # we need to consolidate metabolites... I think
                        if(gene2 is not gene):

                            if("gene_symbol:MDM2" in gene.idList or "gene_symbol:MDM2" in gene2.idList):
                                print("SUBSUME GENE\n")
                               # print(gene.rampId)
                               # print(gene.idList)
                               # print(gene.idDict)

                                print("///\n")
                                #print(gene2.rampId)
                                #print(gene2.idList)
                                #print(gene2.idDict)
                                #print(" ")
                                #print(" ")

                            # keep the original metabolite (met2) and transfer info
                            gene2.subsumeGene(gene)
                            
                            # now we need to point references to metabolite to met2
                            #metabolite = met2
                            for id in gene.idList:
                                self.geneList.addGene(id, gene2)
                            
                        
                    else:
                        gene.addId(altId, source)
                        # lets add the metabolite back in based on the id
                        self.geneList.addGene(altId, gene)
                        # safe add, adds unique source to metabolite
                        gene.addSource(source)
        f.close()



    def loadOntolgies(self):
        """
        loads 5 ontologies per data source, if available.
        This builds the overall onology resource and collects associations
        """
  
        # biofluids
        parentTerm = 'Biofluid and excreta'
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "biofluidLocation.txt"
            
            if not(path.exists(file)):
                continue

            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
            
            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)
                
        # origins
        parentTerm = 'Source'
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "exoEndoDictionary.txt"
            
            if not(path.exists(file)):
                continue
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)

        # cellular location
        parentTerm = 'Subcellular'        
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "cellularLocation.txt"
            
            if not(path.exists(file)):
                continue

            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)

        # tissue location
        parentTerm = 'Tissue and substructures'
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "tissueLocation.txt"
            
            if not(path.exists(file)):
                continue

            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)

        parentTerm = 'Organ and components'
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "organLocation.txt"
            
            if not(path.exists(file)):
                continue

            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)
        
        # metabolite application
        parentTerm = 'Industrial application'      
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteApplication.txt"
            
            if not(path.exists(file)):
                continue
        
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)
                
        # health effect
        parentTerm = 'Health condition'        
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "healthCondition.txt"
            
            if not(path.exists(file)):
                continue
        
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)       
    
        print("ontology size="+str(len(self.ontologyList.getFullOntologyList())))
    
    def recordOntology(self, parentTerm, childTerm, metId):
        if self.resConfig.termIsOnOntologyDenyList(parentTerm, childTerm):
            return
        ontology = self.ontologyList.getOntologyFromParentChild(parentTerm, childTerm)
        if ontology is None:
            ontology = Ontology()
            ontology.ontolParent = parentTerm
            ontology.ontolChild = childTerm
            ontology.ontolRampId = self.generateRampId("OL")
            self.ontologyList.forceAddOntologyRecord(ontology)
        
        # now add ontology record to metabolite
        met = self.metaboliteList.getMetaboliteBySourceId(metId)
        if met is not None:
            met.addOntologyTerm(ontology) 
    
        
    def addGeneCommonNameAndSynonyms(self):
        """
        Adds common name and gene synonyms for all data sources based on geneInfoDictionary files.
        """
        for src in self.sourceList:
#            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "geneInfoDictionary.txt"
            
            if not(path.exists(file)):
                continue
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
                
            for i,row in df.iterrows():
                if row[1] == "common_name" or row[1] == 'protein_name' or row[1] == 'gene_name':
                    gene = self.geneList.getGeneById(row[0])
                    if gene is not None:
                        gene.addCommonNameAndSynonym(row[0], row[2], source, row[1])

        # resolve common name for ids without corresponding common names
        genes = self.geneList.getUniqueGenes()
        
        for gene in genes:
            gene.resolveCommonNames()

    def metaboliteClassConnections(self):
        
        for src in self.sourceList:
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteClass.txt"
                        
            if(path.exists(file) and src.haveChemClassInfo):    
                data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
                df = pd.DataFrame(data)
                df = self.remove_whitespace(df)
                
                for i, row in df.iterrows():               
                    metid = row[0]
                    classLevel = row[1]
                    className = row[2]
                     
                    met = self.metaboliteList.getMetaboliteBySourceId(metid)
                    if met is not None:
                        met.addMetClass(source, metid, classLevel, className)
                  
                                    

    def buildGeneToPathwayConnections(self):
        """
        Constructs gene/pathway connections for all data sources.
        """
        strandedGeneSourceIds = list()
        strandedGeneSourceIds = list()
        noPathwayGenes = 0
        
        for src in self.sourceList:

            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "pathwaysWithGenesDictionary.txt"
                        
            if not(path.exists(file)):
                continue
    
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
            
            self.geneToPathAssocSourceTallies[source] = 0
        
            map = dict()    
            for i,row in df.iterrows():
                
                if row[1] in map:                   
                    map[row[1]].append(row[0])
                else:
                    newlist = list()
                    newlist.append(row[0])
                    map[row[1]] = newlist
            
            i = 0
            assocCount = 0
            for geneId in map.keys():
                gene = self.geneList.getGeneById(geneId)
                if gene is not None and len(map[geneId]) < 25000:
                    for pathId in map[geneId]:
                        pathway = self.pathList.getPathwayBySourceId(pathId)
                        if pathway is not None:
                            gene.addPathway(pathway, source)
                            self.geneToPathAssocSourceTallies[source] = self.geneToPathAssocSourceTallies[source] + 1
                            assocCount = assocCount + 1
 #                           if(assocCount % 10000 == 0):
 #                               print("associations processsed = " + str(assocCount), flush = True)
                                
                
                else:
                    #print("high pathway metab: " + geneId + " pathway count = " + str(len(map[geneId])))
                    noPathwayGenes = noPathwayGenes + 1
#                else:
#                    print("we have a met without a MET")
                i = i + 1     
#                if(i % 1000 == 0):
#                    print("genes processed = " + str(i), flush=True)

    
    def loadMetaboliteToGene(self):
        
        for src in self.sourceList:
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metabolitesLinkedToGenes.txt"
    
            metaboliteToGene = dict()
            geneToProteinType = dict()
            
            if path.exists(file):
#                print ("metaboliteToGene mappings for " + source)
                
                data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False, engine='python')
                df = pd.DataFrame(data)
                df = self.remove_whitespace(df)
                
                for i,row in df.iterrows():
                    met = row[0]
                    gene = row[1]
                    proteinType = row[2]
                    
                    if met not in metaboliteToGene:
                        metaboliteToGene[met] = list()
                    
                    geneToProteinType[gene] = proteinType
                    
                    metaboliteToGene[met].append(gene)

            # now traverse mets 
            for met in metaboliteToGene:
                for gene in metaboliteToGene[met]:
                    
                    targetMet = self.metaboliteList.getMetaboliteBySourceId(met)
                    targetGene = self.geneList.getGeneById(gene)
                    
                    # add protein type to targetGene
                    if targetGene is not None:
                        targetGene.proteinType = geneToProteinType[gene] 
                    
                    if targetMet and targetGene:
                        targetMet.addAssociatedGene(targetGene)
                    
                    
                    
    def processRheaReactions(self):
        
        rheaConfig = None
        
        for source in self.sourceList:
            if source.sourceName == 'rhea':
                rheaConfig = source
        
        if rheaConfig is None:
            return
                
        # rhea output 
        rheaPath = rheaConfig.sourceLocPath
        
        self.buildRxnsFromRhea(rheaPath + "/rhea_primary_records.txt")
        self.appendRxnProteinsFromRhea(rheaPath + "/rhea_uniprot_mapping.txt")
        self.appendRxnParticipantsFromRhea(rheaPath + "/rhea_rxn_to_chebi_and_dir.txt")
        self.dumpReactionToEcEnzymeClass(rheaPath + "/rhea_reaction_to_ec.txt")
                
    def buildRxnsFromRhea(self, path):
        print("Building Rhea Reactions")
        
        records = pd.read_table(path, header = None)
        
        for idx, record in records.iterrows():
            rxn = RheaReaction()
            rxn.rxnRampId = self.generateRampId("R")
            
            rxn.assignPrimaryFields(record)

            self.reactionDict[rxn.rhea_id] = rxn            
        
        
        
    def appendRxnProteinsFromRhea(self, path):
        print("Adding Reaction Proteins")

        # first pull uniprot human accessions
        # this helps to append rhea specific reaction tables and marking human proteins

        records = pd.read_table(path, header=None)

        self.rheaToSwissprotDict = dict()
        self.rheaToTremblDict = dict()

        # just read them in first...
        for idx, record in records.iterrows():
            
            rheaId = record[0]
            uniprot = record[1]
            isReviewed = record[3]

            rxn = self.reactionDict.get(rheaId, None)
            
            if isReviewed == 1:
                idSet = self.rheaToSwissprotDict.get(rheaId, None)
                if idSet is None:
                    idSet = set()
                    self.rheaToSwissprotDict[rheaId] = idSet
                
                idSet.add(uniprot)
            else:
                idSet = self.rheaToTremblDict.get(rheaId, None)
                if idSet is None:
                    idSet = set()
                    self.rheaToTremblDict[rheaId] = idSet
                
                idSet.add(uniprot)
                        
                
            
            if uniprot == "uniprot:Q13574":
                print(" Found the uniprot Q13574 in the FILE.............................................................................")
            if rxn is not None:
                protein = self.geneList.getGeneById(uniprot)
                if protein is not None:
                    if uniprot == "uniprot:Q13574":
                        print(" Found the uniprot Q13574 actual PROTEIN.............................................................................")
                        print("Primary ID in protein = "+protein.sourceId)
                        print("Print all ids in order:")
                        rids = protein.idDict.get("rhea", None)
                        
                        if rids is not None:
                            for id in rids:
                                print(id)
                        
                        print("End Protein Ids List")
                    
                    protein.sourceId = uniprot
                    rxn.proteins.append(protein)
                    protein.isReviewed = isReviewed
                
                else:
                    # we need to make a protein? 
                    # how is the geneList made. Does it include all genes/proteins from Rhea????
                    print("No rxn for uniprot = "+uniprot)
            
            else:
                print("No rxn record = "+rheaId)
                   
        
        
    def appendRxnParticipantsFromRhea(self,path):
        print("Adding Reaction Participants")    
    
        records = pd.read_table(path, header=None)
        
        rheaCofactCnt = 0
        
        for idx, record in records.iterrows():
            rheaId = record[0]
            chebi = record[1]
            rxnSide = record[2]
            chebiCofactor = record[3]
            
            rxn = self.reactionDict.get(rheaId, None)            
            met = self.metaboliteList.getMetaboliteBySourceId(chebi)
            
            if met is None:
                print("hey lack a metabolite for this chebi...: "+chebi)
            else:
                met.isCofactor = chebiCofactor
            
            if rxn is not None and met is not None:
                
                if met.isCofactor == 1:
                    # print("in append rxn members... cofactor = 1 :)")
                    rheaCofactCnt = rheaCofactCnt + 1
                    
                if(rxnSide == 0):
                    rxn.left_comps.append(met)
                    rxn.left_comp_ids.append(chebi)
                else:
                    rxn.right_comps.append(met)    
                    rxn.right_comp_ids.append(chebi)
            else:
                print("in append participants from Rhea... have a None rxn for id: "+rheaId)
        
        print("n/n/n/n/****************************")
        print("Rhea cofact count/est: "+str(rheaCofactCnt))

        
    def dumpReactionToEcEnzymeClass(self, path):

        rxn2EcClassFile = open("../misc/sql/rheaReactionToEcClass.txt", 'w')
        
        with open(path, 'r') as data:
            for line in data:
                #print("reading rh2ec file")
                #print(line)
                sline = line.split("\t")
                rheaId = sline[0]
                #print(rheaId)
                rxn = self.reactionDict.get(rheaId, None)
            
                if rxn is not None:
                    rampRxnId = rxn.rxnRampId
                    if rampRxnId != "":
                        rxn2EcClassFile.write(rampRxnId + "\t" + line)
            
        rxn2EcClassFile.close()
        #print("reaction dict key examples")
        #print(list(self.reactionDict.keys())[0:4])


    def writeAnalyteSource(self):
        """
        Write final files for analyte source
        """
        sourcefile  = open("../misc/sql/analytesource.txt", "w+", encoding='utf-8') 
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
        print("Starting Write of metabolites: size = " + str(len(mets)))
        
        for met in mets:
            s = met.toSourceString()
            
            try:
                sourcefile.write(s + '\n')
            except:
                print("Error writing this record:"+met.rampId)
                print("Error writing this record:"+met.idList[0])
                #print(s)
                pass
         
         
        genes = self.geneList.getUniqueGenes()
        
        print("Starting Write of genes: size = " + str(len(genes)))
        
        for gene in genes:
            s = gene.toSourceString()

            try:
                sourcefile.write(s + '\n')
            except Exception as err:
                print("Error writing this record:"+gene.rampId)
                print("Error writing this record:"+gene.idList[0])
                print(err)
                gene.printGene()
                pass 
           
        sourcefile.close()
        
        
    
    def writeAnalyteToPathway(self):
        """
        write analyte to pathway mappings to final files
        """
        sourcefile  = open("../misc/sql/analytetopathway.txt", "w+", encoding='utf-8')
        
        mets = self.metaboliteList.getAllMetabolites()
        
        for met in mets:
            sourcefile.write(met.toPathwayMapString())
            
        genes = self.geneList.getUniqueGenes()
        
        for gene in genes:
            sourcefile.write(gene.toPathwayMapString())
                    
        sourcefile.close()    
            
            
    def writePathways(self):
        """
        Write pathway object records for all data sources
        """
        sourcefile  = open("../misc/sql/pathway.txt", "w+", encoding='utf-8')

        for pathway in self.pathList.getPathwaysAsList():
            sourcefile.write(pathway.toPathwayString())
            
        sourcefile.close()    


    def writeAnalyte(self):
        """
        Writes analyte list for all data sources.
        """
        sourcefile  = open("../misc/sql/analyte.txt", "w+", encoding='utf-8')

        for met in self.metaboliteList.getUniqueMetabolites():
            sourcefile.write(f"{met.get_insert_format()}\n")
            
        for gene in self.geneList.getUniqueGenes():
            sourcefile.write(f"{gene.get_insert_format()}\n")
                        
        sourcefile.close() 

       
    def writeChemProps(self):
        """
        Writes chemical properties file for all data sources.
        """
        chemPropsFile  = open("../misc/sql/chemProps.txt", "w+", encoding='utf-8')
        mets = self.metaboliteList.getAllMetabolites()
        for met in mets:
            if len(met.chemPropsMolecules) > 0:
                chemPropsFile.write(met.toChemPropsString())

        chemPropsFile.close()


    def writeAnalyteSynonyms(self):
        """
        Writes analyte synonym annotations to file for all data sources
        """
        synonymfile  = open("../misc/sql/analytesynonym.txt", "w+", encoding='utf-8')
        
        mets = self.metaboliteList.getAllMetabolites()
        
        for met in mets:
            s = met.toSynonymsString()        
            synonymfile.write(s)
            
        genes = self.geneList.getUniqueGenes()
        
        for gene in genes:
            s = gene.toSynonymsString()
            synonymfile.write(s)
           
        synonymfile.close()
                

    def writeMetGeneAssociation(self):
        """
        Writes analyte synonym annotations to file for all data sources
        """
        file = open("../misc/sql/catalyzes.txt", "w+", encoding='utf-8')
        
        mets = self.metaboliteList.getAllMetabolites()
        
        for met in mets:
            s = met.toMetToGeneAssociationString()
            if len(s) > 0:
                file.write(s)
        
        file.close()


    def writeOntologies(self):
        """
        Writes analyte synonym annotations to file for all data sources
        """
        file = open("../misc/sql/ontology.txt", "w+", encoding='utf-8')        
        
        ontos = self.ontologyList.getFullOntologyList()
        
        for ont in ontos:
            s = ont.getOntologyString()
            if len(s) > 0:
                file.write(s)
        
        file.close()
    
    
    def writeOntologyAssociations(self):
        mets = self.metaboliteList.getAllMetabolites()

        file = open("../misc/sql/analyteToOntology.txt", "w+", encoding='utf-8')
        for met in mets:
            s = met.toMetaboliteOntologyString()
            if len(s) > 0:
                file.write(s)
            
        file.close()   
        
            
    def writeMetaboliteClass(self):
        mets = self.metaboliteList.getAllMetabolites()

        file = open("../misc/sql/metaboliteClass.txt", "w+", encoding='utf-8')
        for met in mets:
            file.write(met.toMetaboliteClassString())
            
        file.close()    


    def writeReactionEntities(self):
        
        file = open("../misc/sql/reaction.txt", "w+", encoding='utf-8')
        
        for rxnId in self.reactionDict:
            rxn = self.reactionDict[rxnId]
            if rxn.status > 0:
                file.write(rxn.getMainRecordString())
            
        file.close()
        
        cofactorCompCount = 0
        #cofactorCompIdCount = 0
        
        file = open("../misc/sql/reaction_to_metabolite.txt", "w+", encoding='utf-8')

        for rxnId in self.reactionDict:
            rxn = self.reactionDict[rxnId]
            
            cofactorCompCount = cofactorCompCount + rxn.doIHaveACofactorCompoundCheck()

            #cofactorCompIdCount = cofactorCompIdCount + rxn.doIHaveACofactorCompoundIdCheck()
            
            file.write(rxn.getMainReactionToMetString('rhea'))
            
        file.close()
        
        # cofactor checks
        print('n/n/n/n/ cofactor check when writing rhea rxn2met')
        print(cofactorCompCount)
        #print(cofactorCompIdCount)

        file = open("../misc/sql/reaction_to_protein.txt", "w+", encoding='utf-8')
        
        for rxnId in self.reactionDict:
            rxn = self.reactionDict[rxnId]
            swissProtIds = self.rheaToSwissprotDict.get(rxnId, None)
            tremblProtIds = self.rheaToTremblDict.get(rxnId, None)
            if swissProtIds is not None or tremblProtIds is not None:
                file.write(rxn.getMainReactionToProteinStringAllIds('rhea', self.uniprotSecondaryAccessions, swissProtIds, tremblProtIds))
            
        file.close()

        file = open("../misc/sql/reaction_protein_to_metabolite.txt", "w+", encoding='utf-8')
        
        for rxnId in self.reactionDict:
            rxn = self.reactionDict[rxnId]
            file.write(rxn.getReactionProteinToMetString('rhea'))
            
        file.close()


    def remove_whitespace(self, dF):     
        for colName in dF.columns:
            if is_string_dtype(dF[colName]):
                dF[colName] = dF[colName].str.strip()
        return dF
    
    
    
    
    
    def loadChemstry(self, sources):
        """
        Loads chemistry for sources
        """
        cw = ChemWrangler(self.resConfig)
        cw.loadRampChemRecords(sources)
        self.chemSourceRecords = cw.getChemSourceRecords()
        
        
    def resolveChemistry(self, sources):
        """
        Associates molecular entities with parent metabolites
        """
        for source in sources:
            chemRecords = self.chemSourceRecords[source]            
            
            for key in chemRecords:
                met = self.metaboliteList.getMetaboliteBySourceId(key)
                if met is not None:
                    mol = chemRecords[key]
                    met.addChemProps(mol)
        
        self.metaboliteList.printChemPropSummaryStats()
     
    
    
    
    
    
    
    
    def crosscheckChemPropsMW(self, mwTolerance = 0.1, pctOrAbs = 'pct'):
        """
        Utility method to evaluate populated metabolite list and chemical properties.
        This applies a monoisotopic mass cutoff with associated tolerances.
        Any metabolite containing molecules that diverge from the tolerance are added to a problem metabolites list.
        The problem metabolites are returned in a list.
        """
        problemMets = list()
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
        for met in mets:
            dev = met.checkMWParity(mwTolerance, pctOrAbs)
            if dev > 0.0:
                problemMets.append(met)
        
        return problemMets        


    def crosscheckChemPropsInchiBase(self):
        """
        Utility method to verify parity in molecular constituency and connectivity.
        Returns a list of metabolites that contain molecules that have different inchikey base strings.
        The list consists of metabolites with discordant molecule entitities.
        """
        problemMets = list()
        
        mets = self.metaboliteList.getAllMetabolites()
        
        for met in mets:
            uniqueInchiBaseCnt = met.checkInchiBaseParity()
            if uniqueInchiBaseCnt > 1:
                problemMets.append(met)
        
        return problemMets    

    def download_pubchem_molecular_information(self):

        file_path = '../misc/data/chemprops/pubchem/cid_to_mw.tsv'
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        pubchemMW = {}
        if os.path.exists(file_path):
            with open(file_path, 'r', newline='') as tsvfile:
                lines = tsvfile.readlines()
                for row in lines:
                    cid, monoisotopic_mass, inchikey, formula = row.strip().split('\t')
                    pubchemMW[cid] = [monoisotopic_mass, inchikey, formula]
        count = len(pubchemMW)
        print(f'preloading {count} pubchem entries from an earlier run')
        self.loadMetaboList()
        metabolites = self.metaboliteList.getUniqueMetabolites()

        with open(file_path, 'a+') as tsvfile:
            writer: csv = csv.writer(tsvfile, delimiter='\t')
            for met in metabolites:
                for id in met.idList:
                    if id.startswith("pubchem"):
                        if id not in pubchemMW:
                            try:
                                # give pubchem a rest: no more than 5 per second
                                # https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access
                                time.sleep(0.25)
                                c = pcp.Compound.from_cid(id.split(":")[1])
                            except:
                                print("Cid lacks a record at pubchem: " + id)
                                continue
                            count += 1
                            if count % 5000 == 0:
                                print(f"Queried {count} compounds from pubchem")

                            mol_info = [c.monoisotopic_mass, None, None]

                            if c.inchikey is not None:
                                mol_info[1] = c.inchikey
                            if c.molecular_formula is not None:
                                mol_info[2] = c.molecular_formula

                            writer.writerow(
                                [
                                    id,
                                    mol_info[0],
                                    mol_info[1],
                                    mol_info[2]
                                ])
                            pubchemMW[id] = mol_info




    def crossCheckMetaboliteHarmony(self, buildMetAndCompoundProps = True, criteria = "MW", tolerance = 0.1, pctOrAbs = 'pct'):
        """
        A high level utility method that checks for metabolite 'harmony' based on molecular weight or inchi-key prefix.
        The result is an output file containing ramp ids and suspect id pairs in a report.
        Current output is in this directory within the RaMP Backend project:
            ../../misc/resourceConfig/metaboliteMappingIssues.txt
        """
        
        self.loadMetaboList()
        self.addMetaboliteCommonName()
        self.addMetaboliteSynonyms()

        sources = ["hmdb", "chebi", "lipidmaps"]
        self.loadChemstry(sources)
        self.resolveChemistry(sources)
        self.metaboliteList.collapseMetsOnInchiKeyPrefix()

        # load chemistry based on sources, resolveChemistry will attach chem props to metabolites and rampids
        sources = ["hmdb", "chebi", "kegg", "pubchem", "lipidmaps"]
        self.loadChemstry(sources)
        self.resolveChemistry(sources)

        problemMWMets = self.crosscheckChemPropsMW(tolerance, pctOrAbs)

        print("Check mw on mets, problem mets..." + str(len(problemMWMets)))

        totalMismatches = list()
            
        for met in problemMWMets:
            mismatchList = self.getMetaboliteIDMismatchMW(met, tolerance, pctOrAbs)
            
            totalMismatches.extend(mismatchList)

        totalMismatches.sort()
        with open("../misc/metaboliteMappingIssues.txt", "w", encoding='utf-8') as outfile:
            outfile.write("\n".join(totalMismatches))
        outfile.close()


    def get_mw_for_id(self, met, id):
        mws = []
        for source in met.chemPropsMolecules:
            if id in met.chemPropsMolecules[source]:
                mol = met.chemPropsMolecules[source][id]
                if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                    mws.append(float(mol.monoisotopicMass))
        return mws

    def get_common_name_for_id(self, met, id):
        common_names = set()
        for source in met.commonNameDict:
            if id in met.commonNameDict[source]:
                common_names.add(met.commonNameDict[source][id])
        return list(common_names)

    def write_all_ids_and_edges(self, met):
        ids = set(met.idList)
        nodes = []
        relationships = []
        for id in ids:
            nodes.append((id, self.get_mw_for_id(met, id), self.get_common_name_for_id(met, id)))

            if id in self.sourceIdToIDDict:
                linked_ids = self.sourceIdToIDDict[id]
                relationships.extend([(id, other_id) for other_id in linked_ids])

        with open(f"problem_mets_with_no_known_source/{met.rampId}_nodes.tsv", 'w') as tsvfile:
            writer: csv = csv.writer(tsvfile, delimiter='\t')
            writer.writerows(nodes)

        with open(f"problem_mets_with_no_known_source/{met.rampId}_relationships.tsv", 'w') as tsvfile:
            writer: csv = csv.writer(tsvfile, delimiter='\t')
            writer.writerows(relationships)




    #  KJK - refactor to do an all-to-all comparison to find the ID pairs that correspond to very different molecules
    def getMetaboliteIDMismatchMW(self, met, tolerance = 0.1, pctOrAbs = 'pct'):
        """
        Supporting method that looks at all molecules contained in the input Metabolite and finds 
        molecules that fail the test (MW based comparison with a tolerance).
        The returned value is a string representing mismapped molecule pairs.
        The mis-mapped molecules from this method can be output to the curation file.
        """
        
        rampId = met.rampId

        mwDict = dict()
        smileDict = dict()
        rgroupFormulaMolecules = dict()
        idToFormula = dict()

        for source in met.chemPropsMolecules:
            for id in met.chemPropsMolecules[source]:
                mol = met.chemPropsMolecules[source][id]
                smileDict[mol.id] = mol.smiles

                if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                    mwDict[mol.id] = float(mol.monoisotopicMass)

                if mol.formula:
                    idToFormula[id] = mol.formula
                    if "R" in mol.formula:
                        mwDict[mol.id] = float(-1)
                        rgroupFormulaMolecules[id] = mol.formula
                else:
                    idToFormula[id] = ""

        misMatchList = list()
        found_something = False
        for index, id in enumerate(mwDict.keys()):
            mw = mwDict[id]
            smiles = smileDict[id]
            for index2, id2 in enumerate(mwDict.keys()):
                if index >= index2:
                    continue
                mw2 = mwDict[id2]
                smiles2 = smileDict[id2]
                if abs(mw-mw2)/min(mw, mw2) > tolerance or rgroupFormulaMolecules.get(id, False) or rgroupFormulaMolecules.get(id2, False):

                    if self.isaPrimaryIdMapping(id, id2) or self.isaPrimaryIdMapping(id2, id):
                        found_something = True
                        misMatchList.append(
                            f"{rampId}\t{id}\t{mw}\t{id2}\t{mw2}\t{met.toCommonNameJoinString()}\t" +
                            f"{idToFormula.get(id,'')}\t{idToFormula.get(id2,'')}\t{smiles}\t{smiles2}")

        if not found_something:
            self.write_all_ids_and_edges(met)

        return misMatchList

    
class DataSource(object):
    """
    Utility class that holds data source information required for building
    RaMP entities and relations
    """    
    def __init__(self):
        
        self.sourceName = "hmdb"
        self.filePrefix = "hmdb"
        self.sourceLocPath = "../misc/output/hmdb"
        self.haveChemClassInfo = True
        self.exportPath = "../misc/sql"
        
     
class MappingExclusionList(object):        
    '''
     This class holds a collection of improperly associated metabolite ids.
     This list is one means to check associations and skip those that are in error from the source
    '''   
    def __init__(self):
    
        # sourceId : extID List
        self.sourceIdToExtIdDict = dict()
        
    def isMappingProblem(self, sourceId, extId):    

        if sourceId in self.sourceIdToExtIdDict and extId in self.sourceIdToExtIdDict[sourceId]:
            # print("excluded pair (1)"+sourceId+" "+extId)
            return True

        if extId in self.sourceIdToExtIdDict and sourceId in self.sourceIdToExtIdDict[extId]:
            # print("excluded pair (2)"+sourceId+" "+extId)
            return True

        return False
    
    def populateExclusionList(self, filePath):

        data = pd.read_csv(filePath, delimiter=r'\t+', header=0, index_col=None, engine='python')
        df = pd.DataFrame(data)
            
        for i,row in df.iterrows():
            sourceId = row[0]
            extId = row[1]
 
            if sourceId not in self.sourceIdToExtIdDict:
                self.sourceIdToExtIdDict[sourceId] = list()
            if extId not in self.sourceIdToExtIdDict:
                self.sourceIdToExtIdDict[extId] = list()
            
            self.sourceIdToExtIdDict[sourceId].append(extId)   
            self.sourceIdToExtIdDict[extId].append(sourceId)                    
        
        print("Exclusion List Size = " + str(len(list(self.sourceIdToExtIdDict.keys()))))
        

# builder = EntityBuilder()
# builder.crossCheckMetaboliteHarmony(buildMetAndCompoundProps = True, criteria = "MW", tolerance = 0.1, pctOrAbs = 'pct')

# builder.fullBuild()
# print("starting to load metabolites")
# builder.loadMetaboList()
# print("finished metabolite loading")
# print("starting common name")
# builder.addMetaboliteCommonName()
# print("finished common name")
# met = builder.metaboliteList.getMetaboliteBySourceId("hmdb:HMDB0000513")
# 
# 
# if met is not None:
#     print(met.toSourceString())
#     print("have the met!!!!!!!")
# else:
#     print("None met")
# 
# print("\n\n\n\n\n\n\n")
# 
# met = builder.metaboliteList.getMetaboliteBySourceId("hmdb:HMDB0000511")
# 
# 
# if met is not None:
#     print(met.toSourceString())
#     print("have the met!!!!!!!")
# else:
#     print("None met")
#     
#     
# 
# met = builder.metaboliteList.getMetaboliteBySourceId("chebi:30813")
# 
# 
# if met is not None:
#     print(met.toSourceString())
#     print("have the met!!!!!!!")
# else:
#     print("None met")
#     
#     
#     
# met = builder.metaboliteList.getMetaboliteBySourceId("pubchem:2969")
# 
# 
# if met is not None:
#     print(met.toSourceString())
#     print("have the met!!!!!!!")
# else:
#     print("None met")
    
# rConf = RampConfig()
# rConf.loadConfig("../config/external_resource_config.txt")
# eBuild = EntityBuilder(rConf)
# eBuild.fullBuild()
    

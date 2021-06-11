'''
Created on Nov 16, 2020

@author: braistedjc
'''
import os
from os import path
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

from pathlib import Path


import pandas as pd
import numpy as np
from pandas.io.common import file_path_to_url
from pandas.api.types import is_string_dtype
from pandas.io.html import _remove_whitespace
from sqlalchemy.sql.expression import false
from collections import defaultdict

import pubchempy as pcp
import time
#from builtins import True

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

    def __init__(self):
        '''
        Constructor
        '''
        
        # The full metabolite list object
        self.metaboliteList = MetaboliteList()
        
        # The gene list object
        self.geneList = GeneList()
        
        # The pathway list object. Note that pathway entities are also linked to metabolites and genes
        self.pathList = PathwayList()
        
        # ontology list
        self.ontologyList = OntologyList()
        
        # List of DataSource objects. These hold data source configuration.
        self.sourceList = []
        
        # Note: The following populates reactome and wikipathway sources and appends to the default hmdb source
        # This data source list will eventually be populated by config file
        self.source = DataSource()        
        self.sourceList.append(self.source)
    
        self.dataSource2 = DataSource()
        self.dataSource2.sourceName = 'reactome'
        self.dataSource2.filePrefix = 'reactome'
        self.dataSource2.haveChemClassInfo = False
        self.dataSource2.sourceLocPath = '../../misc/output/reactome';
        
        self.sourceList.append(self.dataSource2)
        
        self.dataSource3 = DataSource()        
        self.dataSource3.sourceName = 'wiki'
        self.dataSource3.filePrefix = 'wikipathwayRDF'
        self.dataSource3.haveChemClassInfo = False
        self.dataSource3.sourceLocPath = '../../misc/output/wikiPathwayRDF';        
        
        self.sourceList.append(self.dataSource3)
        
        self.dataSource4 = DataSource()        
        self.dataSource4.sourceName = 'lipidmaps'
        self.dataSource4.filePrefix = 'lipidmaps'
        self.dataSource4.haveChemClassInfo = True
        self.dataSource4.sourceLocPath = '../../misc/output/lipidmaps/';        

        self.sourceList.append(self.dataSource4)
        # End DataSource code
                
        # dictionary that holds data statistics
        self.geneToPathAssocSourceTallies = dict()
        self.metToPathAssocSourceTallies = dict()
        
        # mapping exclusion list and population of the list
        # The population of the exclusion list should be delegated to a method
        self.mappingExclustionList = MappingExclusionList()
        self.mappingExclustionList.populateExclusionList("../../misc/resourceConfig/curation_mapping_issues_list.txt")
    
        # Collection of Molecule objects holding chemical properties.
        self.chemSourceRecords = dict()
        
        # Maps all source ids to their list of alternate ids.
        # Ths captures all primary associations between source id 
        self.sourceIdToIDDict = dict()
        
        
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

        # load genes
        self.loadGeneList()
        self.addGeneCommonNameAndSynonyms()
        self.buildGeneToPathwayConnections()

        # load metabolite list, over all sources and hamonization during build
        self.loadMetaboList()
        self.addMetaboliteCommonName()
        self.addMetaboliteSynonyms()
        self.buildMetaboliteToPathwayConnections()
        
        self.loadOntolgies()
        
        self.loadMetaboliteToGene()        
        self.metaboliteClassConnections()
        
        # load chemistry based on sources, resolveChemistry will attach chem props to metabolites and rampids
        # 1/2021 - currently hmdb and chebi sources
        self.loadChemstry(["hmdb", "chebi", "lipidmaps"])
        self.resolveChemistry(["hmdb", "chebi", "lipidmaps"])      
        
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
                break
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
         
            for i,row in df.iterrows():
            
                if row[1] == "smiles":
                    continue

                currSourceId = str(row[0])                                            
                altId = str(row[2])
                
                # add the sourceId and altId to the support dictionary
                if currSourceId not in self.sourceIdToIDDict:
                    self.sourceIdToIDDict[currSourceId] = list()
                
                self.sourceIdToIDDict[currSourceId].append(altId)
                
                excludeMappingConnection = False
                
                # check exclusion mapping list
                excludeMappingConnection = self.mappingExclustionList.isMappingProblem(currSourceId, altId)
                
                metabolite = self.metaboliteList.getMetaboliteBySourceId(currSourceId)
                altMetabolite = self.metaboliteList.getMetaboliteBySourceId(altId)
                
                # if the source id isn't already captured, make a metabolite
                if(metabolite is None):
                    
                    if(altMetabolite is not None):
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
                    if met2 is not None:
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
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteCommonName.txt"
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
                
            for i,row in df.iterrows():
                met = self.metaboliteList.getMetaboliteBySourceId(row[0])
                if met is not None:
                    met.addCommonName(row[0], row[1], source)
    
        # resolve common name for ids without corresponding common names
        mets = self.metaboliteList.getUniqueMetabolites()
        for met in mets:
            met.resolveCommonNames()
    
    
    def addMetaboliteSynonyms(self):
        """
        Adds all metabolite synonyms for all data sources
        """
        for src in self.sourceList:
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metabolitesWithSynonymsDictionary.txt"
            
            if not os.path.exists(file) or os.path.getsize(file) < 1:
                return
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
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


    def loadPathways(self):
        """
        Loads pathways from all data sources
        """
        for src in self.sourceList:

            #Load Pathway Dictionary First            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "pathwayDictionary.txt"
        
            if not(path.exists(file)):
                break
        
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
                        
            for i,row in df.iterrows():
                pathway = Pathway()
                pathway.pathRampId = self.generateRampId("P")
                pathway.pathSource = source
                pathway.pathSourceId = row[0]
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
                break
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)

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
                break
        
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
            
            print("Number of pathway associations = " + str(data.shape[0]))
        
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

        print("Finished met to path mapping stranded counts (mets and paths)")
        print(str(len(strandedMetSourceIds)))
        print(str(len(strandedPathSourceIds)))
        


    def loadGeneList(self, eqMetric = 0):
        """
        Populates the gene list from all data sources using the <source>geneInfoDictionary files.
        This builds gene entities and merges based on common ids.
        """
        Metabolite.__equalityMetric = eqMetric
        
        for src in self.sourceList:
#            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "geneInfoDictionary.txt"
            
            if not(path.exists(file)):
                break
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
                
            for i,row in df.iterrows():
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

                    gene.addId(currSourceId, source)                        
                    gene.addSource(source)                    
                    self.geneList.addGene(currSourceId, gene)
                
                    # this is a sourceId lets add 
                else:
                    # need to check if the alt id already exists as a key id
                    gene2 = self.geneList.getGeneById(altId)
                    if gene2 is not None:
                        gene2.addId(altId, source)
                        gene2.addSource(source)
                        gene2.addId(currSourceId, source)
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


    def loadOntolgies(self):
        """
        loads 5 ontologies per data source, if available.
        This builds the overall onology resource and collects associations
        """
  
        # biofluids
        parentTerm = 'biofluid'
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "biofluidLocation.txt"
            
            if not(path.exists(file)):
                break

            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
            
            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)
                
        # origins
        parentTerm = 'origins'
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "exoEndoDictionary.txt"
            
            if not(path.exists(file)):
                break
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)

        # cellular location
        parentTerm = 'cellular location'        
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "cellularLocation.txt"
            
            if not(path.exists(file)):
                break

            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)

        # tissue location
        parentTerm = 'tissue location'
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "tissueLocation.txt"
            
            if not(path.exists(file)):
                break

            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)
        
        # metabolite application
        parentTerm = 'application'      
        for src in self.sourceList:
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteApplication.txt"
            
            if not(path.exists(file)):
                break
        
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)

            for i,row in df.iterrows():
                metId = row[0]
                childTerm = row[1]
                self.recordOntology(parentTerm, childTerm, metId)        
    
        print("ontology size="+str(len(self.ontologyList.getFullOntologyList())))
    
    def recordOntology(self, parentTerm, childTerm, metId):
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
                break
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
                
            for i,row in df.iterrows():
                if row[1] == "common_name":
                    gene = self.geneList.getGeneById(row[0])
                    if gene is not None:
                        gene.addCommonNameAndSynonym(row[0], row[2], source)

        # resolve common name for ids without corresponding common names
        genes = self.geneList.getUniqueGenes()
        
        for gene in genes:
            gene.resolveCommonNames()

    def metaboliteClassConnections(self):
        
        for src in self.sourceList:
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteClass.txt"
                        
            if(path.exists(file) and src.haveChemClassInfo):    
                data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)            
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
                break
    
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
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
            
            
            if path.exists(file):
#                print ("metaboliteToGene mappings for " + source)
                
                data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
                df = pd.DataFrame(data)
                df = self.remove_whitespace(df)
                
                for i,row in df.iterrows():
                    met = row[0]
                    gene = row[1]
                    if met not in metaboliteToGene:
                        metaboliteToGene[met] = list()
                    
                    metaboliteToGene[met].append(gene)

            # now traverse mets 
            for met in metaboliteToGene:
                for gene in metaboliteToGene[met]:
                    
                    targetMet = self.metaboliteList.getMetaboliteBySourceId(met)
                    targetGene = self.geneList.getGeneById(gene)
                    
                    if targetMet and targetGene:
                        targetMet.addAssociatedGene(targetGene)
                    
                    
                    
            
                
    
    
    
    
#     def fullBuild(self):
#         """
#         This high level method performs the entire process of entity construction
#         associations and writing. The stages of construction are:
#         1.) Constructing pathways and pathway category
#         2.) Constructing genes, gene common names, synonyms.
#         3.) Building gene/pathway edges.
#         4.) Constructing/harmonizing metabolites.
#         5.) Adding metabolites common names and synonyms.
#         6.) Building metabolite/pathway edges.
#         7.) Load chemical properties - currently (1/2021), HMDB and ChEBI
#         8.) Associate chemical properties (molecules) with with ramp metabolites
#         9.) Write methods: pathways, analyte-source info, synonyms, analyte to pathways, analyte registry, chemical properties.
#         """
#         # load pathways
#         self.loadPathways()
#         self.addPathwayCategory()
# 
#         # load genes
#         self.loadGeneList()
#         self.addGeneCommonNameAndSynonyms()
#         self.buildGeneToPathwayConnections()
# 
#         # load metabolite list, over all sources and hamonization during build
#         self.loadMetaboList()
#         self.addMetaboliteCommonName()
#         self.addMetaboliteSynonyms()
#         self.buildMetaboliteToPathwayConnections()
#         
#         # load chemistry based on sources, resolveChemistry will attach chem props to metabolites and rampids
#         self.loadChemstry()
#         self.resolveChemistry(["hmdb", "chebi"])      
#         
#         # loader file writes
#         self.writePathways()
#         self.writeAnalyteSource()
#         self.writeAnalyteSynonyms()
#         self.writeAnalyteToPathway()
#         self.writeAnalyte()
#         self.writeChemProps()
            


    def writeAnalyteSource(self):
        """
        Write final files for analyte source
        """
        sourcefile  = open("../../misc/sql/analytesource.txt", "w+", encoding='utf-8') 
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
        print("Starting Write of metabolites: size = " + str(len(mets)))
        
        for met in mets:
            s = met.toSourceString()
            
            try:
                sourcefile.write(s)
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
                sourcefile.write(s)
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
        sourcefile  = open("../../misc/sql/analytetopathway.txt", "w+", encoding='utf-8')
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
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
        sourcefile  = open("../../misc/sql/pathway.txt", "w+", encoding='utf-8')

        for pathway in self.pathList.getPathwaysAsList():
            sourcefile.write(pathway.toPathwayString())
            
        sourcefile.close()    


    def writeAnalyte(self):
        """
        Writes analyte list for all data sources.
        """
        sourcefile  = open("../../misc/sql/analyte.txt", "w+", encoding='utf-8')

        for met in self.metaboliteList.getUniqueMetabolites():
            sourcefile.write(met.rampId + "\tcompound\n")
            
        for gene in self.geneList.getUniqueGenes():
            sourcefile.write(gene.rampId + "\tgene\n")
                        
        sourcefile.close() 

       
    def writeChemProps(self):
        """
        Writes chemcial properties file for all data sources.
        """
        chemPropsFile  = open("../../misc/sql/chemProps.txt", "w+", encoding='utf-8')
        mets = self.metaboliteList.getUniqueMetabolites()
        for met in mets:
            if len(met.chemPropsMolecules) > 0:
                chemPropsFile.write(met.toChemPropsString())

        chemPropsFile.close()


    def writeAnalyteSynonyms(self):
        """
        Writes analyte synonym annotations to file for all data sources
        """
        synonymfile  = open("../../misc/sql/analytesynonym.txt", "w+", encoding='utf-8')
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
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
        file = open("../../misc/sql/catalyzes.txt", "w+", encoding='utf-8')
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
        for met in mets:
            s = met.toMetToGeneAssociationString()
            if len(s) > 0:
                file.write(s)
        
        file.close()


    def writeOntologies(self):
        """
        Writes analyte synonym annotations to file for all data sources
        """
        file = open("../../misc/sql/ontology.txt", "w+", encoding='utf-8')        
        
        ontos = self.ontologyList.getFullOntologyList()
        
        for ont in ontos:
            s = ont.getOntologyString()
            if len(s) > 0:
                file.write(s)
        
        file.close()
    
    
    def writeOntologyAssociations(self):
        mets = self.metaboliteList.getUniqueMetabolites()

        file = open("../../misc/sql/analyteToOntology.txt", "w+", encoding='utf-8')
        for met in mets:
            s = met.toMetaboliteOntologyString()
            if len(s) > 0:
                file.write(s)
            
        file.close()   

                
    def writeMetaboliteClass(self):
        mets = self.metaboliteList.getUniqueMetabolites()

        file = open("../../misc/sql/metaboliteClass.txt", "w+", encoding='utf-8')
        for met in mets:
            file.write(met.toMetaboliteClassString())
            
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
        cw = ChemWrangler()
        # sources = ["hmdb","chebi"]
        # sources = ["hmdb","chebi","kegg","pubchem"]
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
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
        for met in mets:
            uniqueInchiBaseCnt = met.checkInchiBaseParity()
            if uniqueInchiBaseCnt > 1:
                problemMets.append(met)
        
        return problemMets    
    
    
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
        
        # load chemistry based on sources, resolveChemistry will attach chem props to metabolites and rampids
        self.loadChemstry(["hmdb", "chebi", "kegg","pubchem", "lipidmaps"])
        self.resolveChemistry(["hmdb", "chebi", "kegg", "pubchem", "lipidmaps"])  

        problemMWMets = self.crosscheckChemPropsMW(tolerance, pctOrAbs)

        print("Check mw on mets, problem mets..." + str(len(problemMWMets)))

        problemInchiMets = self.crosscheckChemPropsInchiBase()

        print("Check inChI base on mets, problem mets..." + str(len(problemInchiMets)))

        # Decided only to use MW criteria, not the more stringent inchikey prefix criteria
        # probMets = list(set(problemMWMets + problemInchiMets))

        # Note - only using molecular weight        
        probMets = problemMWMets
        
        # Note - current method is only using MW criteria. InchiKey prefix is too stringent.
        print("Union of problem mets..." + str(len(probMets)))
        
        moleculeCount = 0;
        for met in probMets:
            for source in met.chemPropsMolecules:
                for id in met.chemPropsMolecules[source]:
                    moleculeCount = moleculeCount + 1
                    
        moleculeCount = 0;
        for met in problemMWMets:
            for source in met.chemPropsMolecules:
                for id in met.chemPropsMolecules[source]:
                    moleculeCount = moleculeCount + 1
                    
        print("Total molecule records (only having MW issue) " + str(moleculeCount))
        
        metCnt = 0
        totalMismatches = list()
        
        if criteria == "MW":
            badMets = problemMWMets
        else:
            badMets = problemInchiMets
            
        for met in badMets:
            
            # s = met.toSourceString()
            # print(s+"\n")
            metCnt = metCnt + 1
            
            if(criteria == "MW"):
                mismatchList = self.getMetaboliteIDMismatchMW(met, tolerance, pctOrAbs)                
            else:
                mismatchList = self.getMetaboliteIDMismatchInchiBase(met)
            
            totalMismatches.extend(mismatchList)
                  
        with open("../../misc/resourceConfig/metaboliteMappingIssues.txt", "w", encoding='utf-8') as outfile:
            outfile.write("\n".join(totalMismatches))
        outfile.close()
    
    
    
        
    def getMetaboliteIDMismatchInchiBase(self, met):
        """
        Returns problem metabolites having compounds with differnt inchikey bases.
        """
        chebiIds = dict()
        hmdbIds = dict()
        pubchemIds = dict()
        rampId = met.rampId
        
        for source in met.chemPropsMolecules:
            for id in met.chemPropsMolecules[source]:
                mol = met.chemPropsMolecules[source][id]
                if mol.id.startswith("hmdb"):
                    if mol.inchiKey and len(mol.inchiKey) > 10:
                        hmdbIds[mol.id] = mol.inchiKey.split("-")[0]
                if mol.id.startswith("chebi"):
                    if mol.inchiKey and len(mol.inchiKey) > 10:
                        chebiIds[mol.id] = mol.inchiKey.split("-")[0]

        for id in met.idList:
            if id.startswith("pubchem"):
                
                if id not in pubchemIds:
                    try:
                        c = pcp.Compound.from_cid(id.split(":")[1])
                    except:
                        continue
                        print("Cid lacks a record at pubchem: " + id)
                    
                    # give pubchem a rest
                    time.sleep(0.5)
                    if c is not None and c.inchikey is not None:
                        pubchemIds[id] = c.inchikey.split("-")[0]
                        print("Retrieved inchikey for pubchem id: " + id)
                
    
        # now reconcile the differences, hmdb to pubchem cid, and chebi to pubchem, and hmdb to chebi 
        misMatchList = list()
        for pid in pubchemIds:
            for hid in hmdbIds:
                if pubchemIds[pid] != hmdbIds[hid]:
                    misMatchList.append(rampId + "\t" + hid + "\t" + hmdbIds[hid] + "\t" + pid + "\t" + pubchemIds[pid] + "\t" + met.toCommonNameJoinString())
            for cid in chebiIds:
                if pubchemIds[pid] != chebiIds[cid]:
                    misMatchList.append(rampId + "\t" + cid + "\t" + chebiIds[cid] + "\t" + pid + "\t" + pubchemIds[pid] + "\t" + met.toCommonNameJoinString())        
        
        # check hmdb to chebi
        for hid in hmdbIds:
            for cid in chebiIds:
                if hmdbIds[hid] != chebiIds[cid]:
                   misMatchList.append(rampId + "\t" + hid + "\t" + hmdbIds[hid] + "\t" + cid + "\t" + chebiIds[cid] + "\t" + met.toCommonNameJoinString())
                                       
        return misMatchList
       
       
       
    def getMetaboliteIDMismatchMW(self, met, tolerance = 0.1, pctOrAbs = 'pct'):
        """
        Supporting method that looks at all molecules contained in the input Metabolite and finds 
        molecules that fail the test (MW based comparison with a tolerance).
        The returned value is a string representing mismapped molecule pairs.
        The mis-mapped molecules from this method can be output to the curation file.
        """
        chebiMW = dict()
        hmdbMW = dict()
        pubchemMW = dict()
        keggMW = dict()
        
        chebiSmiles = dict()
        hmdbSmiles = dict()
        pubchemSmiles = dict()
        keggSmiles = dict()
        
        myFriend = False
        
        rampId = met.rampId
        metPathwayCount = met.getPathwayCount()
        
        for source in met.chemPropsMolecules:
            for id in met.chemPropsMolecules[source]:
                mol = met.chemPropsMolecules[source][id]
                if mol.id.startswith("hmdb"):
                    if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                        hmdbMW[mol.id] = float(mol.monoisotopicMass)
                        hmdbSmiles[mol.id] = mol.smiles
                if mol.id.startswith("chebi"):
                    if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                        chebiMW[mol.id] = float(mol.monoisotopicMass)
                        chebiSmiles[mol.id] = mol.smiles
                    # special deal for kegg R group kegg ids without mass...
                    if "R" in mol.formula:
                        chebiMW[mol.id] = float(-1.0)
                        chebiSmiles[mol.id] = mol.smiles   
                if mol.id.startswith("kegg"):
                    if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                        keggMW[mol.id] = float(mol.monoisotopicMass)
                        keggSmiles[mol.id] = mol.smiles
                    # special deal for kegg R group kegg ids without mass...
                    if "R" in mol.formula:
                        keggMW[mol.id] = float(-1.0)
                        keggSmiles[mol.id] = mol.smiles
                if mol.id.startswith("pubchem"):
                    if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                        pubchemMW[mol.id] = float(mol.monoisotopicMass)
                        pubchemSmiles[mol.id] = mol.smiles
                        
        # a subset of molecules contain R groups.
        # These are special generic molecules that can cause a lot of aggregation
        rgroupFormulaMolecules = dict()
        idToFormula = dict()
        for source in met.chemPropsMolecules:
            for id in met.chemPropsMolecules[source]:
                
                mol = met.chemPropsMolecules[source][id]
                if mol.formula and "R" in mol.formula:                  
                    rgroupFormulaMolecules[id] = mol.formula
                        
                if mol.formula:
                    idToFormula[id] = mol.formula
                else:
                    idToFormula[id] = ""
                    
        
#        12/2021 - temporary code used to evaluate pubchem cids that were associated with
#        The code uses pubchempy and their api to pull in compound attributes used to validate pubchem mapping
#        
#         for id in met.idList:
#              #this will accumulate pubchem monoisotopic masses
#             if id.startswith("pubchem"):
#                  
#                            
#                 if id not in pubchemMW:
#                     try:
#                         c = pcp.Compound.from_cid(id.split(":")[1])
#                     except:
#                         continue
#                         print("Cid lacks a record at pubchem: " + id)
#                      
#                     # give pubchem a rest
#                     time.sleep(0.75)
#                     if c is not None and c.inchikey is not None:
#                         pubchemMW[id] = c.monoisotopic_mass
#                         if c.inchikey is not None:
#                             print(id + "\t" + str(c.monoisotopic_mass) + "\t" + c.inchikey + "\t" + c.molecular_formula)

                   
        # now reconcile the differences, hmdb to pubchem and chebi to pubchem, also hmdb and chebi to kegg
        misMatchList = list()
        for pid in pubchemMW:

            pubchemMass = pubchemMW[pid]
            pubchemSmile = pubchemSmiles.get(pid, "")
            
            for hid in hmdbMW:
                hmdbMass = hmdbMW[hid]
                hmdbSmile = hmdbSmiles.get(hid, "")
                if abs(hmdbMass-pubchemMass)/min(hmdbMass, pubchemMass) > tolerance or rgroupFormulaMolecules.get(pid, False):
                    if self.isaPrimaryIdMapping(hid, pid) or self.isaPrimaryIdMapping(pid, hid):
                        misMatchList.append(rampId + "\t" + hid + "\t" + str(hmdbMass) + "\t" + pid + "\t" + str(pubchemMass) + "\t" + met.toCommonNameJoinString() + "\t" + idToFormula.get(hid, "") + "\t" + idToFormula.get(pid,"") + "\t" + hmdbSmile + "\t" + pubchemSmile)

            for cid in chebiMW:
                chebiMass = chebiMW[cid]
                chebiSmile = chebiSmiles.get(cid, "")
                if abs(chebiMass-pubchemMass)/min(chebiMass, pubchemMass) > tolerance or rgroupFormulaMolecules.get(pid, False):
                    if self.isaPrimaryIdMapping(cid, pid) or self.isaPrimaryIdMapping(pid, cid):
                        misMatchList.append(rampId + "\t" + cid + "\t" + str(chebiMass) + "\t" + pid + "\t" + str(pubchemMass) + "\t" + met.toCommonNameJoinString() + "\t" + idToFormula.get(hid, "") + "\t" + idToFormula.get(cid,"") + "\t" + hmdbSmile  + "\t" + chebiSmile)        
        
        # checking hmdb to chebi    
        for hid in hmdbMW:
            hmdbMass = hmdbMW[hid]
            hmdbSmile = hmdbSmiles.get(hid, "")

            # compare chebi to hmdb
            for cid in chebiMW:
                chebiMass = chebiMW[cid]
                chebiSmile = chebiSmiles.get(cid, "")
                if abs(chebiMass-hmdbMass)/min(chebiMass, hmdbMass) > tolerance or rgroupFormulaMolecules.get(cid, False):
                    if self.isaPrimaryIdMapping(hid, cid) or self.isaPrimaryIdMapping(cid, hid):  
                        misMatchList.append(rampId + "\t" + hid + "\t" + str(hmdbMass) + "\t" + cid + "\t" + str(chebiMass) + "\t" + met.toCommonNameJoinString() + "\t" + idToFormula.get(hid, "") + "\t" + idToFormula.get(cid,"") + "\t" + hmdbSmile + "\t" + chebiSmile)
            # compare kegg to hmdb         
            for keggId in keggMW:
                keggMass = keggMW[keggId]
                keggSmile = keggSmiles.get(keggId, "")
                if abs(keggMass-hmdbMass)/min(keggMass, hmdbMass) > tolerance or rgroupFormulaMolecules.get(keggId, False):
                    if self.isaPrimaryIdMapping(hid, keggId) or self.isaPrimaryIdMapping(keggId, hid):  
                        misMatchList.append(rampId + "\t" + hid + "\t" + str(hmdbMass) + "\t" + keggId + "\t" + str(keggMass) + "\t" + met.toCommonNameJoinString() + "\t" + idToFormula.get(hid, "") + "\t" + idToFormula.get(keggId,"") + "\t" + hmdbSmile + "\t" + keggSmile)
        
        # finish with chebi to kegg                     
        for cid in chebiMW:
            chebiMass = chebiMW[cid]
            chebiSmile = chebiSmiles.get(cid, "")
            for keggId in keggMW:
                keggMass = keggMW[keggId]
                keggSmile = keggSmiles.get(keggId, "")
                if abs(keggMass-chebiMass)/min(keggMass, chebiMass) > tolerance or rgroupFormulaMolecules.get(keggId, False):
                    if self.isaPrimaryIdMapping(cid, keggId) or self.isaPrimaryIdMapping(keggId, cid):  
                        misMatchList.append(rampId + "\t" + cid + "\t" + str(chebiMass) + "\t" + keggId + "\t" + str(keggMass) + "\t" + met.toCommonNameJoinString() + "\t" + idToFormula.get(cid, "") + "\t" + idToFormula.get(keggId,"") + "\t" + chebiSmile + "\t" + keggSmile)
     
                          
        return misMatchList
    
    
    
    def utilCheckHMDBMappingValidity(self):
        '''
        This method checks all mappings from HMDB ids to KEGG and ChEBI Ids.
        In this method to report is meant to capture the number of HMDB IDs that map to old/stale KEGG or ChEBI ids.
        The precondition is that we have loaded metabolites from all sources AND we have loaded chemistry from all sources.
        Each HMDB ID to alternate ID is checked to verify that the ID is valid, in particular against chebi and kegg where
        we have complete records for chemistry.
        '''
        keggTotalMappings = 0
        keggBadMappings = 0
        uniqueBadKegg = dict()

        chebiTotalMappings = 0
        chebiBadMappings = 0
        uniqueBadChebi = dict()
        
        allChebiIds = dict()
        allstarChebiIds = dict()
        uniqueInvalidChebi = dict()
        
        uniqueKegg = dict()
        uniqueChebi = dict()
        
        allChebiFile = "C:\\Tools\\git_projects\\ramp\\RaMP-BackEnd\\misc\\data\\chemprops\\compounds_3star.tsv"
        with open(allChebiFile, "r", encoding='utf-8') as chebiFile:
            for line in chebiFile:
                chebiId = line.split("\t")[0]
                chebiId = chebiId.strip()
                chebiId = "chebi:"+chebiId
                allChebiIds[chebiId] = chebiId
                                
        chebiFile.close()
        
        allChebiFile = "C:\\Tools\\git_projects\\ramp\\RaMP-BackEnd\\misc\\data\\chemprops\\compounds_allstar.tsv"
        with open(allChebiFile, "r", encoding='utf-8') as chebiFile:
            for line in chebiFile:
                chebiId = line.split("\t")[0]
                chebiId = chebiId.strip()
                chebiId = "chebi:"+chebiId
                allstarChebiIds[chebiId] = chebiId
                                
        chebiFile.close()
        
        for id in self.sourceIdToIDDict:
            if id.startswith("hmdb:"):
                for altId in self.sourceIdToIDDict[id]:
                    
                    if altId.startswith("kegg"):
                        keggTotalMappings = keggTotalMappings + 1
                        uniqueKegg[altId] = altId
                        if altId not in self.chemSourceRecords["kegg"]:
                            keggBadMappings = keggBadMappings + 1
                            uniqueBadKegg[altId] = altId
                            # print("bad kegg..." + altId)
 
                    if altId.startswith("chebi"):
                        chebiTotalMappings = chebiTotalMappings + 1
                        uniqueChebi[altId] = altId
                        if altId not in self.chemSourceRecords["chebi"]:                            
                            # check that the chebi isn't just a chebi that lacks a structure in sdf
                            if altId not in allChebiIds:
                                chebiBadMappings = chebiBadMappings + 1
                                uniqueBadChebi[altId] = altId
                                
                                
                            if altId not in allstarChebiIds:
                                uniqueInvalidChebi[altId] = altId    
                                # print("bad chebi..." + altId)
                              
                    # Consider adding pubchem -                 


                                
        print("\n\nCheck of HMDB to KEGG and CHEBI, are the alt ids valid?")
        print("Total mappings to KEGG IDs: " + str(keggTotalMappings))
        print("Unique KEGG ids:" + str(len(uniqueKegg)))
        print("Total mapping to invalid KEGG IDs: "+str(keggBadMappings))
        print("Unique invalid KEGG IDs: "+str(len(uniqueBadKegg)))

        print("\nTotal mappings to ChEBI IDs: " + str(chebiTotalMappings))
        print("Unique Chebi IDs: "+ str(len(uniqueChebi)))
        print("Total mapping to non-3-star ChEBI IDs: "+str(chebiBadMappings))
        print("Unique non-3-star ChEBI IDs: "+str(len(uniqueBadChebi)))
        print("Unique invalid ChEBI IDs: "+str(len(uniqueInvalidChebi)))


                        
    
class DataSource(object):
    """
    Utility class that holds data source information required for building
    RaMP entities and relations
    """    
    def __init__(self):
        
        self.sourceName = "hmdb"
        self.filePrefix = "hmdb"
        self.sourceLocPath = "../../misc/output/hmdb"
        self.haveChemClassInfo = True
        self.exportPath = "../../misc/sql"
        
     
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
            return True

        if extId in self.sourceIdToExtIdDict and sourceId in self.sourceIdToExtIdDict[extId]:
            return True

        return False
    
    def populateExclusionList(self, filePath):

        data = pd.read_csv(filePath, delimiter=r'\t+', header=0, index_col=None)
        df = pd.DataFrame(data)
            
        for i,row in df.iterrows():
            sourceId = row[1]
            extId = row[3]
 
            if sourceId not in self.sourceIdToExtIdDict:
                self.sourceIdToExtIdDict[sourceId] = list()
            if extId not in self.sourceIdToExtIdDict:
                self.sourceIdToExtIdDict[extId] = list()
            
            self.sourceIdToExtIdDict[sourceId].append(extId)   
            self.sourceIdToExtIdDict[extId].append(sourceId)                    
        
        print("Exclusion List Size = " + str(len(list(self.sourceIdToExtIdDict.keys()))))
        

builder = EntityBuilder()
#builder.loadOntolgies()
##builder.writeOntologies()
#builder.crossCheckMetaboliteHarmony(True, "MW", 0.1, 'pct')
#builder.utilCheckHMDBMappingValidity()
builder.fullBuild()

# builder.loadMetaboList()
# builder.addMetaboliteCommonName()
# builder.addMetaboliteSynonyms()
# builder.metaboliteClassConnections()
# builder.writeMetaboliteClass()
# 
# mets = builder.metaboliteList.getUniqueMetabolites()
# 
# metClassCnt = 0
# printedDouble = false
# totAltIds = 0
# uniqueIdDict = dict()
# for met in mets:
#     hmdbMC = met.metClasses.get("hmdb", None)
#     if(hmdbMC is not None):
#         metClassCnt = metClassCnt + 1
#         if metClassCnt == 100:
#             print(met.toMetaboliteClassString())
#         if len(hmdbMC) > 1 and printedDouble == False:
#             print("multi hmdb classes")
#             print(met.toChemPropsString())
#             printedDouble = True
#     
#     totAltIds = totAltIds + len(met.idList)
#     for altId in met.idList:
#         uniqueIdDict[altId] = 1
#     
#        
# print("Met Class Count = "+ str(metClassCnt))
# print("tot id count = "+ str(totAltIds))
# print("unique id count = "+ str(len(uniqueIdDict)))
# 
# met = builder.metaboliteList.getMetaboliteBySourceId("hmdb:HMDB0005402")
# if met is not None:
#     print(met.toMetaboliteClassString())
# 




#         
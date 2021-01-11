'''
Created on Nov 16, 2020

@author: braistedjc
'''
import os
from rampEntity.Gene import Gene
from rampEntity.GeneList import GeneList
from rampEntity.Metabolite import Metabolite
from rampEntity.MetaboliteList import MetaboliteList
from rampEntity.Pathway import Pathway
from rampEntity.PathwayList import PathwayList
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
    classdocs
    '''

    __rampMetStartId = 0
    __rampGeneStartId = 0
    __rampPathStartId = 0

    def __init__(self):
        '''
        Constructor
        '''
        self.metaboliteList = MetaboliteList()
        
        self.geneList = GeneList()
        
        self.pathList = PathwayList()
        
        self.sourceList = []
        
        self.source = DataSource()
        
        self.sourceList.append(self.source)
    
        self.dataSource2 = DataSource()
        
        self.dataSource2.sourceName = 'reactome'
        self.dataSource2.filePrefix = 'reactome'
        self.dataSource2.sourceLocPath = '../../misc/output/reactome';
        
        self.sourceList.append(self.dataSource2)
        
        self.dataSource3 = DataSource()
        
        self.dataSource3.sourceName = 'wiki'
        self.dataSource3.filePrefix = 'wikipathwayRDF'
        self.dataSource3.sourceLocPath = '../../misc/output/wikiPathwayRDF';
        
        self.sourceList.append(self.dataSource3)
        
        self.geneToPathAssocSourceTallies = dict()
        self.metToPathAssocSourceTallies = dict()
        
        self.mappingExclustionList = MappingExclusionList()
        self.mappingExclustionList.populateExclusionList("../../misc/resourceConfig/curation_mapping_issues_list.txt")
    
        self.chemSourceRecords = dict()
        
        self.sourceIdToIDDict = dict()
        
    def loadMetaboList(self, eqMetric = 0):

        Metabolite.__equalityMetric = eqMetric
        
        for src in self.sourceList:
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteIDDictionary.txt"
            
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
                
                # if the source id isn't already captured, make a metabolite
                if(metabolite is None):
                    metabolite = Metabolite()
                    metabolite.sourceId = currSourceId
                    metabolite.addSource(source)
                    metabolite.addId(currSourceId, source)
                    
                    metabolite.rampId = self.generateRampId("C")
                    self.metaboliteList.addMetaboliteByAltId(currSourceId, metabolite)
                    
                    # check that the id mapping is not to be excluded
                    if not excludeMappingConnection:            
                        metabolite.addId(altId, source)
                        self.metaboliteList.addMetaboliteByAltId(altId, metabolite)
                
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
                    
        
        print(" primary id cnt in id mapping dict"+str(len(self.sourceIdToIDDict)))
        for id in self.sourceIdToIDDict["hmdb:HMDB0010338"]:
            print("checking mapped ids for hmdb:HMDB0010338: "+id)
        
        # refine the list...
        # metaboliteList.idBasedMetaboliteMerge()
        
#         print("Finished hmdb test")
#         print("list size = "+str(self.metaboliteList.length()))
#         print("unique met list size = " + str(len(self.metaboliteList.getUniqueMetabolites())))
#                 
#         met = self.metaboliteList.getMetaboliteBySourceId("hmdb:HMDB0038650")
#         if met is not None:
#             met.printMet()
#         
#         met = self.metaboliteList.getMetaboliteBySourceId("pubchem:13592840")
#         if met is not None:
#             met.printMet()
#         
#         met = self.metaboliteList.getMetaboliteBySourceId("CAS:103541-16-8")
#         if met is not None:
#             met.printMet()
#         
#         met = self.metaboliteList.getMetaboliteBySourceId("pubchem:278")
#         if met is not None:
#             met.printMet()
#         else:
#             print("Not a distinct source id:" + "pubchem:278")
#             
#         print("record based on wikidata..."+"\n")
#         met = self.metaboliteList.getMetaboliteBySourceId("wikidata:Q424409")
#         if met is not None:
#             met.printMet()
#         else:
#             print("Not a distinct source id:" + "wikidata:Q424409")
            
    def isaPrimaryIdMapping(self, sourceId, altId):
        return sourceId in self.sourceIdToIDDict and altId in self.sourceIdToIDDict[sourceId]
         
        
    def addMetaboliteCommonName(self):
        
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
        if(type == "C"):
            self.__rampMetStartId = self.__rampMetStartId + 1
            return "RAMP_C_" + (str(self.__rampMetStartId)).zfill(9)
        elif(type == "G"):
            self.__rampGeneStartId = self.__rampGeneStartId + 1
            return "RAMP_G_" + (str(self.__rampGeneStartId)).zfill(9)
        elif(type == "P"):
            self.__rampPathStartId = self.__rampPathStartId + 1
            return "RAMP_P_" + (str(self.__rampPathStartId)).zfill(9)



    def loadPathways(self):
        print("loading pathways")

        for src in self.sourceList:

            #Load Pathway Dictionary First
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "pathwayDictionary.txt"
        
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
        
        
        print("Hey we have pathwaaays... how many?.. "+str(self.pathList.length()))        
        
#         p = self.pathList.getPathwayBySourceId("WP554")        
#         if p is not None:
#             p.printPathway()    
    
    def addPathwayCategory(self):

        for src in self.sourceList:
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "pathwayCategory.txt"
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
            #df = pd.DataFrame(data)
            #df = self.remove_whitespace(df)
            
            for i,row in data.iterrows():
                pathway = self.pathList.getPathwayBySourceId(row[0])
                if pathway is not None:
                    if row[1] != np.nan:
                        pathway.pathCategory = row[1]
                    else:
                        pathway.pathCategory = "NA"
        
        
    def buildMetaboliteToPathwayConnections(self):
        
        strandedMetSourceIds = list()
        strandedPathSourceIds = list()

        for src in self.sourceList:

            #Load Pathway Dictionary First
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metabolitesWithPathwaysDictionary.txt"
    
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
                            if(assCount % 10000 == 0):
                                print("associations processsed = " + str(assCount), flush = True)
                                
                
                else:
                    print("high pathway metab: " + metId + " pathway count = " + str(len(map[metId])))
#                else:
#                    print("we have a met without a MET")
                i = i + 1     
                if(i % 1000 == 0):
                    print("metabolites processed = " + str(i), flush=True)


                

#                 met = self.metaboliteList.getMetaboliteBySourceId(row[0].strip())
#                 if(met is None):
#                     strandedMetSourceIds.append(row[0])
#                 else:
#                     pathway = self.pathList.getPathwayBySourceId(row[1].strip())
#                     if(pathway is None):
#                         strandedPathSourceIds.append(row[1])
#                     else:
#                         # Now we have an association between a metabolite and a pathway
#                         met.addPathway(pathway)
        
        
        print("Finished met to path mapping stranded counts (mets and paths)")
        print(str(len(strandedMetSourceIds)))
        print(str(len(strandedPathSourceIds)))
        


    def loadGeneList(self, eqMetric = 0):

        Metabolite.__equalityMetric = eqMetric
        
        for src in self.sourceList:
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "geneInfoDictionary.txt"
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
                
            for i,row in df.iterrows():
                currSourceId = row[0]
                altId = row[2]
                gene = self.geneList.getGeneById(currSourceId)
                if(gene is None):
                    gene = Gene()
                    gene.sourceId = currSourceId
                    gene.addSource(source)
                    gene.addId(altId, source)
                    gene.rampId = self.generateRampId("G")
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



    def addGeneCommonNameAndSynonyms(self):
        
        for src in self.sourceList:
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "geneInfoDictionary.txt"
            
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


#     def addGeneSynonyms(self):
#         
#         for src in self.sourceList:
#             print(src.sourceName);
#             
#             source = src.sourceName
#             file = src.sourceLocPath + "/" + src.filePrefix + "geneInfoDictionary.txt"
#             
#             data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
#             df = pd.DataFrame(data)
#                 
#             for i,row in df.iterrows():
#                 if row[1] == "common_name":
#                     gene = self.geneList.getGeneById(row[0])
#                     if gene is not None:
#                         gene.addSynony
            

    def buildGeneToPathwayConnections(self):
        
        strandedGeneSourceIds = list()
        strandedGeneSourceIds = list()
        noPathwayGenes = 0
        
        for src in self.sourceList:

            #Load Pathway Dictionary First
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "pathwaysWithGenesDictionary.txt"
    
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None, na_filter = False)
            df = pd.DataFrame(data)
            df = self.remove_whitespace(df)
            
            print("Number of pathway associations = " + str(data.shape[0]))
            
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
                            if(assocCount % 10000 == 0):
                                print("associations processsed = " + str(assocCount), flush = True)
                                
                
                else:
                    #print("high pathway metab: " + geneId + " pathway count = " + str(len(map[geneId])))
                    noPathwayGenes = noPathwayGenes + 1
#                else:
#                    print("we have a met without a MET")
                i = i + 1     
                if(i % 1000 == 0):
                    print("genes processed = " + str(i), flush=True)
                
        print("Pathways are Done for Genes. "+str(noPathwayGenes)+ " genes have no associated pathways")


    def fullBuild(self):
        
        # load pathways
        self.loadPathways()
        self.addPathwayCategory()

        # load genes
        self.loadGeneList()
        self.addGeneCommonNameAndSynonyms()
        self.buildGeneToPathwayConnections()
        
        gene = self.geneList.getGeneById("TFAP2E")
        if gene is not None:
            gene.printGene()

        # load metabolite list, over all sources and hamonization during build
        self.loadMetaboList()
        self.addMetaboliteCommonName()
        self.addMetaboliteSynonyms()
        self.buildMetaboliteToPathwayConnections()
        
        # load chemistry based on sources, resolveChemistry will attach chem props to metabolites and rampids
        self.loadChemstry()
        self.resolveChemistry(["hmdb", "chebi"])      
        
        # loader file writes
        self.writePathways()
        self.writeAnalyteSource()
        self.writeAnalyteSynonyms()
        self.writeAnalyteToPathway()
        self.writeAnalyte()
        self.writeChemProps()


#         metSourceSummary = self.metaboliteList.generateMetaboliteSourceStats(self.sourceList)
#         
#         geneSourceSummary = self.geneList.generateGeneSourceStats(self.sourceList)
#         
#         pathwaySourceSummary = self.pathList.gereratePathwaySourceSummaryStats(self.sourceList)
#                 
#         for source in metSourceSummary.keys():
#             print(source + " metabolite count " + str(metSourceSummary[source]))
#         
#         for source in geneSourceSummary.keys():
#             print(source + " gene count " + str(geneSourceSummary[source]))
# 
#         for source in pathwaySourceSummary.keys():
#             print(source + " pathway count " + str(pathwaySourceSummary[source]))
#             
#         for source in self.metToPathAssocSourceTallies.keys():
#             print(source + " metabolite to pathway association count " + str(self.metToPathAssocSourceTallies[source]))
# 
#         for source in self.geneToPathAssocSourceTallies.keys():
#             print(source + " gene to pathway association count " + str(self.geneToPathAssocSourceTallies[source]))
# 
#         gene = self.geneList.getGeneById("GAPDH")
#         if gene is not None:
#             gene.printGene()
#         
#         print("check GAPDH via main uniprot")
#         gene = self.geneList.getGeneById("uniprot:P04406")
#         if gene is not None:
#             gene.printGene()
#             
#         gene = self.geneList.getGeneById("GAPDXYZ")
#         if gene is not None:
#             gene.printGene()
            


    def writeAnalyteSource(self):

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
        
        sourcefile  = open("../../misc/sql/analytetopathway.txt", "w+", encoding='utf-8')
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
        for met in mets:
            sourcefile.write(met.toPathwayMapString())
            
        genes = self.geneList.getUniqueGenes()
        
        for gene in genes:
            sourcefile.write(gene.toPathwayMapString())
                    
        sourcefile.close()    
            
            
    def writePathways(self):
        
        sourcefile  = open("../../misc/sql/pathway.txt", "w+", encoding='utf-8')

        for pathway in self.pathList.getPathwaysAsList():
            sourcefile.write(pathway.toPathwayString())
            
        sourcefile.close()    


    def writeAnalyte(self):
        sourcefile  = open("../../misc/sql/analyte.txt", "w+", encoding='utf-8')

        for met in self.metaboliteList.getUniqueMetabolites():
            sourcefile.write(met.rampId + "\tcompound\n")
            
        for gene in self.geneList.getUniqueGenes():
            sourcefile.write(gene.rampId + "\tgene\n")
                        
        sourcefile.close() 

       
    def writeChemProps(self):
        
        chemPropsFile  = open("../../misc/sql/chemProps.txt", "w+", encoding='utf-8')
        mets = self.metaboliteList.getUniqueMetabolites()
        for met in mets:
            if len(met.chemPropsMolecules) > 0:
                chemPropsFile.write(met.toChemPropsString())

        chemPropsFile.close()

    def writeAnalyteSynonyms(self):
        
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
                

    def remove_whitespace(self, dF):     
        for colName in dF.columns:
            if is_string_dtype(dF[colName]):
                dF[colName] = dF[colName].str.strip()
                print("fixing column...")
        return dF
    
    
    def loadChemstry(self):
        cw = ChemWrangler()
        sources = ["hmdb","chebi","kegg","pubchem"]
        cw.loadRampChemRecords(sources)
        self.chemSourceRecords = cw.getChemSourceRecords()
        
    def resolveChemistry(self, sourceOrder):
        for source in sourceOrder:
            chemRecords = self.chemSourceRecords[source]            
            
            for key in chemRecords:
                met = self.metaboliteList.getMetaboliteBySourceId(key)
                if met is not None:
                    mol = chemRecords[key]
                    met.addChemProps(mol)
                    
                    if met.rampId == "RAMP_C_000044109" and key == "kegg:C13828":
                        print("\n*\n*\n*\nHeyyyyy in resolve chemistry RAMP 44109 and kegg:C13828\n*\n*\n*\n")
        
        self.metaboliteList.generateChemPropSummaryStats()
     
    
    
    def crosscheckChemPropsMW(self, mwTolerance = 0.1, pctOrAbs = 'pct'):
        
        problemMets = list()
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
        for met in mets:
            dev = met.checkMWParity(mwTolerance, pctOrAbs)
            if dev > 0.0:
                problemMets.append(met)
        
        return problemMets        


    def crosscheckChemPropsInchiBase(self):

        problemMets = list()
        
        mets = self.metaboliteList.getUniqueMetabolites()
        
        for met in mets:
            uniqueInchiBaseCnt = met.checkInchiBaseParity()
            if uniqueInchiBaseCnt > 1:
                problemMets.append(met)
        
        return problemMets    
    
    

    def crossCheckMetaboliteHarmony(self, buildMetAndCompoundProps = True, criteria = "MW", tolerance = 0.1, pctOrAbs = 'pct'):
        
        self.loadMetaboList()
        self.addMetaboliteCommonName()
        self.addMetaboliteSynonyms()
        
        # load chemistry based on sources, resolveChemistry will attach chem props to metabolites and rampids
        self.loadChemstry()
        self.resolveChemistry(["hmdb", "chebi", "kegg", "pubchem"])  

        problemMWMets = self.crosscheckChemPropsMW(tolerance, pctOrAbs)

        print("Check mw on mets, problem mets..." + str(len(problemMWMets)))

        problemInchiMets = self.crosscheckChemPropsInchiBase()

        print("Check inChI base on mets, problem mets..." + str(len(problemInchiMets)))

        probMets = list(set(problemMWMets + problemInchiMets))
        
        # probMets = problemMWMets
        
        print("Union of problem mets..." + str(len(probMets)))
        
        moleculeCount = 0;
        for met in probMets:
            for source in met.chemPropsMolecules:
                for id in met.chemPropsMolecules[source]:
                    moleculeCount = moleculeCount + 1
                    
        print("Total molecule records " + str(moleculeCount))
        
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
            
            if metCnt % 25 == 0:
                print("problem metabolites processed = " + str(metCnt))
    
        
    
        with open("../../misc/resourceConfig/metMappingIssuesReport_withKegg_ONLY_primaryMappings_nowPCQ.txt", "w", encoding='utf-8') as outfile:
            outfile.write("\n".join(totalMismatches))
        outfile.close()
    
        
    def getMetaboliteIDMismatchInchiBase(self, met):
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
        chebiMW = dict()
        hmdbMW = dict()
        pubchemMW = dict()
        keggMW = dict()
        myFriend = False
        
        rampId = met.rampId
        
        if "hmdb:HMDB0010338" in met.idList:
          #if met.rampId == "RAMP_C_000056372": # old case RAMP_C_000044109
            print("Evaluating my test case") 
            print(met.toSourceString())
            myFriend = True
        else:
            myFriend = False
        
        for source in met.chemPropsMolecules:
            for id in met.chemPropsMolecules[source]:
                mol = met.chemPropsMolecules[source][id]
                if mol.id.startswith("hmdb"):
                    if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                        hmdbMW[mol.id] = float(mol.monoisotopicMass)
                if mol.id.startswith("chebi"):
                    if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                        chebiMW[mol.id] = float(mol.monoisotopicMass)
                if mol.id.startswith("kegg"):
                    if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                        keggMW[mol.id] = float(mol.monoisotopicMass)
                        if myFriend:
                            print("have my bad case, getting my kegg mass " + mol.id + " " + str(mol.monoisotopicMass))
                if mol.id.startswith("pubchem"):
                    if mol.monoisotopicMass and len(mol.monoisotopicMass) > 0:
                        pubchemMW[mol.id] = float(mol.monoisotopicMass)
                        
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
#                             print(id + "\t" + str(c.monoisotopic_mass) + "\t" + c.inchikey)
                        
                    
                    
        # now reconcile the differences, hmdb to pubchem and chebi to pubchem, also hmdb and chebi to kegg
        misMatchList = list()
        for pid in pubchemMW:

            pubchemMass = pubchemMW[pid]
            
            for hid in hmdbMW:
                hmdbMass = hmdbMW[hid] 
                if myFriend:
                    print("Test case comparing hmdb to pubchem " + hid + " " + pid )               
                if abs(hmdbMass-pubchemMass)/min(hmdbMass, pubchemMass) > tolerance:
                    if myFriend:
                        print("hmdbmass to pubchem mass (passed cutoffs)" + str(hmdbMass) + " " + str(pubchemMass))                
                    if self.isaPrimaryIdMapping(hid, pid) or self.isaPrimaryIdMapping(pid, hid):
                        if myFriend:
                            print("\n****passed cutoffs AND is a primary mapping***\n")
                        misMatchList.append(rampId + "\t" + hid + "\t" + str(hmdbMass) + "\t" + pid + "\t" + str(pubchemMass) + "\t" + met.toCommonNameJoinString())

            for cid in chebiMW:
                chebiMass = chebiMW[cid]
                if abs(chebiMass-pubchemMass)/min(chebiMass, pubchemMass) > tolerance:
                    if self.isaPrimaryIdMapping(cid, pid) or self.isaPrimaryIdMapping(pid, cid):
                        misMatchList.append(rampId + "\t" + cid + "\t" + str(chebiMass) + "\t" + pid + "\t" + str(pubchemMass) + "\t" + met.toCommonNameJoinString())        
        
        # checking hmdb to chebi    
        for hid in hmdbMW:
            hmdbMass = hmdbMW[hid]
            
            # compare chebi to hmdb
            for cid in chebiMW:
                chebiMass = chebiMW[cid]
                if myFriend:
                    print("comparing " + hid + " to " + cid)
                    print(str(hmdbMass) + " " + str(chebiMass))
                        
                if abs(chebiMass-hmdbMass)/min(chebiMass, hmdbMass) > tolerance:
                    if myFriend:
                        print("test case meets tolerance")
                    if self.isaPrimaryIdMapping(hid, cid) or self.isaPrimaryIdMapping(cid, hid):  
                        misMatchList.append(rampId + "\t" + hid + "\t" + str(hmdbMass) + "\t" + cid + "\t" + str(chebiMass) + "\t" + met.toCommonNameJoinString())
                        if myFriend:
                            print("test case is a primary mapping, meets criteria, being added")
            # compare kegg to hmdb         
            for keggId in keggMW:
                keggMass = keggMW[keggId]
                if myFriend:
                    print("comparing " + hid + " to " + keggId)
                    print(str(hmdbMass) + " " + str(keggMass))
                if abs(keggMass-hmdbMass)/min(keggMass, hmdbMass) > tolerance:
                    if myFriend:
                        print("*****Should be adding hmdb to kegg problem case to bad list.")
                    if self.isaPrimaryIdMapping(hid, keggId) or self.isaPrimaryIdMapping(keggId, hid):  
                        misMatchList.append(rampId + "\t" + hid + "\t" + str(hmdbMass) + "\t" + keggId + "\t" + str(keggMass) + "\t" + met.toCommonNameJoinString())
        
        # finish with chebi to kegg                     
        for cid in chebiMW:
            chebiMass = chebiMW[cid]
            
            for keggId in keggMW:
                keggMass = keggMW[keggId]

                if abs(keggMass-chebiMass)/min(keggMass, chebiMass) > tolerance:
                    if self.isaPrimaryIdMapping(cid, keggId) or self.isaPrimaryIdMapping(keggId, cid):  
                        misMatchList.append(rampId + "\t" + cid + "\t" + str(chebiMass) + "\t" + keggId + "\t" + str(keggMass) + "\t" + met.toCommonNameJoinString())
     
                          
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
    
    def __init__(self):
        
        self.sourceName = "hmdb"
        self.filePrefix = "hmdb"
        self.sourceLocPath = "../../misc/output/hmdb"
        self.exportPath = "../../misc/sql"
        
'''
 This class holds a collection of improperly associated metabolite ids.
 This list is one means to check associations and skip those that are in error from the source
'''        
class MappingExclusionList(object):        
    
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
        
            print("exclusions: " + sourceId + "\t" + extId)
        
            if sourceId not in self.sourceIdToExtIdDict:
                self.sourceIdToExtIdDict[sourceId] = list()
            if extId not in self.sourceIdToExtIdDict:
                self.sourceIdToExtIdDict[extId] = list()
            
            self.sourceIdToExtIdDict[sourceId].append(extId)   
            self.sourceIdToExtIdDict[extId].append(sourceId)                    
        
        print("Exclusion List Size = " + str(len(list(self.sourceIdToExtIdDict.keys()))))
        
        
builder = EntityBuilder()
builder.crossCheckMetaboliteHarmony(True, "MW", 0.1, 'pct')
builder.utilCheckHMDBMappingValidity()

#builder.fullBuild()


# met = builder.metaboliteList.getMetaboliteBySourceId("hmdb:HMDB0128442")
# if met is not None:
#     met.printMet()

# builder.writeMetaboliteSource()
# builder.writeMetaboliteToPathway()
# builder.writePathways()
# builder.writeAnalyte()

# met = builder.metaboliteList.getMetaboliteBySourceId("hmdb:HMDB0002306")
#  
# if met is not None:
#     met.printMet()
# else:
#     print("No metabolite HCl record for hmdb:HMDB0002306")
# 
# met = builder.metaboliteList.getMetaboliteBySourceId("hmdb:HMDB0029225")
#  
# if met is not None:
#     met.printMet()
# else:
#     print("No metabolite Coumeric acid record for HMDB0029225")
# 
# 
# met = builder.metaboliteList.getMetaboliteBySourceId("hmdb:HMDB0000122")
#  
# if met is not None:
#     met.printMet()
# else:
#     print("No metabolite record for glucose hmdb:HMDB0000122")
#     
    

# builder.loadGeneList()
# builder.addGeneCommonName()
# print(str(builder.geneList.length()))
# print(str(len(builder.geneList.getUniqueGenes())))
# builder.loadPathways()
# builder.buildGeneToPathwayConnections()
# 
# 
# gene = builder.geneList.getGeneById("GAPDH")
# if gene is not None:
#     gene.printGene()
# 
# builder.loadMetaboList()
# builder.addMetaboliteCommonName()
#       
# builder.addMetaboliteSynonyms()
# builder.buildMetaboliteToPathwayConnections()
# met = builder.metaboliteList.getMetaboliteBySourceId("chebi:13705")
# 
# if met is not None:
#     met.printMet()
# else:
#     print("Hey... no chebi:13705 metabolite...")
#      
#  
# met = builder.metaboliteList.getMetaboliteBySourceId("hmdb:HMDB0000060")
# if met is not None:
#     met.printMet()
# else:
#     print("Hey... no hmdb:HMDB0000060 metabolite...")
#      
# met = builder.metaboliteList.getMetaboliteBySourceId("kegg:C00164")
# if met is not None:
#     met.printMet()
# else:
#     print("Hey... no kegg:C00164 metabolite...")
#           
#       
# met = builder.metaboliteList.getMetaboliteByAltId("chebi:13705")
# if met is not None:
#     print("have metabolite by alt id query")
#     met.printMet()
# else:
#     print("Hey... no chebi:13705 metabolite...")
#         
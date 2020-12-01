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

import pandas as pd
from pathlib import Path

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
    
    def loadMetaboList(self, eqMetric = 0):

        Metabolite.__equalityMetric = eqMetric
        
        for src in self.sourceList:
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteIDDictionary.txt"
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
            df = pd.DataFrame(data)
                
            for i,row in df.iterrows():
                currSourceId = row[0]
                altId = row[2]
                metabolite = self.metaboliteList.getMetaboliteBySourceId(currSourceId)
                if(metabolite is None):
                    metabolite = Metabolite()
                    metabolite.sourceId = currSourceId
                    metabolite.addSource(source)
                    metabolite.addId(altId, source)
                    metabolite.rampId = self.generateRampId("C")
                    self.metaboliteList.addMetabolite(metabolite)
                
                    # this is a sourceId lets add 
                else:
                    # need to check if the alt id already exists as a key id
                    met2 = self.metaboliteList.getMetaboliteBySourceId(altId)
                    if met2 is not None:
                        met2.addId(altId, source)
                        met2.addSource(source)
                        met2.addId(currSourceId, source)
                        #metaboliteList.addMataboliteByAltId(altId, met2)
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
                        metabolite.addId(altId, source)
                        # lets add the metabolite back in based on the id
                        self.metaboliteList.addMetaboliteByAltId(altId, metabolite)
                        # safe add, adds unique source to metabolite
                        metabolite.addSource(source)
                    
        
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
            

        
        
        
    def addMetaboliteCommonName(self):
        
        for src in self.sourceList:
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metaboliteCommonName.txt"
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
            df = pd.DataFrame(data)
                
            for i,row in df.iterrows():
                met = self.metaboliteList.getMetaboliteBySourceId(row[0])
                if met is not None:
                    met.addCommonName(row[0], row[1], source)
    
    
    
    def addMetaboliteSynonyms(self):
        
        for src in self.sourceList:
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metabolitesWithSynonymsDictionary.txt"
            
            if not os.path.exists(file) or os.path.getsize(file) < 1:
                return
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
            df = pd.DataFrame(data)
                
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
        
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
            df = pd.DataFrame(data)
                        
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
    
        
        
    def buildMetaboliteToPathwayConnections(self):
        
        strandedMetSourceIds = list()
        strandedPathSourceIds = list()

        for src in self.sourceList:

            #Load Pathway Dictionary First
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "metabolitesWithPathwaysDictionary.txt"
    
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
            df = pd.DataFrame(data)
            
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
    
    def loadCommonName(self):
        print("loading common name")
        


    def loadGeneList(self, eqMetric = 0):

        Metabolite.__equalityMetric = eqMetric
        
        for src in self.sourceList:
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "geneInfoDictionary.txt"
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
            df = pd.DataFrame(data)
                
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



    def addGeneCommonName(self):
        
        for src in self.sourceList:
            print(src.sourceName);
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "geneInfoDictionary.txt"
            
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
            df = pd.DataFrame(data)
                
            for i,row in df.iterrows():
                if row[1] == "common_name":
                    gene = self.geneList.getGeneById(row[0])
                    if gene is not None:
                        gene.addCommonName(row[0], row[2], source)

        # resolve common name for ids without corresponding common names
        genes = self.geneList.getUniqueGenes()
        for gene in genes:
            gene.resolveCommonNames()




    def buildGeneToPathwayConnections(self):
        
        strandedGeneSourceIds = list()
        strandedGeneSourceIds = list()
        noPathwayGenes = 0
        
        for src in self.sourceList:

            #Load Pathway Dictionary First
            
            source = src.sourceName
            file = src.sourceLocPath + "/" + src.filePrefix + "pathwaysWithGenesDictionary.txt"
    
            data = pd.read_csv(file, delimiter=r'\t+', header=None, index_col=None)
            df = pd.DataFrame(data)
            
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
        
        self.loadPathways()
        
        self.loadGeneList()
        self.addGeneCommonName()
        self.buildGeneToPathwayConnections()

        self.loadMetaboList()
        self.addMetaboliteCommonName()      
        self.addMetaboliteSynonyms()
        self.buildMetaboliteToPathwayConnections()

        metSourceSummary = self.metaboliteList.generateMetaboliteSourceStats(self.sourceList)
        
        geneSourceSummary = self.geneList.generateGeneSourceStats(self.sourceList)
        
        pathwaySourceSummary = self.pathList.gereratePathwaySourceSummaryStats(self.sourceList)
                
        for source in metSourceSummary.keys():
            print(source + " metabolite count " + str(metSourceSummary[source]))
        
        for source in geneSourceSummary.keys():
            print(source + " gene count " + str(geneSourceSummary[source]))

        for source in pathwaySourceSummary.keys():
            print(source + " pathway count " + str(pathwaySourceSummary[source]))
            
        for source in self.metToPathAssocSourceTallies.keys():
            print(source + " metabolite to pathway association count " + str(self.metToPathAssocSourceTallies[source]))

        for source in self.geneToPathAssocSourceTallies.keys():
            print(source + " gene to pathway association count " + str(self.geneToPathAssocSourceTallies[source]))

class DataSource(object):
    
    def __init__(self):
        
        self.sourceName = "hmdb"
        self.filePrefix = "hmdb"
        self.sourceLocPath = "../../misc/output/hmdb"
            
builder = EntityBuilder()
builder.fullBuild()

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
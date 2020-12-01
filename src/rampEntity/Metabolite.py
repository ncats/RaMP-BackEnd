'''
Created on Nov 6, 2020

@author: braistedjc
'''
import pandas as pd

class Metabolite(object):
    
    # Class variable to indicate the equality policy
    # 0 = ID based, 1 = full lychi-based, 2 = lychi H3 based, 3 = full InchiKey based, 4 = InchiKey prefix match
    __metaboliteEqualityMetric = 0
    
    def __init__(self):
        
        self.sourceId = ""
        
        self.rampId = ""
        
        # uniuqe list of ids
        self.idList = list()

        # source: id dictionary
        self.idDict = dict()
        
        # keys are source, values are dictiontaries of source id to commonName
        self.commonNameDict = dict()
        
        self.synonymDict = dict()
                
        self.primarySource = ""

        self.pathways = dict()
                
        self.smiles = ""
        
        self.inchi = ""
        
        self.inchikey = ""
        
        self.sources = list()
        
        self.lychi = ""
               
    def __eq__(self, other):
        if self.rampId and other.rampId:
            return self.rampId == other.rampId
        return len(set(self.idList).intersection(set(other.idList))) > 0
    
    def shareAltIds(self, other):
            return len(set(self.idList).intersection(set(other.idList))) > 0
       
    def __hash__(self):
        if self.rampId:
            return hash(str(self.rampId))
        else:
            return hash("&".join(self.idList))
        
    def addId(self, id, source):
        # keep this a unique id list
        if id not in self.idList:
            self.idList.append(id)
        if source not in self.idDict:
            self.idDict[source] = list()
        if id not in self.idDict[source]:
            self.idDict[source].append(id)
                
    def printMet(self):
        print("METABOLITE RECORD\n")
        s= "RampID: " + self.rampId + "\n"
        s= s + "sourceId " + self.sourceId + "\n"
        
        print(s)
        
        s = ""
        for src in self.idDict:
            for id in self.idDict[src]:
                s = s + "source: " + src + " altId: " + id + "\n"
        print(s)
        
        s=""
        for src in self.commonNameDict:
            for id in self.commonNameDict[src]:
                s = s + "source: " + src + " id: " + id + " commonName: " + self.commonNameDict[src][id] + "\n"
            
        print(s)
        
        s=""
        for src in self.synonymDict:
            for syn in self.synonymDict[src]:
                s = s + "source: " + src + " syn: " + syn + "\n"
        
        print("pathways")
        for source in self.pathways.keys():
            print(source + " pathways")
            for pathway in self.pathways[source]:
                pathway.printPathway()
            print("")
        
    def addSource(self, sourceName):
        # maintain as a unique list
        if sourceName not in self.sources:
            self.sources.append(sourceName)
        
    def getSortedSources(self):
        srcs = set(self.sources)
        list(srcs).sort()
        return srcs
    
    def mergeMets(self, otherMet):
        for id in otherMet.idList:
            self.idList.append(id)
        for src in otherMet.sources:
            self.sources.append(src)
            
    def addPathway(self, pathway, source):
        if source not in self.pathways:
            self.pathways[source] = list()
        if pathway not in self.pathways[source]:
            self.pathways[source].append(pathway)
            
    def addCommonName(self, id, commonName, source):
        if source not in self.commonNameDict:            
            self.commonNameDict[source] = dict()
        self.commonNameDict[source][id] = commonName
        
    def addSynonym(self, synonym, source):
        if source not in self.synonymDict:
            self.synonymDict[source] = list()
        if synonym not in self.synonymDict[source]:
            self.synonymDict[source].append(synonym)
            
    def subsumeMetabolite(self, metabolite):
        # copy ids
        for source in metabolite.idDict:
            for id in metabolite.idDict[source]:
                self.addId(id, source)
        for source in metabolite.commonNameDict:
            for id in metabolite.commonNameDict[source]:
                self.addCommonName(id, metabolite.commonNameDict[source][id] ,source)
        for source in metabolite.sources:
            self.addSource(source)
        for source in metabolite.synonymDict:
            for syn in metabolite.synonymDict[source]:
                self.addSynonym(syn, source)
                
    
    
                
        
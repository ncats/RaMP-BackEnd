'''
Created on Nov 24, 2020

@author: braistedjc
'''

import pandas as pd

class Gene(object):
    '''
    classdocs
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.sourceId = ""
        
        self.rampId = ""
        
        # uniuqe list of ids
        self.idList = list()

        # source: id dictionary
        self.idDict = dict()
        
        # keys are source, values are dictiontaries of source id to commonName
        self.commonNameDict = dict()
        
        self.synonymList = list()
                
        self.primarySource = ""

        self.pathways = list()
                        
        self.sources = list()
        
    
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
            
    def addSource(self, sourceName):
        # maintain as a unique list
        if sourceName not in self.sources:
            self.sources.append(sourceName)
            
    def addCommonName(self, id, commonName, source):
        if source not in self.commonNameDict:            
            self.commonNameDict[source] = dict()
        self.commonNameDict[source][id] = commonName
        self.addSynonym(commonName)
        
    def addSynonym(self, synonym):
        if synonym not in self.synonymList:
            self.synonymList.append(synonym)
            
    def addPathway(self, pathway):
        if pathway not in self.pathways:
            self.pathways.append(pathway)        
        
    def printGene(self):
        s = "rampId: " + self.rampId + "\n"        
        for source in self.idDict:
            for id in self.idDict[source]:
                s = s + "id: " + id + ", source: " + source + "\n"
        for source in self.commonNameDict:
            for id in self.commonNameDict[source].keys():
                s = s + "id: " + id + ", source: " + source + ", name:" + self.commonNameDict[source][id]+ "\n"
        for syn in self.synonymList:
            s = s + "syn: " + syn + "\n"
        
        print(s)
        
        print("pathway count: " + str(len(self.pathways)))
        
        for pathway in self.pathways:
            pathway.printPathway()
            
             
    
    
    '''
    The original code in writeToSQL would steal the commonName for a gene
    from another common name when one instance of the gene (from the same source) lacked a common name.
    It's a bit of a hack but for genes it will complete this field for the source table.
    Also, all common names will have distinct rows with the same rampId.
    This means that common name will resolve to a ramp id from any common name available.
    '''
    def resolveCommonNames(self):
        for source in self.idDict:
            for id in self.idDict[source]:
                if source in self.commonNameDict and id not in self.commonNameDict[source]:
                    # now we know we have a common name dictionary for the source
                    # and our id doesn't have a commmon name entry.
                    
                    #grab a key
                    keyId = list(self.commonNameDict[source].keys())[0]
                    self.commonNameDict[source][id] = self.commonNameDict[source][keyId]
        
    
    
    def subsumeGene(self, gene):
        # copy ids
        for source in gene.idDict:
            for id in gene.idDict[source]:
                self.addId(id, source)
        for source in gene.commonNameDict:
            for id in gene.commonNameDict[source]:
                self.addCommonName(id, gene.commonNameDict[source][id] ,source)
        for source in gene.sources:
            self.addSource(source)
        for syn in gene.synonymList:
            self.addSynonym(syn)        
            
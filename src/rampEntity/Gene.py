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
        
        self.synonymDict = dict()
                
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
        for source in gene.synonymDict:
            for syn in gene.synonymDict[source]:
                self.addSynonym(syn, source)        
            
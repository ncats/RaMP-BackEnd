'''
Created on Nov 24, 2020

@author: braistedjc
'''

class GeneList(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.geneList = dict()
        
        self.sourceSummary = dict()
        
    def addGene(self, id, gene):
        self.geneList[id] = gene
    
    def getGeneById(self, id):
        return self.geneList.get(id)  
    
    def length(self):
        return len(self.geneList)
    
    def getUniqueGenes(self):
        return list(set(self.geneList.values()))
    
    
    def generateGeneSourceStats(self, sourceList):
        
        for source in sourceList:
            self.sourceSummary[source.sourceName] = 0
        
        genes = self.getUniqueGenes()
        
        for gene in genes:
            for source in gene.sources:
                self.sourceSummary[source] = self.sourceSummary[source] + 1
    
        return self.sourceSummary
    
    
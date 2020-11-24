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
        
    def addGene(self, id, gene):
        self.geneList[id] = gene
    
    def getGeneById(self, id):
        return self.geneList.get(id)  
    
    def length(self):
        return len(self.geneList)
    
    def getUniqueGenes(self):
        return list(set(self.geneList.values()))
    
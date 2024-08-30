'''
Created on Nov 24, 2020

@author: braistedjc
'''

class GeneList(object):
    '''
    The GeneList class holds RaMP Gene objects indexed by source id.
    The list objects had helper methods for access and list statistics
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.geneList = dict()
        
        self.sourceSummary = dict()
        
    def addGene(self, id, gene):
        """Adds a gene for a given id"""
        self.geneList[id] = gene
    
    def getGeneById(self, id):
        """Returns the gene based on the supplied id, else None"""
        return self.geneList.get(id)  
    
    def length(self):
        """Returns the length of the list"""
        return len(self.geneList)
    
    def getUniqueGenes(self):
        """Returns a list object of all genes"""
        return list(set(self.geneList.values()))
    
    
    def generateGeneSourceStats(self, sourceList):
        """A utility method to return gene count for each data source"""
        for source in sourceList:
            self.sourceSummary[source.sourceName] = 0
        
        genes = self.getUniqueGenes()
        
        for gene in genes:
            for source in gene.sources:
                self.sourceSummary[source] = self.sourceSummary[source] + 1
    
        return self.sourceSummary
    
    def determineBestNames(self):
        for gene in self.getUniqueGenes():
            most_common_uniprot = gene.get_most_common_value('sourceId', 'uniprot')
            if most_common_uniprot is not None:
                gene.representativeName=most_common_uniprot
            else:
                gene.representativeName=gene.get_most_common_value('sourceId')
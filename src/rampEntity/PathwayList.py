'''
Created on Nov 6, 2020

@author: braistedjc
'''

class PathwayList(object):
    '''
    Pathway list object. Holds a collection of pathways indexed by pathway ids.
    For convenience the collection exists as a single list as well indexed by pathway source.
    '''


    def __init__(self):
        '''
        Constructor
        '''
        # general source id to pathway dictionary
        self.pathwayDict = dict()
        
        # The same pathway objects organized by source
        self.pathwayBySourceDict = dict()
        
        # Holds summary statistics by data source
        self.sourceSummary = dict()
        
    def addPathway(self, id, pathway):
        """
        Adds a pathway by source id.
        Note that the pathway is also indexed by it's pathway::source 
        """
        self.pathwayDict[id] = pathway
        if self.pathwayBySourceDict.get(pathway.pathSource) is None:
            self.pathwayBySourceDict[pathway.pathSource] = list()
        self.pathwayBySourceDict[pathway.pathSource].append(pathway)
        
        
    def getPathwayBySourceId(self, id):
        """
        Utility method to return a pathway by its source id.
        """
        return self.pathwayDict[id]
        
    def length(self):
        """
        Returns the number of distinct pathways
        """
        return (len(self.pathwayDict))
    
    def getPathwaysAsList(self):
        """
        Returns all pathways as a list object
        """
        return self.pathwayDict.values()
    
    def gereratePathwaySourceSummaryStats(self, sourceList):
        """
        Generates and returns a hash of source to pathway count
        """
        for source in self.pathwayBySourceDict.keys():
            self.sourceSummary[source] = len(self.pathwayBySourceDict[source])
            
        return self.sourceSummary    
        
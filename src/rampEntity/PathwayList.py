'''
Created on Nov 6, 2020

@author: braistedjc
'''

class PathwayList(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.pathwayDict = dict()
        
        self.pathwayBySourceDict = dict()
        
        self.pathwaySummaryStats = dict()
        
        self.sourceSummary = dict()
    
    def addPathway(self, id, pathway):
        self.pathwayDict[id] = pathway
        if self.pathwayBySourceDict.get(pathway.pathSource) is None:
            self.pathwayBySourceDict[pathway.pathSource] = list()
        self.pathwayBySourceDict[pathway.pathSource].append(pathway)
        
    def getPathwayBySourceId(self, id):
        return self.pathwayDict[id]
        
    def length(self):
        return (len(self.pathwayDict))
    
    def getPathwaysAsList(self):
        return self.pathwayDict.values()
    
    def gereratePathwaySourceSummaryStats(self, sourceList):
           
        for source in self.pathwayBySourceDict.keys():
            self.sourceSummary[source] = len(self.pathwayBySourceDict[source])
            
        return self.sourceSummary    
        
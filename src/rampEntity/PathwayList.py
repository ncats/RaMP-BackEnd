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
        
    
    def addPathway(self, id, pathway):
        self.pathwayDict[id] = pathway
        
    def getPathwayBySourceId(self, id):
        return self.pathwayDict[id]
        
    def length(self):
        return (len(self.pathwayDict))
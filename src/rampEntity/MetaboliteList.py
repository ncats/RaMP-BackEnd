'''
Created on Nov 6, 2020

@author: braistedjc
'''
from rampEntity.Metabolite import Metabolite
from itertools import permutations 

class MetaboliteList(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.metaboliteSourceIdList = dict()
        
        self.sourceSummary = dict()
                
    def contains(self, metabolite):
        return (self.metaboliteSourceIdList[metabolite.sourceId] != None)

    def getMetaboliteBySourceId(self, sourceId):
        return self.metaboliteSourceIdList.get(sourceId)   
    
    def addMetabolite(self, metabolite):
        self.metaboliteSourceIdList[metabolite.sourceId] = metabolite
        
    def addMetaboliteByAltId(self, id, metabolite):
        self.metaboliteSourceIdList[id] = metabolite
        
    def length(self):
        return len(self.metaboliteSourceIdList)
    
    def getUniqueMetabolites(self):
        return list(set(self.metaboliteSourceIdList.values()))
    
#     def generateMetaboliteSourceStats(self):
#         mets = self.getUniqueMetabolites()
#         sourceDict = dict()
#         
# 
#         for met in mets:
#             sources = ','.join(met.getSortedSources())
#             
#             if sources in sourceDict:
#                 sourceDict[sources] = sourceDict[sources] + 1
#             else:
#                 sourceDict[sources] = 1
#                 
#         combos = list(sourceDict.keys())
#         combos.sort()
#         
#         for combo in combos:
#             print(combo + "  " + str(sourceDict[combo]) + "\n")
        
    
    def generateMetaboliteSourceStats(self, sourceList):
        
        for source in sourceList:
            self.sourceSummary[source.sourceName] = 0
        
        mets = self.getUniqueMetabolites()
        
        for met in mets:
            for source in met.sources:
                self.sourceSummary[source] = self.sourceSummary[source] + 1
    
        return self.sourceSummary
                        
                    
    def generateChemPropSummaryStats(self):

        mets = self.getUniqueMetabolites()
        
        haveMolCount = 0
        molRecords = 0;

        for met in mets:
            if len(met.chemPropsMolecules) > 0:
                haveMolCount = haveMolCount + 1
                
                for source in met.chemPropsMolecules:
                    molRecords = molRecords + len(met.chemPropsMolecules[source])
                    
        print("\nChemistry Property Stats")
        print("Tot Mets: " + str(len(mets)))
        print("Tot Mets with ChemProps: " + str(haveMolCount))
        print("Tot Mets: " + str(molRecords))
        
        
                    
        
    
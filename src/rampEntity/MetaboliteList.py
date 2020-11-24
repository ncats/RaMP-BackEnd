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
    
    def generateMetaboliteSourceStats(self):
        mets = self.getUniqueMetabolites()
        sourceDict = dict()
        

        for met in mets:
            sources = ','.join(met.getSortedSources())
            
            if sources in sourceDict:
                sourceDict[sources] = sourceDict[sources] + 1
            else:
                sourceDict[sources] = 1
                
        combos = list(sourceDict.keys())
        combos.sort()
        
        for combo in combos:
            print(combo + "  " + str(sourceDict[combo]) + "\n")
        
    
#     def idBasedMetaboliteMerge(self):
#         mets = self.getUniqueMetabolites()
#         
#         iters = len(mets)*len(mets)
#         print("met squared..." + str(iters))
#         checkCnt = 0
#         for met in mets:
#             for met2 in mets:
#                 checkCnt = checkCnt + 1
#                 if(checkCnt % 10000000 == 0):
#                     print("iter: "+str(checkCnt)+ " %done: " + str((checkCnt/iters)*100))
#                 if(met2.shareAltIds(met)):
#                     if(met2.rampId != met.rampId):
#                         met.mergeMets(met2)
#                         self.metaboliteSourceIdList[met2.sourceId] = met
                        
                    
        
        
        
    
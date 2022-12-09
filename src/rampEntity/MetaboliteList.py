'''
Created on Nov 6, 2020

@author: braistedjc
'''
from rampEntity.Metabolite import Metabolite
from itertools import permutations 

class MetaboliteList(object):
    '''
    Container list class holding metabolites. The list class supports access and set methods as well as
    utility methods to extract metabolite records. 
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.metaboliteSourceIdList = dict()
        
        self.inchikeyPrefixToMetab = dict()
        
        self.sourceSummary = dict()
        
                
    def contains(self, metabolite):
        """Utility method returning true if the metabolite exists in the list or false if 
        the metabolite does not exist in the list based on source id. 
        """
        return (self.metaboliteSourceIdList[metabolite.sourceId] != None)


    def getMetaboliteBySourceId(self, sourceId):
        """
        returns a metabolite based on the source id
        """
        return self.metaboliteSourceIdList.get(sourceId)   
    
    
    def addMetabolite(self, metabolite):
        """
        Adds a metabolite to the list
        """
        self.metaboliteSourceIdList[metabolite.sourceId] = metabolite
    
        
    def addMetaboliteByAltId(self, id, metabolite):
        """
        Adds a metabolite connected to a supplied id.
        This permits multiple ids to point to a shared metabolite.        
        """
        self.metaboliteSourceIdList[id] = metabolite
        
        
    def length(self):
        """
        Returns the number of source ids, not distinct metabolites.
        """
        return len(self.metaboliteSourceIdList)
    
    def getUniqueMetabolites(self):
        """
        Returns the set of unique metabolite entities as a list.
        """
        return list(set(self.metaboliteSourceIdList.values()))
    
            
    def generateMetaboliteSourceStats(self, sourceList):
        """
        Returns a map of source to count of metabolites from each source.
        """
        for source in sourceList:
            self.sourceSummary[source.sourceName] = 0
        
        mets = self.getUniqueMetabolites()
        
        for met in mets:
            for source in met.sources:
                self.sourceSummary[source] = self.sourceSummary[source] + 1
    
        return self.sourceSummary
                        
                    
    def printChemPropSummaryStats(self):
        """
        Utility to print chemical property summary statistics for chemical properties. 
        """
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
        print("Tot Molecules records: " + str(molRecords))
        
        
    
    def buildInchiKeyPrefixToMetaboliteMapping(self):
        
        self.inchikeyPrefixToMetab
        
        inchiPrefMap = dict()
        noInchiMet = list()
        
        mets = self.getUniqueMetabolites()
        
        # build a mapping from inchi prefix to metabolites
        for met in mets:
            haveInchiPrefix = False
            inchiKeyPrefixes = met.getInchiPrefixes()
            for prefix in inchiKeyPrefixes:
                haveInchiPrefix = True
                metList = self.inchikeyPrefixToMetab.get(prefix, None)
                if metList is None:
                    metList = list()
                    metList.append(met)
                    self.inchikeyPrefixToMetab[prefix] = metList
                else:
                    metList.append(met)
#             if len(met.chemPropsMolecules) > 0:
#                 for source in met.chemPropsMolecules:
#                     
#                     molDict = met.chemPropsMolecules[molname]
#                     for sourceId in molDict:
#                         mol = molDict[sourceId]
#                         if len(mol.inchiKeyPrefix) > 0:
#                             haveInchiPrefix = True
#                             metList = self.inchikeyPrefixToMetab.get(mol.inchiKeyPrefix, None)
#                             if metList is None:
#                                 metList = list()
#                                 metList.append(met)
#                                 self.inchikeyPrefixToMetab[mol.inchiKeyPrefix] = metList
#                             else:
#                                 metList.append(met)
            if not haveInchiPrefix:
                noInchiMet.append(met)
                
    # if a metabolite has chemical properties, and the chem props include the inchi prefix of interest
    # return the list of metabolites                
    def collapseMetaboliteListOnInchiPrefix(self):
                
        # now we collapse, we need to reach out to others
        inchiPrefixesHit = []
        inchiPrefixesToRemove = []
        
        rampIdSets = dict()
        
        for inchiPrefix in self.inchikeyPrefixToMetab:
            currMets = []
            inchiPrefixesHit = []
            metList = self.inchikeyPrefixToMetab[inchiPrefix]
            rampId = ""
            
            # Just keep the same metabolite objects, but set the ramp id to be the same :)
            if len(metList) > 1:
                for met in metList:
                    for met2 in metList:
                        if met != met2:
                            met.addInchiNeighbor(met2)
            
        # now each metabolite linked by an inchi has a list of it's neighbors
        # neighbors propagate through neighbors in each addition so that neighbors meet neighbors
        # neighbors are connected by shared inchi prefixes
        # now a smart traversal should set common ramp ids across all neighbors
        mets = self.getUniqueMetabolites()
        
        touchedMets = list()
        for met in mets:
            inchiNeigbors = met.getInchiNeighborhood()    
            
            rampId = met.rampId
            for neighbor in inchiNeigbors:
                # if we've not examined a neighborhood, then set the common ramp id
                if neighbor not in touchedMets:
                    touchedMets.append(neighbor)
                    neighbor.rampId = rampId                    
                    
            
    
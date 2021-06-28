'''
Created on Jun 3, 2021

@author: braistedjc
'''
from rampEntity.Ontology import Ontology

class OntologyList(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.simpleOntolList = list()
        
        self.ontolDict = dict()
        
        self.ontolDict["origins"] = dict()
        self.ontolDict["biofluid"] = dict()
        self.ontolDict["cellular location"] = dict()
        self.ontolDict["tissue location"] = dict()
        self.ontolDict["application"] = dict()
        self.ontolDict["health effect"] = dict()
                                     
    def addOntologyRecord(self, ontology):
        if ontology not in self.simpleOntolList:
            print("append ontology")
            self.simpleOntolList.append(ontology)
            self.ontolDict[ontology.ontolParent][ontology.ontolChild] = ontology
        
    def forceAddOntologyRecord(self, ontology):
        self.simpleOntolList.append(ontology)
        self.ontolDict[ontology.ontolParent][ontology.ontolChild] = ontology    
                
    def getOntologyFromParentChild(self, parentTerm, childTerm):
        ontology = None
        if self.ontolDict[parentTerm] is not None:
            ontology = self.ontolDict[parentTerm].get(childTerm, None)
        return ontology
        
    def getOntologyRaMPId(self, parentTerm, childTerm):
        ontolId = None
        if self.ontolDict[parentTerm] is not None:
            ontology = self.ontolDict[parentTerm].get(childTerm, None)
            if ontology is not None:
                ontolId = ontology.ontolRampId
        return ontolId
    
    def getFullOntologyList(self):
        return self.simpleOntolList
    
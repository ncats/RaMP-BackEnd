'''
Created on Jun 3, 2021

@author: braistedjc
'''

class Ontology(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.ontolRampId = ""
        
        self.ontolParent = ""
        
        self.ontolChild = ""
        
    def __eq__(self, other):
        #print("ont check eq" + str((other.ontolParent == self.ontolParent and other.ontolChild == self.ontolChild)))          
        return other.ontolParent == self.ontolParent and other.ontolChild == self.ontolChild
    
    def __hash__(self):
        return hash(self.ontolParent + self.ontolChild)
    
    def getOntologyString(self):
        # note, the DB is set up as ID, Value, Key, so the dump is ordered in that way. Patch 20210910
        return self.ontolRampId + "\t" + self.ontolChild + "\t" + self.ontolParent + "\n"     
    
    
    
    
    
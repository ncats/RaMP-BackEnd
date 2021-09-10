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
        return self.ontolRampId + "\t" + self.ontolChild + "\t" + self.ontolParent + "\n"     
    
    
    
    
    
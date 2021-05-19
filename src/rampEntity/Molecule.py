'''
Created on Dec 7, 2020

@author: braistedjc
'''

class Molecule(object):
    '''
    Container object to hold chemical properties
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.id = ""
        
        self.source = ""
        
        self.smiles = ""
        
        self.inchiKey = ""
        
        self.inchiKeyPrefix = ""
                
        self.inchi = ""        
        
        self.mw = ""
        
        self.formula = ""
        
        self.monoisotopicMass = ""
        
        self.name = ""
        
        self.idDict = dict()
        
        self.classDict = dict()
        
        self.nameDict = dict()
        
    def toChemPropsString(self):
        """
        Ths utility method is used to format an export format for tab delimited files.
        """
        s =  self.source + "\t" + self.id + "\t" + self.smiles + "\t" + self.inchiKeyPrefix + "\t" + self.inchiKey + "\t" + self.inchi + "\t" 
        s = s + self.mw + "\t" + self.monoisotopicMass + "\t" + self.name + "\t" + self.formula+ "\n"
        return s
        
    def toSourceString(self):
        s = ""
        for idType in self.idDict:
            s = s + self.id + "\t" + idType + "\t" + self.idDict[idType] + "\n"
        return s
    
    def toCommonNameString(self): 
        return self.id + "\t" + self.name + "\n"
    
    def toClassString(self):
        s = ""        
        for classLevel in self.classDict:
            s = s + self.id + "\t" + classLevel + "\t" + self.classDict[classLevel] + "\n"            
        return s            
            
            
            
            
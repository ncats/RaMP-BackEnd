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
        
        self.names = [] 
        
        self.idDict = dict()
        
        self.classDict = dict()
        
        self.nameDict = dict()
        
        
    def addName(self, name):
        if name not in self.names:
            self.names.append(name)    
        
    def toChemPropsString(self):
        """
        Ths utility method is used to format an export format for tab delimited files.
        """
        if len(self.names) > 0:
            name = self.names[0]
        s =  self.source + "\t" + self.id + "\t" + self.smiles + "\t" + self.inchiKeyPrefix + "\t" + self.inchiKey + "\t" + self.inchi + "\t" 
        s = s + self.mw + "\t" + self.monoisotopicMass + "\t" + name + "\t" + self.formula+ "\n"
        return s
        
    def toSourceString(self):
        s = ""
        for idType in self.idDict:
            s = s + self.id + "\t" + idType + "\t" + self.idDict[idType] + "\n"
        return s
    
    def toCommonNameString(self):
        s= ""
        for name in self.names:
            self.id + "\t" + name + "\n"            
        return s
    
    def toClassString(self):
        s = ""        
        for classLevel in self.classDict:
            s = s + self.id + "\t" + classLevel + "\t" + self.classDict[classLevel] + "\n"            
        return s            
            
            
            
            
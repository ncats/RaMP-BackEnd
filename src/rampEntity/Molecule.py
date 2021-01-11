'''
Created on Dec 7, 2020

@author: braistedjc
'''

class Molecule(object):
    '''
    classdocs
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
        
        
    def toChemPropsString(self):
        s =  self.source + "\t" + self.id + "\t" + self.smiles + "\t" + self.inchiKeyPrefix + "\t" + self.inchi + "\t" + self.inchiKey + "\t" 
        s = s + self.mw + "\t" + self.monoisotopicMass + "\t" + self.name + "\t" + self.formula+ "\n"
        return s
        

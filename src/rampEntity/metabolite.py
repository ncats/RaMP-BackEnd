'''
Created on Nov 6, 2020

@author: braistedjc
'''
import pandas as pd

class metabolite:
    
    def __init__(self):

        self.sourceId = ""
        
        self.commonName = ""
        
        self.rampId = ""
        
        self.idList = list()        
        
        self.primarySource = ""
        
        self.sourceList = list()
        
        self.smiles = ""
        
        self.inchi = ""
        
        self.inchikey = ""
        
        self.synonyms = ""
        
        
        
    def __eq__(self, other):
        return len(set(self.idList).intersection(other.idList)) > 0
        

    def mergeMetabolite(self, other):
        
        
        
        
        
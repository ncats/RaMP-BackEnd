'''
Created on Jul 5, 2022

@author: braistedjc
'''

class RheaCompound(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.rampId = ""
        
        self.chebiId = ""
        
        self.name = ""
        
        self.htmlName = ""
        
        self.formula = ""
        
        self.isCofactor = 0
      
    def rheaCompoundToRecordString(self):
        s = self.chebiId + "\t" + self.name + "\t" + self.htmlName + "\t" + self.formula + "\t" + str(self.isCofactor) + "\n" 
        return s
          
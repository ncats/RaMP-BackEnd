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
        self.chebiId = ""
        
        self.name = ""
        
        self.htmlName = ""
        
        self.formula = ""
        
      
    def rheaCompoundToRecordString(self):
        s = self.chebiId + "\t" + self.name + "\t" + self.htmlName + "\t" + self.formula + "\n"
        return s
          
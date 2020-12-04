'''
Created on Nov 6, 2020

@author: braistedjc
'''

class Pathway(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.pathRampId = ""
        
        self.pathSourceId = ""
 
        self.pathSource = ""
        
        self.pathName = ""
               
        self.pathCategory = ""
        
    def __eq__(self, otherPathway):
        return self.pathRampId == otherPathway.pathRampId
       
    def printPathway(self):
        s = ""
        s = s + "rampId: " + self.pathRampId + "\n"
        s = s + "source: " + self.pathSource + "\n"
        s = s + "sourceId: " + self.pathSourceId + "\n"
        s = s + "pathwayName: " + self.pathName
        print(s)
    
    def toPathwayString(self):
        s = self.pathRampId + "\t" + self.pathSourceId + "\t" + self.pathSource + "\t" + str(self.pathCategory) + "\t" + self.pathName + "\n"
        return s
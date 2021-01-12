'''
Created on Nov 6, 2020

@author: braistedjc
'''

class Pathway(object):
    '''
    Container class for pathway information.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        # ramp pathway id
        self.pathRampId = ""
        
        # pathway source id
        self.pathSourceId = ""
 
        # pathway source
        self.pathSource = ""
        
        # pathway name
        self.pathName = ""
               
        # pathway category.
        # this is populated when data source adopts a pathway from another resource
        self.pathCategory = ""
        
    def __eq__(self, otherPathway):
        """
        Returns true if it shares a pathway ramp id.
        Currently (12/2021) we do not collapse pathways. 
        """
        return self.pathRampId == otherPathway.pathRampId
       
        
    def printPathway(self):
        """
        Utility method to print a pathway summary report to standard output.
        """
        s = ""
        s = s + "rampId: " + self.pathRampId + "\n"
        s = s + "source: " + self.pathSource + "\n"
        s = s + "sourceId: " + self.pathSourceId + "\n"
        s = s + "pathwayName: " + self.pathName
        print(s)
    
    
    def toPathwayString(self):
        """
        returns a tab-delimited string representing the pathway.
        The format is suitable for loading into the pathway table and is exported to a file.
        """
        s = self.pathRampId + "\t" + self.pathSourceId + "\t" + self.pathSource + "\t" + str(self.pathCategory) + "\t" + self.pathName + "\n"
        return s
    
    
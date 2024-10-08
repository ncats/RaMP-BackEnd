'''
Created on Jun 30, 2022

@author: braistedjc
'''

class Protein(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.rampId = ""
        
        self.uniprotAcc = ""
        
        self.secondaryAccs = []

        self.id = ""
        
        self.recName = ""
        
        self.geneName = ""
        
        self.ec = ""
        
        self.hgncSymbol = ""
        
        # status of review by uniprot.
        # TrEMBL uniprot entries are not fully curated (0)
        # SwissProt uniprot entries are completely curated (1)
        self.isReviewed = 0
        
    
    def getPrimaryRecord(self):
        s = self.uniprotAcc + "\t" + self.id + "\t" + self.geneName + "\t" + self.recName + "\t" + self.hgncSymbol + '\t' + str(self.isReviewed) + "\n"
        return s
    
    def getSecondaryToPrimaryRecord(self):
        s = ""
        for sec in self.secondaryAccs:
            s = s + sec + "\t" + self.uniprotAcc + "\n"
        return s
    
    def getPrimaryToSecondaryRecord(self):
        s = ""
        for sec in self.secondaryAccs:
            s = s + self.uniprotAcc + "\t" + sec + "\n"
        return s
    
    def printProtein(self):
        s = "ACC: " + self.uniprotAcc +"\nSecondaryACC: "
        print(str(len(self.secondaryAccs)))
        for sec in self.secondaryAccs:
            s = s + sec + ","
        s = s + "\n"
        s = s + "ID: " + self.id + "\n"
        s = s + "recName: " + self.recName + "\n"
        s = s + "geneName: " + self.geneName + "\n"
        s = s + "hgncSymbol: " + self.hgncSymbol + "\n"
        
        print(s)
        
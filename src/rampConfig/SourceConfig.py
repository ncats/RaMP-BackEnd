'''
Created on Nov 6, 2020

@author: braistedjc
'''

class SourceConfig(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor    
        '''
        self.resourceName = ""
        
        self.sourceFetchMethod = ""
        
        self.sourceURL = ""
        
        self.sourceFileName = ""
        
        self.compressType = ""

        self.extractFileName = ""
                
        self.localDir = ""
        
        self.resourceType = ""
        
    def printConfiguration(self):
        s = "resourceName: " + self.resourceName + "\n"
        s = s + "fetchMethod: " + self.sourceFetchMethod + "\n"
        s = s + "url: " + self.sourceURL + "\n"
        s = s + "sourceFile: " + self.sourceFileName + "\n"
        s = s + "compression: " + self.compressType + "\n"
        s = s + "extractFileName: " + self.extractFileName + "\n"
        s = s + "localDir: " + self.localDir + "\n"
        s = s + "resourceType: " + self.resourceType
        print(s)
        
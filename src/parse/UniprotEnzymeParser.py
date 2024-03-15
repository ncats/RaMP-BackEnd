'''
Created on Mar 15, 2024

@author: braistedjc
'''
import sys, os
from os.path import exists

from parse.MetabolomicsData import MetabolomicsData
from rampConfig.RampConfig import RampConfig
from rampConfig.SourceConfig import SourceConfig

class UniprotEnzymeParser(MetabolomicsData):
    '''
    classdocs
    '''

    def __init__(self, resConfig):
        
                # relative dir, use '../' for testing, use "" for production calls
        self.relDir = ""
        
        self.resourceConfig = resConfig
        
        self.keyFields = ['ID', "DE"]
        
        self.ecToDescription = dict()
        


    def parseEnzyme(self):

        # need to load uniprot TrEMBL and SwissProt human

        # TrEMBL
        enzymeConfig = self.resourceConfig.getConfig("expasy_enzyme_dat")
        file_proteins = enzymeConfig.sourceFileName
        proteins_url = enzymeConfig.sourceURL        
        localDir = enzymeConfig.localDir
        extractFile = enzymeConfig.extractFileName 
        remoteFile = enzymeConfig.sourceFileName


        # make the data dir if needed...
        if not exists(self.relDir + localDir):
            os.mkdir(self.relDir + localDir)

        if not exists(self.relDir + localDir + extractFile):

            self.download_files(proteins_url, self.relDir + localDir + remoteFile)
            
        else:
            print("Uniprot/Expasy enzyme.dat exists. Using cached copy.")
                
                

        print("starting to parse uniprot expasy enzyme.dat")

        self.parseEnzymeDatFile(self.relDir + localDir + extractFile)

        print("number of ec records")
        ecCount = len(self.ecToDescription)
        print(str(ecCount))

        # now add SwissProt human
  


    def parseEnzymeDatFile(self, filePath):
        enzymeDB = open(filePath, 'r+', encoding="utf-8")
        
        protein = None
        
        start = 1
        
        while True:
               
            line = enzymeDB.readline()

            if start == 1:
                ec = ""
                description = ""
         
            if len(line) == 0:
                print("line is none...")
                break    
                
            prefix = line[0:2]
            
            if prefix in self.keyFields:
                
                start = 2
                
                line = line.strip()
                
                print(line)
                
                if prefix == 'ID':
                    splitLine = line.split(" ")
                    ec = [s for s in splitLine if s][1]
                if prefix == 'DE':
                    if description == "":
                        description = line.replace('DE   ', '')
                        description = description.strip()
                        
                    
                      
            if(prefix == "//"):
                
                # check transfered entry condition
                if description.find('Transferred entry') == -1:                
                    if ec != "":
                        self.ecToDescription[ec] = description
                start = 1

        return self.ecToDescription
    

# testing code

# config = RampConfig()
# source = SourceConfig()
# 
# source.resourceName = "expasy_enzyme_dat"        
# source.sourceFetchMethod = "http"
#         
# source.sourceURL = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"
#         
# source.sourceFileName = "enzyme.dat"
#         
# source.compressType = "None"
# 
# source.extractFileName = "enzyme.dat"
#                 
# source.localDir = "C:/Users/braistedjc/Downloads/"
#         
# source.resourceType = "enzyme info"
# 
# config.configDict['expasy_enzyme_dat'] = source
# 
# uep = UniprotEnzymeParser(config)
# uep.parseEnzyme()
# 
# print(uep.ecToDescription['7.6.2.9'])
             

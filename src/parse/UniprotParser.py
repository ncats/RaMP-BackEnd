'''
Created on Jun 30, 2022

@author: braistedjc
'''
import sys, os
from os.path import exists
import gzip
import shutil
from parse.MetabolomicsData import MetabolomicsData
from rampEntity.Protein import Protein
from rampConfig.RampConfig import RampConfig

class UniprotParser(MetabolomicsData):
    '''
    classdocs
    '''

    def __init__(self, resConfig):
        '''
        Constructor
        '''
        
        # relative dir, use '../' for testing, use "" for production calls
        self.relDir = ""
        
        self.resourceConfig = resConfig
        
        self.uniprotRecords = dict()
        
        self.keyFields = ['AC', "DE", "GN", "ID", "DR"]
        
        self.parsingReviewedProteins = 0

    
    def parseUniprot(self):

        # need to load uniprot TrEMBL and SwissProt human

        # TrEMBL
        proteinConfig = self.resourceConfig.getConfig("uniprot_human")
        file_proteins = proteinConfig.sourceFileName
        proteins_url = proteinConfig.sourceURL        
        localDir = proteinConfig.localDir
        extractFile = proteinConfig.extractFileName 
        remoteFile = proteinConfig.sourceFileName


        # make the data dir if needed...
        if not exists(self.relDir + localDir):
            os.mkdir(self.relDir + localDir)

        if not exists(self.relDir + localDir + extractFile):

            self.download_files(proteins_url, self.relDir + localDir + remoteFile)
            
            with gzip.open(self.relDir + localDir + remoteFile, 'rb') as f_in:
                with open(self.relDir + localDir + extractFile, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            print("Uniprot (Human) exists. Using cached copy.")
                
                

        print("starting to parse uniprot trembl dat file")
        print(extractFile)                
                
        self.parsingReviewedProteins = 0
        self.parseUniprotFile(self.relDir + localDir + extractFile, 0)

        print("number of uniprot trembl records")
        tremblCount = len(self.uniprotRecords)
        print(str(tremblCount))


        # now add SwissProt human
        proteinConfig = self.resourceConfig.getConfig("swissprot_human")
        file_proteins = proteinConfig.sourceFileName
        proteins_url = proteinConfig.sourceURL        
        localDir = proteinConfig.localDir
        extractFile = proteinConfig.extractFileName 
        remoteFile = proteinConfig.sourceFileName


        # make the data dir if needed...
        if not exists(self.relDir + localDir):
            os.mkdir(self.relDir + localDir)

        if not exists(self.relDir + localDir + extractFile):

            self.download_files(proteins_url, self.relDir + localDir + remoteFile)
            
            with gzip.open(self.relDir + localDir + remoteFile, 'rb') as f_in:
                with open(self.relDir + localDir + extractFile, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            print("Uniprot (Human) exists. Using cached copy.")
            
            
            
        print("starting to parse uniprot swissprot dat file")
        print(extractFile)
                
        self.parsingReviewedProteins = 1
        self.parseUniprotFile(self.relDir + localDir + extractFile, 1)

        print("number of uniprot trembl PLUSE swissprot records")
        tremblCount = len(self.uniprotRecords)
        print(str(tremblCount))
        
        self.exportUniprotIntermediatFiles()
    
    
    
    
    def parseUniprotFile(self, filePath, isProteinReviewed):
        proteinDB = open(filePath, 'r+', encoding="utf-8")
        
        protein = None
        
        start = 1
        
        while True:
            
            line = proteinDB.readline()

            if len(line) == 0:
                print("line is none...")
                break    
            
            if(start == 1):
                protein = Protein()
                protein.isReviewed = isProteinReviewed
                start = 0
                
            prefix = line[0:2]
            
            if prefix in self.keyFields:
                # print("processing data...")
                self.processData(prefix, line, proteinDB, protein)
                      
            if(prefix == "//"):
                self.uniprotRecords[protein.uniprotAcc] = protein
#                 print("added record protein = "+protein.uniprotAcc+" secondaryAccCount="+str(len(protein.secondaryAccs)))
#                 for a in protein.secondaryAccs:
#                     print("secondary:"+a)
                
                #print("uniprot key: "+protein.uniprotAcc)
                protein = Protein()
                protein.isReviewed = isProteinReviewed
            
            
            
    def processData(self, prefix, line, proteinDB, protein):
        
        #protein.isReviewed = self.parsingReviewedProteins
        
        line = (line[3:(len(line)-1)]).strip()
        
        ['AC', "DE", "GN", "ID", "DR"]
        
        if(prefix == 'AC'):
            
            if(line[len(line)-1] == ';'):
                line = line[0:len(line)-1]
            
            accs =  line.split(';')
                
            for i in range(len(accs)):
                accs[i] = 'uniprot:' + accs[i].strip()
            
            if protein.uniprotAcc == "":
                protein.uniprotAcc = accs[0]
                protein.secondaryAccs = accs
                
            else:
                for acc in accs:
                    protein.secondaryAccs.append(acc)

            
        elif(prefix == 'DE'):
            line = line.strip()
            if protein.recName == "":
                if(line[0:7] == 'RecName'):
                    line = line.replace("RecName: Full=", "")
                    recName = line.split('{')[0].strip()
                    if recName[-1] == ";":
                        recName = recName[0:(len(recName)-1)]
                    protein.recName = recName
                if(line[0:7] == 'SubName'):
                    line = line.replace("SubName: Full=", "")
                    recName = line.split('{')[0].strip()
                    if recName[-1] == ";":
                        recName = recName[0:(len(recName)-1)]
                    protein.recName = recName
                    
            
        elif(prefix == 'GN'):
            # print("GN process")
            # print(line)
            if line[0:4] == 'Name':
                geneName = line.split("=")[1].strip()
                geneName = geneName.split(";")[0].strip() 
                geneName = geneName.split(" ")[0].strip() 
                             
                # geneName = geneName[0:(len(geneName)-1)]
                protein.geneName = geneName
                
        elif(prefix == 'ID'):
            id = line.split(" ")
            id = id[0].strip()
            protein.id = id
                
        elif(prefix == 'DR'):
            drLine = line.split(" ")
            if(drLine[0] == 'HGNC;'):
                symbol = drLine[2].strip()
                symbol = symbol[0:(len(symbol)-1)]
                protein.hgncSymbol = symbol


    def exportUniprotIntermediatFiles(self):
        
        dir = self.relDir + "../misc/output/uniprot_human/"
        
        if not exists(dir):
            os.mkdir(dir)
        
        recordsFile = "uniprot_records.txt"
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for acc in self.uniprotRecords:
            protein = self.uniprotRecords[acc]
            recordOut.write(protein.getPrimaryRecord())

        recordOut.close()
        
        mappingFile = "uniprot_acc_mapping.txt"
        
        mappingOut = open(dir + mappingFile, 'w', encoding="utf-8")
        for acc in self.uniprotRecords:
            protein = self.uniprotRecords[acc]
            mappingOut.write(protein.getSecondaryToPrimaryRecord())

        mappingOut.close()
            

# rConf = RampConfig()
# rConf.loadConfig("../../config/external_resource_config.txt")
#                         
# up = UniprotParser(rConf)
# up.parseUniprot()
# #up.parseUniprotFile("C:/Users/braistedjc/Desktop/Analysis/RaMP/Rhea/human_uniprot/uniprot_sprot_human.dat")
# print(str(len(up.uniprotRecords)))
# 
# uniProtDict = up.uniprotRecords
# print("\n")
# p = uniProtDict['Q9C0D3']
# 
# p.printProtein()
# print("\n")
# 
# p = uniProtDict['O43149']
# 
# p.printProtein()
# print("\n")

'''
Created on Dec 7, 2020

@author: braistedjc
'''
import sys, os
from os.path import exists
from rampEntity.Molecule import Molecule
from parse.MetabolomicsData import MetabolomicsData
import zipfile
import gzip
import shutil
import re
import pandas as pd
import numpy as np
from rampConfig.RampConfig import RampConfig
import itertools

from rampEntity.Metabolite import Metabolite

from ftplib import FTP
from bz2 import compress

class ChemWrangler(object):
    '''
    The ChemWrangler class has production and utility methods for working with molecular information
    associated with RaMP Metabolites.
    '''
    def __init__(self, resConfig):
        '''
        Constructor        
        '''        
        self.resourceConfig = resConfig
        
        self.chemLibDict = dict()


    """
    Fetch compound properties
    """    
    def fetchCompoundPropertiesFiles(self, sources):
                
        metData = MetabolomicsData()
                
        if "hmdb" in sources:        
            
            conf = self.resourceConfig.getConfig('hmdb_met_sdf')
            localDir = conf.localDir
            url = conf.sourceURL
            remoteFile = conf.sourceFileName
            extractFile = conf.extractFileName
            
            # makes the dir if needed
            metData.check_path(localDir)
            
            if not exists(localDir + extractFile):
                metData.download_files(url, localDir+remoteFile)
                        
                with zipfile.ZipFile(localDir+remoteFile,"r") as zip_ref:
                    zip_ref.extractall(localDir)
            else:
                print("HMDB SDF extists. Using cached copy of SDF.")

            
        if "chebi" in sources:
            
            conf = self.resourceConfig.getConfig('chebi_met_sdf')
            localDir = conf.localDir
            url = conf.sourceURL
            remoteFile = conf.sourceFileName
            extractFile = conf.extractFileName

            # makes the dir if needed
            metData.check_path(localDir)   

            if not exists(localDir + extractFile):

                metData.download_files(url, localDir+remoteFile)
            
                with gzip.open(localDir+remoteFile, 'rb') as f_in:
                    with open(localDir+extractFile, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            else:
                print("Chebi SDF extists. Using cached copy of SDF.")
                
                
        if "lipidmaps" in sources:        
            
            conf = self.resourceConfig.getConfig('lipidmaps_met')
            localDir = conf.localDir
            url = conf.sourceURL
            remoteFile = conf.sourceFileName
            extractFile = conf.extractFileName
            
            # makes the dir if needed
            metData.check_path(localDir)
            
            if not exists(localDir + extractFile):
                metData.download_files(url, localDir+remoteFile)
                        
                with zipfile.ZipFile(localDir+remoteFile,"r") as zip_ref:
                    zip_ref.extractall(localDir)
            else:
                print("LipidMaps SDF extists. Using cached copy of SDF.")                
    
#     """
#     Fetch compound properties
#     """    
#     def fetchFile(self, source, outputDir, fetchURL, remoteFile, compressionMode = 'gzip'):
#         
#         print(os.getcwd())
#         
#         metData = MetabolomicsData()
#         url = fetchURL
#         targetDest = outputDir+"/"+remoteFile
#         metData.download_files(url, targetDest)
#         
#         if compressionMode == 'zip':
#             with zipfile.ZipFile(targetDest,"r") as zip_ref:
#                 zip_ref.extractall(outputDir)
#         
#         if compressionMode == 'gzip':
#             extractTarget = targetDest.replace(".gz", '')      
#             with gzip.open(targetDest, 'rb') as f_in:
#                 with open(extractTarget, 'wb') as f_out:
#                     shutil.copyfileobj(f_in, f_out)
                
                
    
    
        
    def readHMDBSDF(self, source, filePath):
        """
        Utility method to read HMDB SDF file for specific data types to populate Molecule objects.
        """
        print("HMDB SDF")
        
        i = 0
        sdfDB = open(filePath, 'r+', encoding="utf-8")
        molDict = dict()
        
        mol = Molecule()
        mol.source = source
            
        while True:
            line = sdfDB.readline()

            if len(line) == 0:
                print("line is none...")
                break
            
            line = line.strip()            
            if line == '$$$$':
                i = i + 1
               # print("processing structure " + str(i))
                molDict[mol.id] = mol
                mol = Molecule()
                mol.source = source
            if line == '> <DATABASE_ID>':
                mol.id = "hmdb:" + sdfDB.readline().strip()
            if line == '> <SMILES>':
                mol.smiles = sdfDB.readline().strip()
            if line == '> <INCHI_KEY>':
                mol.inchiKey = sdfDB.readline().strip()
                mol.inchiKeyPrefix = mol.inchiKey.split("-")[0]                
            if line == '> <INCHI_IDENTIFIER>':
                mol.inchi = sdfDB.readline().strip()                                          
            if line == '> <MOLECULAR_WEIGHT>':
                mol.mw = sdfDB.readline().strip()
            if line == '> <EXACT_MASS>':
                mol.monoisotopicMass = sdfDB.readline().strip()
            if line == '> <GENERIC_NAME>':
                mol.addName(sdfDB.readline().strip())
            if line == '> <FORMULA>':
                mol.formula = sdfDB.readline().strip()

        self.chemLibDict[source] = molDict
        print("Finished HMDB chemprops: size="+str(len(molDict)))
       
        
    def readChebiSDF(self, source, filePath):
        """
        Utility method to read Chebi SDF file for specific data types to populate Molecule objects.
        """        
        print("ChEBI SDF")

        i = 0
        sdfDB = open(filePath, 'r+', encoding="utf-8")
        molDict = dict()
        mol = Molecule()
        mol.source = source
        
        while True:
            line = sdfDB.readline()

            if len(line) == 0:
                print("line is none...")
                break

            line = line.strip()            
            if line == '$$$$':
                i = i + 1
               # print("processing structure " + str(i))
                molDict[mol.id] = mol
                mol = Molecule()
                mol.source = source
            if line == '> <ChEBI ID>':
                mol.id = sdfDB.readline().strip()
                mol.id = mol.id.replace("CHEBI", "chebi")
            if line == '> <SMILES>':
                mol.smiles = sdfDB.readline().strip()
            if line == '> <InChIKey>':
                mol.inchiKey = sdfDB.readline().strip()
                mol.inchiKeyPrefix = mol.inchiKey.split("-")[0]                
            if line == '> <InChI>':
                mol.inchi = sdfDB.readline().strip()                                      
            if line == '> <MASS>':
                mol.mw = sdfDB.readline().strip()
            if line == '> <Monoisotopic Mass>':
                mol.monoisotopicMass = sdfDB.readline().strip()
            if line == '> <ChEBI Name>':
                mol.addName(sdfDB.readline().strip())
            if line == '> <Formulae>':
                mol.formula = sdfDB.readline().strip()

        self.chemLibDict[source] = molDict
        print("Finished ChEBI chemprops: size="+str(len(molDict)))


    def readKEGGCompound(self, source, filePath):
        """
        Utility method to read KEGG 2021 Compound file for specific data types to populate Molecule objects.
        """
        i = 0
        sdfDB = open(filePath, 'r+', encoding="utf-8")
        molDict = dict()
        mol = Molecule()
        mol.source = source
            
        linekey = ""
        while True:
            line = sdfDB.readline()

            lineToks = re.split("\s{1,}", line)
            if(len(lineToks) > 1):
                linekey = lineToks[0]
                # print(linekey + ".... linekey")

            if len(line) == 0:
                print("line is none...")
                break

            line = line.strip()            
            if line == '///':
                i = i + 1
                # print("processing structure " + str(i))
                molDict[mol.id] = mol
                mol = Molecule()
                mol.source = source
            if linekey == 'ENTRY':
                mol.id = "kegg:"+lineToks[1].strip()                           
            if linekey == 'MOL_WEIGHT':
                mol.mw = lineToks[1].strip()
            if linekey == 'EXACT_MASS':
                mol.monoisotopicMass = lineToks[1].strip()
            if linekey == 'NAME':
                name = lineToks[1].strip()
                mol.addName(name[:-1])
            if linekey == 'FORMULA':
                mol.formula = lineToks[1].strip()

        # add to to full dictionary
        self.chemLibDict[source] = molDict


    def readLipidMapsSDF(self, source, filePath):
        """
        Utility method to read Chebi SDF file for specific data types to populate Molecule objects.
        """        
        
        print("Lipidmaps SDF")

        idDict = {
        "> <PUBCHEM_CID>":"pubchem",
        "> <KEGG_ID>":"kegg",
        "> <HMDB_ID>":"hmdb",
        "> <CHEBI_ID>":"chebi",
        "> <SWISSLIPIDS_ID>":"swisslipids",
        "> <LIPIDBANK_ID>":"lipidbank",
        "> <PLANTFA_ID>":"plantfa"
        }
        
        classDict = { "> <CATEGORY>":"LipidMaps_category", "> <MAIN_CLASS>":"LipidMaps_main_class", "> <SUB_CLASS>":"LipidMaps_sub_class"}

        i = 0
        sdfDB = open(filePath, 'r+', encoding = 'utf-8')

        molDict = dict()
        mol = Molecule()
        mol.source = source

        while True:
            line = sdfDB.readline()

            if len(line) == 0:
                print("line is none...")
                break

            line = line.strip()            
            if line == '$$$$':
                i = i + 1
               # print("processing structure " + str(i))
                molDict[mol.id] = mol
                mol = Molecule()
                mol.source = source
            if line == '> <LM_ID>':
                mol.id = sdfDB.readline().strip()
                mol.id = "LIPIDMAPS:" + mol.id 
            if line == '> <SMILES>':
                mol.smiles = sdfDB.readline().strip()
            if line == '> <INCHI_KEY>':
                mol.inchiKey = sdfDB.readline().strip()
                mol.inchiKeyPrefix = mol.inchiKey.split("-")[0]                
            if line == '> <INCHI>':
                mol.inchi = sdfDB.readline().strip()                                      
            if line == '> <MASS>':
                mol.mw = sdfDB.readline().strip()
            if line == '> <EXACT_MASS>':
                mol.monoisotopicMass = sdfDB.readline().strip()
            if line == '> <NAME>':
                name = sdfDB.readline().strip()
                mol.addName(name)
                mol.nameDict[source] = name
            if line == '> <COMMON_NAME>':
                name = sdfDB.readline().strip()
                mol.addName(name)
                mol.nameDict[source] = name               
            if line == '> <SYSTEMATIC_NAME>':
                name = sdfDB.readline().strip()
                mol.addName(name)
                mol.nameDict[source] = name
            if line == '> <SYNONYMS>':
                vals = line.split(";")
                for val in vals:
                    mol.addName(val.strip())
                    mol.nameDict[source] = val.strip()
            if line == '> <FORMULA>':
                mol.formula = sdfDB.readline().strip()
            if line in idDict:
                mol.idDict[idDict[line]] = idDict[line]+":"+sdfDB.readline().strip()
            if line in classDict:
                mol.classDict[classDict[line]] = sdfDB.readline().strip()

        self.chemLibDict[source] = molDict
        print("Finished Lipidmaps chemprops: size="+str(len(molDict)))

            
    def readPubchemTabIdMiInchikey(self, source, file):
        """
        Utility method to read a small subset of pubchem to evaluate harmonization
        This includes pubchem id, inchikey, monoisotopic mass and formula (12/2021)
        """
        sdfDB = open(file, 'r+', encoding="utf-8")
        
        molDict = dict()
        mol = Molecule()
        mol.source = source
        
        while True:
            line = sdfDB.readline()

            if line is None or len(line) == 0:
                break
            
            mol = Molecule()
            mol.source = source
        
            lineToks = line.split("\t")
            mol.id = lineToks[0].strip()
            mol.monoisotopicMass = lineToks[1].strip()
            mol.inchiKey = lineToks[2].strip()
            mol.inchiKeyPrefix = mol.inchiKey.split("-")[0]
        
            molDict[mol.id] = mol
            
        sdfDB.close()
        
        self.chemLibDict[source] = molDict



    def readSDF(self, source, file):
        """
        Reads SDF for specified data source and file
        """
        if source == 'hmdb':
            self.readHMDBSDF(source, file)
        if source == 'chebi':
            self.readChebiSDF(source, file)
        if source == 'kegg':
            self.readKEGGCompound(source, file)
        if source == 'pubchem':
            self.readPubchemTabIdMiInchikey(source, file)
        if source == 'lipidmaps':
            self.readLipidMapsSDF(source, file)          
            

    def loadRampChemRecords(self, sources):
        """
        Populates the chemical record list for the list of passed sources.
        """
        
        print("Fetching Compound Properties SDF Files")

        
        self.fetchCompoundPropertiesFiles(sources)
        
        print("Finished fetching compound property files.")
        print("Processing compound property source files. Building molecule list.")
        
        for source in sources:
            if source == 'hmdb':
                sdfConfig = self.resourceConfig.getConfig('hmdb_met_sdf')
                file = sdfConfig.localDir + sdfConfig.extractFileName                
                self.readSDF('hmdb', file)
            if source == 'chebi':
                sdfConfig = self.resourceConfig.getConfig('chebi_met_sdf')
                file = sdfConfig.localDir + sdfConfig.extractFileName                
                self.readSDF('chebi', file)
            if source == 'kegg':
                file = "../../misc/data/chemprops/kegg_compound.txt"
                self.readSDF('kegg', file)
            if source == 'pubchem':
                file = "../../misc/data/chemprops/pubchem_id_mi_inchikey_issue_set.txt"
                self.readSDF('pubchem', file) 
            if source == 'lipidmaps':
                sdfConfig = self.resourceConfig.getConfig('lipidmaps_met')
                file = sdfConfig.localDir + sdfConfig.extractFileName 
                self.readSDF("lipidmaps", file)      
        
        print("Finished processing compound property source files.")
   
           
         
    def getChemSourceRecords(self):
        """
        Returns the chemical library dictionary. Molecule objects associated with input data sources.
        """
        return self.chemLibDict
       

                               
    def evaluateMolecularCollisions(self, source):
        '''
        This is a utility method developed to report on constituative isomers
        This simply reviews compuonds from a source, collecting inchikey prefixes that map to compound formula keys.
        The constructed has will have formula keys and a collection of inchikey prefixes. 
        A second hash holds collections of source ids for each formula.
        A file is output that reports on each distinct formula that has a least one pair of isomers.
        '''        
        sources = list()
        sources.append(source)
        self.loadRampChemRecords(sources)
        
        # two ways to evaluate Mi collisions
        # first look for the same formula, different inchi prefix
        # these are the constituent isomers
        
        # build has hash of formula to a list of inchi prefixes
        chemList = list()
        formulaToInchiPrefixListDict = dict()
        formulaToIdDict = dict()
        distinctInchiPrefixList = list()

        currList = None
                
        chemList = self.chemLibDict[source]
        compCount = 0
        for key in chemList:
            mol = chemList.get(key)
            compCount = compCount + 1
            if mol.inchiKeyPrefix is not None and mol.formula is not None:
                
                if(compCount % 10000 == 0):
                    print(mol.formula)
                    
                currList = formulaToInchiPrefixListDict.get(mol.formula, None)
                if currList is not None:
                    if mol.inchiKeyPrefix not in currList:
                        currList.append(mol.inchiKeyPrefix)
                        if mol.inchiKeyPrefix not in distinctInchiPrefixList:
                            distinctInchiPrefixList.append(mol.inchiKeyPrefix)
                else:
                    currList = list()
                    currList.append(mol.inchiKeyPrefix)
                    formulaToInchiPrefixListDict[mol.formula] = currList
            
            currList = formulaToIdDict.get(mol.formula, None)
            if currList is not None:
                if mol.id not in currList:
                    currList.append(mol.id)
            else:
                currList = list()
                currList.append(mol.id)
                formulaToIdDict[mol.formula] = currList
                                      
        file  = open("../../misc/data/chemprops/constituativeIsomerList.txt", "w+", encoding='utf-8')
         
        for key in formulaToInchiPrefixListDict:
            if len(formulaToInchiPrefixListDict[key]) > 1:
                s = key + "\t" + str(len(formulaToInchiPrefixListDict[key])) + "\t" 
                for inchiPre in formulaToInchiPrefixListDict[key]:
                    s = s + inchiPre + ", "
                s = s[:-2]
                s = s + "\t"
                for id in formulaToIdDict[key]:
                    s = s + id + ", "
                s = s[:-2]
                s = s + "\n"
                file.write(s)
          
        file.close()

    
    def castMoleculesToMetabolites(self, moleculeDict, source):
        
        for molId in moleculeDict:
            mol = moleculeDict[molId]
            met = Metabolite()
            
            # first ids 
            met.sourceId = mol.sourceId
            met.idDict[source] = mol.sourceId
            for idType in mol.idDict:
                met.idDict[idType] = mol.idDict[idType]
            
            # now common name
            met.commonNameDict[source] = moleculeDict[source]
            
            # now add source
            met.addSource(source)
    
            # now add classes
            for classLevel in mol.classDict:
                met.addMetClass(source, met.sourceId, classLevel, mol.classDict)
                
        
    def getCatalyticDistances(self, catalyzesFile):
        # prototype method exploration
        # collecting metablolite distances based on catalyzed table
            
        print("computing catalytic distances")
        
        catMat = pd.read_csv(catalyzesFile, sep='\t', header=None)
        catMat.columns = ["compound","protein"]

        lowestDist = dict()
        print("g to c shape")
        print(catMat.shape)
        #print(catMat)
        
        #print(catMat["protein"].to_list())
        #print(catMat["compound"].to_list())
              
        # build base connectivity
        # gene by compound
        conMat = pd.DataFrame(index=np.unique(catMat["compound"].to_list()), columns=np.unique(catMat["protein"].to_list()))
        print("connectivity shape")
        print(conMat.shape)
        
        conMat = conMat.fillna(0)
#         
        print("connectivity shape")
        print(conMat.shape)
        
        #print(conMat)
#
        #print(conMat.index)
        #print(conMat.columns)
        # base connetivity
#         for i, row in catMat.iterrows():
#             if(i % 10000 == 0):
#                 print(str(i))
#             conMat.loc[row[0], row[1]] = 1
        
         
        rowCnt = catMat.shape[0] - 1    
        for i in range(0, rowCnt):
            if(i % 10000 == 0):
                print(str(i))
            conMat.loc[catMat.iloc[i,0], catMat.iloc[i,1]] = 1
            
                

        print("connectivity shape")
        print(conMat.shape)
#       
        #let's transpose, then collect tuple sets
        conMat = conMat.transpose()
        adjMat = list()
        overallCount = 0
        for ind, row in conMat.iterrows():
            cList = list()
            for i in range(0,len(row)-1):
                overallCount = overallCount + 1
                if(overallCount % 1000000 == 0):
                    print(str(overallCount))
                if row[i] == 1:
                    cList.append(conMat.columns[i])
            adjMat.append(itertools.combinations(cList,2))        
            
        print("length adjMat = " + str(len(adjMat)))
        
        adjMat = pd.DataFrame(adjMat)
        
        print("adjMat shape = " + str(adjMat.shape))
        adjMat = adjMat.drop_duplicates()
        
        print("adjMat shape = " + str(adjMat.shape))
                
        
        # level 1 connectivity, adjacency matrix
#         #adjMat = conMat.dot(conMat.transpose())
#         adjMat = np.matmul(conMat, conMat.transpose())
#         
#         # have adjacency matrix 
#         # cast back into a sparse representation.
#         
#         
#         np.fill_diagonal(adjMat.values, 0)
#         
#         print("adjMat shape")
#         print(adjMat.shape)
#         
#         self.collectNonZeroEntries(adjMat, lowestDist, 1)
#         
#         levelMat = adjMat
#         
#         for level in range(2,10):
#             levelMat = levelMat.dot(adjMat)
#             np.fill_diagonal(levelMat.values, 0)
#             self.collectNonZeroEntries(levelMat, lowestDist, level)
        
        # set self distance to 0
        #conMat2.values[[np.arange(conMat2.shape[0])]*2] = 0
        
        
        
#         print("one step")
#         print(conMat2)
# 
#         conMat3 = conMat2.dot(conMat2)
#         #conMat3 = conMat2.dot(conMat)
#         np.fill_diagonal(conMat3.values, 0)
# 
#         self.collectNonZeroEntries(conMat3, lowestDist, 2)
#         
#         print("two steps")        
#         print(conMat3)
# 
#         conMat4 = conMat3.dot(conMat2)
#         np.fill_diagonal(conMat4.values, 0)
#         
#         self.collectNonZeroEntries(conMat4, lowestDist, 3)
# 
#         print("three steps")
#         print(conMat4)
#         
#         conMat5 = conMat4.dot(conMat2)
#         np.fill_diagonal(conMat5.values, 0)
#         
#         print("four steps")
#         print(conMat5)
#         
#         self.collectNonZeroEntries(conMat5, lowestDist, 4)
# 
#         print(conMat5.values.nonzero()[1])
        
#         for key in lowestDist:
#             print(key + " " + str(lowestDist[key]))

    def collectNonZeroEntries(self, matrix, lowestDist, level):
        
        newPath = False        
        nonzero = matrix.values.nonzero()
        print("heres the non zero")
        print(nonzero)
        for i in range(0,len(nonzero[0])-1):
            print(i)
            comp1 = matrix.index[nonzero[0][i]]
            comp2 = matrix.columns[nonzero[1][i]]
            compKey = comp1 + ":" + comp2
            compKey2 = comp2 + ":" + comp1            
            if compKey not in lowestDist and compKey2 not in lowestDist:                
                lowestDist[compKey] = level
                newPath = True
        
        return newPath        

# resourceConfFile = "../../config/external_resource_config.txt" 
# resourceConf = RampConfig()
# cw = ChemWrangler(resourceConf)
# 
# #cw.getCatalyticDistances("C:/Tools/git_projects/ramp/RaMP-BackEnd/src/testCatal.txt")
# cw.getCatalyticDistances("C:/Tools/git_projects/ramp/RaMP-BackEnd/misc/sql/catalyzes.txt")

#cw.loadRampChemRecords(["hmdb","chebi","lipidmaps"])

# Test    
# resourceConfFile = "../../config/external_resource_config.txt" 
# resourceConf = RampConfig()
# resourceConf.loadConfig(resourceConfFile)
# cw = ChemWrangler(resourceConf)
# sources = ["hmdb", "chebi", "lipidmaps"]
# cw.fetchCompoundPropertiesFiles(sources) 
# cw.loadRampChemRecords(sources)
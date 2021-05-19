import urllib.request
#from xml.etree.ElementTree import ElementTree
from lxml import etree as ET
import zipfile
import time
import os
from parse.MetabolomicsData import MetabolomicsData
import pandas as pd

from chemprop.ChemWrangler import ChemWrangler


class lipidmapsChemData(MetabolomicsData):
    '''This class parses lipidmaps SDF file to capture lipidmaps structures specifically for the following tables
    source, met_classes and chem_props. The goal here is to increase chemical knowledge and expand ids.
    This data source does not conform to typical MetabolimicsData but is better suited specifically as a collection of Molecule entries.
    '''
    
    def __init__(self):

        self.moleculeList = list()
        
        self.sourceName = "LIPIDMAPS"
        
        self.sourceFile = "../../misc/data/chemprops/lipid_maps/structures.sdf"
           
    
    def parseLipidMaps(self, writeToFile=False):
        chemist = ChemWrangler()
        chemist.fetchFile("LIPIDMAPS", "../../misc/data/chemprops/lipid_maps/", "https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip", "LMSD.sdf.zip", "zip")
        
        chemist.readLipidMapsSDF(self.sourceName, self.sourceFile)
        lipidMapMolecules = chemist.chemLibDict[self.sourceName]
        
        if writeToFile:
            self.write_myself_files('lipidmaps')
            self.writeFiles(lipidMapMolecules)
    
    def writeFiles(self, molDict):
        classFile = "lipidmapsmetaboliteClass.txt"
        metIdFile = "lipidmapsmetaboliteIDDictionary.txt"
        commonNameFile = "lipidmapsmetaboliteCommonName.txt"
        
        lipidMapsOutputDir = "../../misc/output/lipidmaps/"
        
        try:
            os.makedirs(lipidMapsOutputDir)
        except FileExistsError:
            # directory already exists
            pass
        
        classFileHandle  = open(lipidMapsOutputDir+classFile, "w+", encoding='utf-8')

        for id in molDict :
            classFileHandle.write(molDict[id].toClassString())
        
        classFileHandle.close()

        idFileHandle  = open(lipidMapsOutputDir+metIdFile, "w+", encoding='utf-8')

        for id in molDict :
            idFileHandle.write(molDict[id].toSourceString())
        
        idFileHandle.close()

        commonNameFileHandle  = open(lipidMapsOutputDir+commonNameFile, "w+", encoding='utf-8')

        for id in molDict :
            commonNameFileHandle.write(molDict[id].toCommonNameString())
        
        commonNameFileHandle.close()
        
        
    def getEverything(self, writeToFile = False):
        self.parseLipidMaps(writeToFile)    
    
#     print(len(lipidMapMolecules))
#     
#     mol = lipidMapMolecules["LIPIDMAPS:LMGP10070002"]
#     
#     print(mol.toChemPropsString())
#         
#     print(mol.toSourceString())
#     
#     print(mol.toCommonNameString())
#     
#     print(mol.toClassString())
    
    
        
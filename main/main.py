import sys
sys.path.append('../src')
from rampConfig.RampConfig import RampConfig
from parse.wikipathwayRDF import WikipathwaysRDF
from parse.hmdbData import hmdbData
from parse.reactomeData import reactomeData
from parse.lipidmapsChemData import lipidmapsChemData
from parse.KeggData import KeggData
from parse.RheaParser import RheaParser
from util.EntityBuilder import EntityBuilder
from getStatistics import getStatistics
from writeToSQL import writeToSQL
import os
import time

class Main():

    def runEverything(self, resourceConfigFile, getDatabaseFiles = True):

        start = time.time()

        sql = writeToSQL()
        
        # build the ramp resource config
        resourceConf = RampConfig()
        resourceConf.loadConfig(resourceConfigFile)
        
        stat = getStatistics()
        hmdb = hmdbData(resourceConf)
        wikipathways = WikipathwaysRDF(resourceConf)
        reactome = reactomeData(resourceConf)
        kegg = KeggData()
        lipidmaps = lipidmapsChemData(resourceConf)
        rhea = RheaParser(resourceConf)
        
        # works based on your computer, setup working directory
        os.chdir('../main/')

        #kegg.getEverything(False)
        #print("KEGG Wonder")
        print("Getting hmdb...")
        hmdb.getEverything(True)
        print("Getting wiki...")
        wikipathways.getEverything(True)
        print("Getting reactome...")
        reactome.getEverything(True)
        
        # This parses and writes lipid maps
        # sql write will be handled by EntityBuilder
        print("Getting LipidMaps...")
        lipidmaps.getEverything(True)
        
        print("Getting Rhea info...")
        rhea.processRhea()
        
        # constructs the entity builder
        builder = EntityBuilder(resourceConf)
        
        # performs a full build of entities for loading
        # the input are files in /misc/output
        # the result are files for DB loading in /misc/sql
          
        builder.fullBuild()

        print(time.time() - start)
        

        # Database loading is handled as a separate, un-coupled step.
            
resourceConfFile = "../config/external_resource_config.txt"                
main = Main()
main.runEverything(resourceConfigFile = resourceConfFile)







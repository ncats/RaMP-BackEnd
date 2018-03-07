from wikipathwaysData import wikipathwaysData
from writeToSQL import writeToSQL
from IDconversion import IDconversion
from getStatistics import getStatistics
import time
import csv
import unittest
import random
import time
class TestWikipathwaysMain(unittest.TestCase):

    def testMain(self):
        
        sql = writeToSQL()
        idconvert = IDconversion()
        stat = getStatistics()
        
        wikipathways = wikipathwaysData()
        wikipathways.getDatabaseFiles()
        
        
        
        
        wikipathways.getEverything()
        r1 =random.choice(list(wikipathways.geneInfoDictionary.keys()))
        r2 =random.choice(list(wikipathways.geneInfoDictionary.keys()))
        r3 =random.choice(list(wikipathways.geneInfoDictionary.keys()))
        print(wikipathways.geneInfoDictionary[r1])
        print(wikipathways.geneInfoDictionary[r2])
        print(wikipathways.geneInfoDictionary[r3])
        print(wikipathways.geneInfoDictionary['ENSG00000139977'])
        print(wikipathways.geneInfoDictionary["path:hsa04530"])
        #time.sleep(10)
        #time.sleep(30)
        wikipathways.getCommonNameForChebi()
        
        #idconvert.GeneConvert(wikipathways.geneInfoDictionary, "wiki")
        wikipathways.write_myself_files('wiki')
        wikicompoundnum = sql.createRampCompoundID(wikipathways.metaboliteIDDictionary, "wiki", 0)
        wikigenenum = sql.createRampGeneID(wikipathways.geneInfoDictionary, "wiki", 0)
        
        print("Write to file...")
        wikipathwaysnumbers = sql.write(
                 wikipathways.metaboliteCommonName,
                 wikipathways.pathwayDictionary, 
                 wikipathways.pathwayCategory,
                 wikipathways.metabolitesWithPathwaysDictionary,
                 wikipathways.metabolitesWithSynonymsDictionary,
                 wikipathways.metaboliteIDDictionary,
                 wikipathways.pathwaysWithGenesDictionary,
                 wikipathways.metabolitesLinkedToGenes,
                 wikipathways.geneInfoDictionary,
                 wikipathways.biofluidLocation,
                 wikipathways.biofluid,
                 wikipathways.cellularLocation,
                 wikipathways.cellular,
                 wikipathways.pathwayOntology,
                 wikipathways.exoEndoDictionary,
                 wikipathways.exoEndo,
                 wikipathways.tissueLocation,
                 wikipathways.tissue,
                 "wiki",
                 0, 0)
        
        print("Pathways number is " + str(len(wikipathways.pathwayDictionary)))
        print("metabolites number is " + str(len(wikipathways.metaboliteIDDictionary)))
        print('genes number is '+ str(len(wikipathways.geneInfoDictionary)))
        
        print("Compound:") 
        stat.analyteOverlaps(sql.rampCompoundIdInWhichDatabases, sql.rampCompoundIDdictionary, "Compound")
        print("\n")
        print("Gene:") 
        stat.analyteOverlaps(sql.rampGeneIdInWhichDatabases, sql.rampGeneIDdictionary, "Gene")
        
        
        

if __name__ == "__main__":
    unittest.main()

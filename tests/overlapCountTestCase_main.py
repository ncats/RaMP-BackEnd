from hmdbData import hmdbData
from writeToSQL import writeToSQL
from getStatistics import getStatistics
from IDconversion import IDconversion
from wikipathwaysData import wikipathwaysData
from reactomeData import reactomeData
from KeggData import KeggData
import time
import csv
import unittest
import os

class TestHMDBMain(unittest.TestCase):
    def testMain(self):
        sql = writeToSQL()
        hmdb = hmdbData()
        hmdb.getDatabaseFiles()
        idconvert = IDconversion()
        stat = getStatistics()
        wiki = wikipathwaysData()
        react = reactomeData()
        kegg = KeggData()
        print('Running overlap plot test case...')
        
        hmdb.getGenes()
        hmdb.getPathwaysLinkedToGene()
        
        wiki.getEverything()
        wiki.getCommonNameForChebi()
        react.getGenes()
        react.getCommonNameFromUniprot()
        kegg.getPathways()
        kegg.getMetabolites()
        kegg.getGenes()
        kegg.getGeneInfo()
        kegg.getPathwayLinkedToGene()
        #idconvert.GeneConvert(hmdb.geneInfoDictionary, "hmdb")
        
        hmdbgenenum = sql.createRampGeneID(hmdb.geneInfoDictionary, "hmdb", 0)
        keggnum = sql.createRampGeneID(kegg.geneInfoDictionary,'kegg',hmdbgenenum)
        wikinum = sql.createRampGeneID(wiki.geneInfoDictionary,'wiki',keggnum)
        reactnum = sql.createRampGeneID(react.geneInfoDictionary,'reactome',wikinum)
        
        
        stat.analyteOverlaps(sql.rampGeneIdInWhichDatabases, sql.rampGeneIDdictionary, 'Gene')   
        
        #sql.writeIdInWhichdatabase()
        # print("Compound:")
        # stat.analyteOverlaps(sql.rampCompoundIdInWhichDatabases, sql.rampCompoundIDdictionary, "Compound")

        # print("\n")
        # print("Gene:")
        # stat.analyteOverlaps(sql.rampGeneIdInWhichDatabases, sql.rampGeneIDdictionary, "Gene")

        # stat.databaseContent(hmdb.pathwayDictionary,
        #         hmdb.pathwayCategory,
        #         hmdb.metabolitesWithPathwaysDictionary,
        #         hmdb.metabolitesWithSynonymsDictionary,
        #         hmdb.metaboliteIDDictionary,
        #         hmdb.pathwaysWithGenesDictionary,
        #         hmdb.geneInfoDictionary,
        #         hmdb.biofluidLocation,
        #         hmdb.biofluid,
        #         hmdb.cellularLocation,
        #         hmdb.cellular,
        #         hmdb.pathwayOntology,
        #         hmdb.exoEndoDictionary,
        #         "hmdb")
        
if __name__ == "__main__":
    unittest.main()
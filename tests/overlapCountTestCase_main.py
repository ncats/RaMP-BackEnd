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

      
        hmdb.write_myself_files('hmdb')
        wiki.write_myself_files('wiki')
        
        hmdbgenenum = sql.createRampGeneID(hmdb.geneInfoDictionary, "hmdb", 0)
        wikinum = sql.createRampGeneID(wiki.geneInfoDictionary,'wiki',hmdbgenenum)
        reactnum = sql.createRampGeneID(react.geneInfoDictionary,'reactome',wikinum)
        keggnum = sql.createRampGeneID(kegg.geneInfoDictionary,'kegg',reactnum)
        
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
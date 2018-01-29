from hmdbData import hmdbData
from writeToSQL import writeToSQL
from getStatistics import getStatistics
from IDconversion import IDconversion

import time
import csv
import unittest
import os

class TestHMDBMain(unittest.TestCase):
    def testMain(self):
        sql = writeToSQL()
        hmdb = hmdbData()
        print(os.getcwd())
        # If does not have database file
        a= hmdb.getSMPDB_Category()
        print(a[0:5])
        print(len(a))
        hmdb.getDatabaseFiles()
        idconvert = IDconversion()
        stat = getStatistics()
        
        #print(hmdb.pathwaysWithGenesDictionary['SMP00172'])
        #print(len(hmdb.pathwaysWithGenesDictionary['SMP00172']))
        #print(hmdb.pathwaysWithGenesDictionary['00240'])
        #print(len(hmdb.pathwaysWithGenesDictionary['00240']))
        
        
       
        print("Getting HMDB Metabolites...")
        tree = hmdb.getMetaboliteOtherIDs()
        print("Getting HMDB pathways and synonyms...")
        hmdb.getPathwaysandSynonyms(tree)
        
        
        print('How many pathways relationship ...')
        #print(str(len(hmdb.metabolitesWithPathwaysDictionary)))
        print("Getting HMDB genes...")
        
        
        hmdb.getGenes(tree)

        print("Getting HMDB biofluid and cellular locations...")
        hmdb.getBiofluidCellularLocationDisease(tree)
        
        hmdb.getPathwaysLinkedToGene()
        print("Writing to files...")
        print("output each metabolites has how many pathways ...")
        

        #idconvert.GeneConvert(hmdb.geneInfoDictionary, "hmdb")

        print("hmdb compounds...")
        #sql.checkForWithinDatabaseDuplicatesCompound(hmdb.metaboliteIDDictionary, "hmdb")
        print("hmdb genes...")
        #sql.checkForWithinDatabaseDuplicatesGene(hmdb.geneInfoDictionary, "hmdb")
        hmdb.write_myself_files('hmdb')
        hmdbcompoundnum = sql.createRampCompoundID(hmdb.metaboliteIDDictionary, "hmdb", 0)
        hmdbgenenum = sql.createRampGeneID(hmdb.geneInfoDictionary, "hmdb", 0)
        
        
           
        sql.write(hmdb.metaboliteCommonName,
                  hmdb.pathwayDictionary,
                  hmdb.pathwayCategory,
                  hmdb.metabolitesWithPathwaysDictionary,
                  hmdb.metabolitesWithSynonymsDictionary,
                  hmdb.metaboliteIDDictionary,
                  hmdb.pathwaysWithGenesDictionary,
                  hmdb.metabolitesLinkedToGenes,
                  hmdb.geneInfoDictionary,
                  hmdb.biofluidLocation,
                  hmdb.biofluid,
                  hmdb.cellularLocation,
                  hmdb.cellular,
                  hmdb.pathwayOntology,
                  hmdb.exoEndoDictionary,
                  hmdb.exoEndo,
                  hmdb.tissueLocation,
                  hmdb.tissue,
                  "hmdb",
                  0, 0)
        
        print('Compound number is ' + str(len(hmdb.metaboliteIDDictionary)))
        print('Gene number is ' + str(len(hmdb.geneInfoDictionary)))
        print('Pathway number is ' + str(len(hmdb.pathwayDictionary)))
        print('Pathway that have Gene with it ' + str(len(hmdb.pathwaysWithGenesDictionary)))
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
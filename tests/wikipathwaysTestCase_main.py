from wikipathwaysData import wikipathwaysData
from writeToSQL import writeToSQL
from IDconversion import IDconversion
from getStatistics import getStatistics
import time
import csv
import unittest

class TestWikipathwaysMain(unittest.TestCase):

    def testMain(self):
        
        sql = writeToSQL()
        idconvert = IDconversion()
        stat = getStatistics()
        
        wikipathways = wikipathwaysData()
        #wikipathways.getDatabaseFiles()
        wikipathways.getEverything()
        print(wikipathways.setOfType)
        wikipathways.getCommonNameForChebi()
        file = open("../misc/output/wikiGenesID.txt","wb")
        for key in wikipathways.geneInfoDictionary:
            file.write(key.encode("utf-8") +b"\n")
        file.close()
        print(len(wikipathways.metaboliteIDDictionary))
        print(len(wikipathways.metabolitesWithPathwaysDictionary))
        print(len(wikipathways.geneInfoDictionary))
        print(len(wikipathways.pathwaysWithGenesDictionary))
        time.sleep(10)
        '''
        outputToFile = open("../misc/output/figure/WIKImetaWithPaths.csv","w",newline="")
        fieldname = ["metabolites","numOfPathway"]
        writer = csv.DictWriter(outputToFile,fieldnames = fieldname)
        writer.writeheader()
        for key in wikipathways.metabolitesWithPathwaysDictionary:
            pathways = wikipathways.metabolitesWithPathwaysDictionary[key]
            if(len(pathways) >0):
                writer.writerow({fieldname[0]:key,fieldname[1]:len(pathways)})
                '''
        idconvert.GeneConvert(wikipathways.geneInfoDictionary, "wiki")
        sql.checkForWithinDatabaseDuplicatesCompound(wikipathways.metaboliteIDDictionary, "wiki")
        sql.checkForWithinDatabaseDuplicatesGene(wikipathways.geneInfoDictionary, "wiki")
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
        '''
        print("Compound:") 
        stat.analyteOverlaps(sql.rampCompoundIdInWhichDatabases, sql.rampCompoundIDdictionary, "Compound")
        print("\n")
        print("Gene:") 
        stat.analyteOverlaps(sql.rampGeneIdInWhichDatabases, sql.rampGeneIDdictionary, "Gene")
        
        print("Compound:") 
        stat.analyteOverlaps(sql.rampCompoundIdInWhichDatabases, sql.rampCompoundIDdictionary, "Compound")
        print("\n")
        print("Gene:") 
        stat.analyteOverlaps(sql.rampGeneIdInWhichDatabases, sql.rampGeneIDdictionary, "Gene")
        '''




from reactomeData import reactomeData
from writeToSQL import writeToSQL
import time

import unittest

class TestReactomeMain(unittest.TestCase):

    def testMain(self):
        
        
        
        sql = writeToSQL()
        reactome = reactomeData()
        reactome.getEverything(True)
        reactome.getDatabaseFiles()
        print("Getting genes...")
        reactome.getGenes()
        print("Getting metabolites...")
        reactome.getMetabolites()
        
        print("Getting common names...")
        reactome.getCommonNameForChebi()
        
        
        print("Getting common names for genes ...")
        
        reactome.getGenes()
        reactome.downloadCommonNameFromUniprot()
        reactome.getCommonNameFromUniprot()
        reactome.write_myself_files('reactome')
        reactomecompoundnum = sql.createRampCompoundID(reactome.metaboliteIDDictionary, "reactome", 0)
        reactomegenenum = sql.createRampGeneID(reactome.geneInfoDictionary, "reactome", 0)
        
        reactomenumbers = sql.write(
                 reactome.metaboliteCommonName,
                 reactome.pathwayDictionary, 
                 reactome.pathwayCategory,
                 reactome.metabolitesWithPathwaysDictionary,
                 reactome.metabolitesWithSynonymsDictionary,
                 reactome.metaboliteIDDictionary,
                 reactome.pathwaysWithGenesDictionary,
                 reactome.metabolitesLinkedToGenes,
                 reactome.geneInfoDictionary,
                 reactome.biofluidLocation,
                 reactome.biofluid,
                 reactome.cellularLocation,
                 reactome.cellular,
                 reactome.pathwayOntology,
                 reactome.exoEndoDictionary,
                 reactome.exoEndo,
                 reactome.tissueLocation,
                 reactome.tissue,
                 "reactome",
                 0,0)
        print("Pathways number is " + str(len(reactome.pathwayDictionary)))
        print("metabolites number is " + str(len(reactome.metaboliteIDDictionary)))
        print('genes number is '+ str(len(reactome.geneInfoDictionary)))
        
if __name__ == "__main__":
    unittest.main()  
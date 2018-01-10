
from reactomeData import reactomeData
from writeToSQL import writeToSQL
import time

import unittest

class TestReactomeMain(unittest.TestCase):

    def testMain(self):
        
        
        
        sql = writeToSQL()
        reactome = reactomeData()
        #reactome.getDatabaseFiles()
        print("Getting genes...")
        reactome.getGenes()
        print("Getting metabolites...")
        reactome.getMetabolites()
        '''
        file = open("../misc/output/reactomeMetabolitesId.txt","wb")
        for key in reactome.metaboliteIDDictionary:
            file.write(key.encode("utf-8") +b"\n")
        file.close()
        '''
        print("Getting common names...")
        reactome.getCommonNameForChebi()
        file = open("../misc/output/reactomeGenesId.txt","wb")
        for key in reactome.geneInfoDictionary:
            file.write(key.encode("utf-8") +b"\n")
        file.close()
        print("Update gene Uniprot files")
        ids = reactome.getCommonNameForGenes1(hmdbdict=None,
                                             keggdict=None,
                                             )
        reactome.downloadCommonNameFromUniprot(ids)
        print("Getting common names for genes ...")
        reactome.getCommonNameForGenes2()
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
        
        
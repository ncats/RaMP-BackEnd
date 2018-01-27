from hmdbData import hmdbData


from MetabolomicsData import MetabolomicsData
import unittest
from numpy import var

class MetabolomicsDataWritenToFileTestCase(unittest.TestCase):
    def testMain(self):
        hmdb = hmdbData()
        tree = hmdb.getMetaboliteOtherIDs()
        
        hmdb.getPathwaysandSynonyms(tree)
        hmdb.getGenes(tree)
        hmdb.getBiofluidCellularLocationDisease(tree)
        hmdb.getPathwaysLinkedToGene()
        
        
        hmdb.write_myself_files(database='hmdb')
        
        
        
        
        
        
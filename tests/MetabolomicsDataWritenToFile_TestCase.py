from hmdbData import hmdbData

import unittest

class MetabolomicsDataWritenToFileTestCase(unittest.TestCase):
    def testMain(self):
        hmdb = hmdbData()
        tree = hmdb.getMetaboliteOtherIDs()
        
        hmdb.getPathwaysandSynonyms(tree)
        hmdb.getGenes(tree)
        hmdb.getPathwaysLinkedToGene()
        
        
        hmdb.write_myself_files(database='hmdb')
        
        
        
        
        
        
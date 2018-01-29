from update.RaMPDatabase import RaMPDatabase
from update.RaMPFixer import RaMPFixer
from hmdbData import hmdbData
import unittest
import timeit
class RaMPDatabaseTest(unittest.TestCase):
    def testMain(self):
        ramp = RaMPDatabase()
        
        fixer = RaMPFixer()
    
        hmdb=hmdbData()
        hmdb.getPathwaysandSynonyms(dir = 'sweat_metabolites.xml')
        print(hmdb.pathwayCategory)
        #print(df.loc[13,"rampId"])
        #print(type(df.loc[13,"rampId"]))

if __name__ == "__main__":
    unittest.main()
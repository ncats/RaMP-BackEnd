from update.RaMPDatabase import RaMPDatabase
from update.RaMPFixer import RaMPFixer
import unittest
import timeit
class RaMPDatabaseTest(unittest.TestCase):
    def testMain(self):
        ramp = RaMPDatabase()
        
        fixer = RaMPFixer()
        '''
        fixer.drop_database("mathelabramp2")
        fixer.create_new_db()
        fixer.create_tbs()
        '''
        df = fixer.remove_wrong_items_from_tb("source")
        
        print(df)
        print(df[df["commonName"] == "NA"])
        #print(df.loc[13,"rampId"])
        #print(type(df.loc[13,"rampId"]))

if __name__ == "__main__":
    unittest.main()
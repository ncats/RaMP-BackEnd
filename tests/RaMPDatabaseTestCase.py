from update.RaMPDatabase import RaMPDatabase
from update.RaMPFixer import RaMPFixer
import unittest
import timeit
class RaMPDatabaseTest(unittest.TestCase):
    def testMain(self):
        ramp = RaMPDatabase()
        df = ramp.connectToRaMP()
        fixer = RaMPFixer()
        fixer.hello_world()
        print(df)
        print(df.loc[13,"commonName"] == "")
        print(type(df.loc[13,"commonName"]))
        print(len(df))
        
        

if __name__ == "__main__":
    unittest.main()
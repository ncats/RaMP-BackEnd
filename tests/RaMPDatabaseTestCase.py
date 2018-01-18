from update.RaMPDatabase import RaMPDatabase
from update.RaMPFixer import RaMPFixer
import unittest
import timeit
class RaMPDatabaseTest(unittest.TestCase):
    def testMain(self):
        ramp = RaMPDatabase()
        
        fixer = RaMPFixer()
        fixer.create_new_db()
        
        
        

if __name__ == "__main__":
    unittest.main()
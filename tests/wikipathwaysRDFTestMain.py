import unittest
from wikipathwayRDF import WikipathwaysRDF



class TestWikipathwaysMain(unittest.TestCase):

    def testMain(self):
        wp = WikipathwaysRDF()
        #wp.getDatabaseFile()
        #wp.displayRDFfile(3)
        wp.getMetabolitesID()
        wp.write_myself_files('wikipathwayRDF')
        
        
        

if __name__ == '__main__':
    unittest.main()
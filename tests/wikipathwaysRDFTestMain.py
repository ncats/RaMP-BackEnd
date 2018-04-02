import unittest
from wikipathwayRDF import WikipathwaysRDF
from writeToSQL import writeToSQL
import pickle as pk

class TestWikipathwaysMain(unittest.TestCase):
    def testMain(self):
        wp = WikipathwaysRDF()
        sql = writeToSQL()
        #wp.getDatabaseFile()
        #wp.displayRDFfile(3)
        wp.getIDMapingWithPathways()
        wp.write_myself_files('wikipathwayRDF')
        sql.createRampCompoundID(wp.metaboliteIDDictionary, 'wiki', 0)
        sql.createRampGeneID(wp.geneInfoDictionary, 'wiki', 0)
        print('Unique metabolites: {} ||| Unique genes: {}'.format(len(sql.rampCompoundIDdictionary),
                                                                   len(sql.rampGeneIDdictionary)))
        wp_pk = open('../misc/output/wikipathwayRdfPk.pkl','wb')
        pk.dump(wp,wp_pk,pk.HIGHEST_PROTOCOL)
        wp_pk.close()
        del wp
        

if __name__ == '__main__':
    unittest.main()
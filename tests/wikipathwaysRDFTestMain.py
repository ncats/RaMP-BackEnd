import unittest
from wikipathwayRDF import WikipathwaysRDF
from writeToSQL import writeToSQL
import pickle as pk

class TestWikipathwaysMain(unittest.TestCase):
    def testMain(self):
        wp = WikipathwaysRDF()
        sql = writeToSQL()
        wp.getEverything(True)
        '''
        wp.getDatabaseFile()
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
        pkf = open('../misc/output/wikipathwayRdfPk.pkl','rb')
        wp = pk.load(pkf)
        pkf.close()
        rp = RampUpdater(wp)
        rp.checkNewAnalyteEntry('compound')
        rp.checkNewAnalyteEntry('gene')
        rp.checkNewPathwayEntry()
        pkf2 = open('../misc/output/updateObject328.pkl','wb')
        pk.dump(rp,pkf2,pk.HIGHEST_PROTOCOL)
        pkf2.close()
        pkf = open('../misc/output/updateObject328.pkl','rb')
        rp = pk.load(pkf)
        print('Unique metabolites: {} \nUnique genes: {}'.format(len(set(rp.newRampCompound.values())),
                                                                  len(set(rp.newRampGene.values()))))
        
        time.sleep(1)
        pkf.close()
        rp2 = RampUpdater(rp)
        rp2.newRampCompound = rp.newRampCompound
        rp2.newRampPathway = rp.newRampPathway
        rp2.newRampGene = rp.newRampGene
        rp2.writeToRamp('wiki')
        '''

if __name__ == '__main__':
    unittest.main()
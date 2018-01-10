from hmdbData import hmdbData
from KeggData import KeggData
from reactomeData import reactomeData
from wikipathwaysData import wikipathwaysData
from writeToSQL import writeToSQL
from getStatistics import getStatistics
from IDconversion import IDconversion
import random
import time
import unittest

class testForFindingDB(unittest.TestCase):
    def testMain(self):
        # hmdb = hmdbData()
        # kegg = KeggData()
        wiki = wikipathwaysData()
        #react = reactomeData()
        print("Get metabolites ...")
        # hmdb.getMetaboliteOtherIDs()
        #kegg.getPathways()
        #kegg.getMetabolites()
        #kegg.getSynonymsAndCHEBI()
        wiki.getEverything()
        wiki.getCommonNameForChebi()
        #react.getMetabolites()
        #react.getCommonNameForChebi()
        chebilist = []
        hmdbHasChebiId = 0
        sizeHmdb = len(wiki.metaboliteIDDictionary)
        rkey = random.choice(list(wiki.metaboliteIDDictionary.keys()))
        idmapping = wiki.metaboliteIDDictionary[rkey]
        idCountDict = dict()
        for key in idmapping:
            idCountDict[key] = 0
        print(idCountDict)
        for key in wiki.metaboliteIDDictionary:
            idMap = wiki.metaboliteIDDictionary[key]
            print(idMap)
           
            for key in idMap:
                if idMap[key] is not 'NA':
                    idCountDict[key] +=1
        
        print(len(wiki.metaboliteIDDictionary))    
        print(idCountDict)
        print(wiki.setOfType)
        print(wiki.idSetFromGeneProducts.keys())
        print(len(wiki.idSetFromGeneProducts['Entrez Gene']))
        print(wiki.idSetFromProtein.keys())
        print(len(wiki.idSetFromProtein['Entrez Gene']))
        print(len(wiki.idSetFromGeneProducts['Entrez Gene'] & wiki.idSetFromProtein['Entrez Gene']))
        
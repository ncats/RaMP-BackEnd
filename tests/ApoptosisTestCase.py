from reactomeData import reactomeData
from KeggData import KeggData
from wikipathwaysData import wikipathwaysData
from getStatistics import getStatistics
from writeToSQL import writeToSQL

import unittest

class apoptosisTest(unittest.TestCase):
     
    def testMain(self):
        
        wiki = wikipathwaysData()
        reactome = reactomeData()
        kegg = KeggData()
        sql = writeToSQL()
        stat = getStatistics()
        
        wiki.pathwaysWithGenesDictionary["WP254"] = ["geneA", "geneB"]
        reactome.pathwaysWithGenesDictionary["R-HSA-109581"] = ["geneA", "geneD"]
        kegg.pathwaysWithGenesDictionary["04210"] = ["geneC", "geneB"]
        
        sql.rampGeneIDdictionary["geneA"] = "RAMP00001"
        sql.rampGeneIDdictionary["geneB"] = "RAMP00002"
        sql.rampGeneIDdictionary["geneC"] = "RAMP00003"
        sql.rampGeneIDdictionary["geneD"] = "RAMP00001"
        
        stat.Apoptosis(sql.rampGeneIDdictionary, wiki.pathwaysWithGenesDictionary, kegg.pathwaysWithGenesDictionary, reactome.pathwaysWithGenesDictionary)
        
        
        
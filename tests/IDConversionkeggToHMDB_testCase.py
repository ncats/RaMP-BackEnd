from hmdbData import hmdbData
from KeggData import KeggData
from IDconversion import IDconversion
import unittest
from writeToSQL import writeToSQL

class TestIDconvsion(unittest.TestCase):
     
     def testKeggToHMDB(self):
         
         hmdb = hmdbData()
         kegg = KeggData()
         sql = writeToSQL()
         
         #metabolite mapping for hmdb
         hmdb.metaboliteIDDictionary["HMDB00001"] = {"chebi_id": "34131", 
                           "drugbank_id": "NA", 
                           "drugbank_metabolite_id": "NA", 
                           "phenol_explorer_compound_id": "NA", 
                           "phenol_explorer_metabolite_id": "NA", 
                           "foodb_id": "FDB012119", 
                           "knapsack_id": "NA", 
                           "chemspider_id": "83153",
                           "kegg_id": "NA",
                           "biocyc_id": "CPD-1823",
                           "bigg_id": "NA",
                           "wikipidia": "NA",
                           "nugowiki": "NA",
                           "metagene": "NA",
                           "metlin_id": "3741",
                           "pubchem_compound_id": "92105",
                           "het_id": "HIC",
                           "hmdb_id": ["HMDB00001"],
                           "CAS": "NA"}
         
         
         #metabolite mapping for kegg        
         kegg.metaboliteIDDictionary["C14814"] = {"chebi_id": ["34131"], 
                    "drugbank_id": "NA", 
                    "drugbank_metabolite_id": "NA", 
                    "phenol_explorer_compound_id": "NA", 
                    "phenol_explorer_metabolite_id": "NA", 
                    "foodb_id": "NA", 
                    "knapsack_id": "NA", 
                    "chemspider_id": "NA",
                    "kegg_id": "NA",
                    "biocyc_id": "NA",
                    "bigg_id": "NA",
                    "wikipidia": "NA",
                    "nugowiki": "NA",
                    "metagene": "NA",
                    "metlin_id": "NA",
                    "pubchem_compound_id": "NA",
                    "het_id": "NA",
                    "hmdb_id": "NA",
                    "CAS": "NA"}
         
         IDconversion.MetaboliteChebiToHMDB(self, kegg.metaboliteIDDictionary, hmdb.metaboliteIDDictionary, "kegg")
         
        
            #Pathway names
         kegg.pathwayDictionary["00010"] = "Glycolysis / Gluconeogenesis"
         kegg.pathwayDictionary["00020"] = "Citrate cycle (TCA cycle)"
         kegg.pathwayDictionary["00520"] = "Fake Pathway Name One"
         kegg.pathwayDictionary["00524"] = "Fake Pathway Name Two"
         kegg.pathwayDictionary["00540"] = "Fake Pathway Name Three"
         kegg.pathwayDictionary["00550"] = "Fake Pathway Name Four"
         kegg.pathwayDictionary["00030"] = "Fake Pathway Name Five"
         kegg.pathwayDictionary["00040"] = "Fake Pathway Name Six"
         kegg.pathwayDictionary["00053"] = "Fake Pathway Name Seven"
         kegg.pathwayDictionary["00250"] = "Fake Pathway Name Eight"
         kegg.pathwayDictionary["00260"] = "Fake Pathway Name Nine"
            
         #Pathway categories
         kegg.pathwayCategory["00010"] = "Metabolism"
         kegg.pathwayCategory["00020"] = "Human Diseases"
         kegg.pathwayCategory["00520"] = "Cellular Processes"
         kegg.pathwayCategory["00524"] = "Human Diseases"
         kegg.pathwayCategory["00540"] = "Human Diseases"
         kegg.pathwayCategory["00550"] = "Metabolism"
         kegg.pathwayCategory["00030"] = "Human Diseases"
         kegg.pathwayCategory["00040"] = "Cellular Processes"
         kegg.pathwayCategory["00053"] = "Metabolism"
         kegg.pathwayCategory["00250"] = "Human Diseases"
         kegg.pathwayCategory["00260"] = "Cellular Processes"
            
         #metabolites linked with pathways 
         kegg.metabolitesWithPathwaysDictionary["C00043"] = ["00520", "00524", "00540", "00550"]
         kegg.metabolitesWithPathwaysDictionary["C00022"] = ["00010", "00020", "00030", "00040", "00053", "00250", "00260"]
         kegg.metabolitesWithPathwaysDictionary["C14814"] = ["00020"]
         kegg.metabolitesWithPathwaysDictionary["C14865"] = ["00540"]
            
            
         #metabolites linkes with synonyms
         kegg.metabolitesWithSynonymsDictionary["C00043"] = ["UDP-N-acetyl-alpha-D-glucosamine", "UDP-N-acetyl-D-glucosamine", "UDP-N-acetylglucosamine"]
         kegg.metabolitesWithSynonymsDictionary["C00022"] = ["Pyruvate", "Pyruvic acid", "2-Oxopropanoate", "2-Oxopropanoic acid", "Pyroracemic acid"]
         kegg.metabolitesWithSynonymsDictionary["C14814"] = ["MetaboliteSynonym1"]
         kegg.metabolitesWithSynonymsDictionary["C14865"] = ["MetaboliteSynonym2"]
         
     
        
         #pathway to gene id
         kegg.pathwaysWithGenesDictionary["00010"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00020"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00520"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00524"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00540"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00550"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00030"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00040"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00053"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00250"] = ["geneA", "geneB"]
         kegg.pathwaysWithGenesDictionary["00260"] = ["geneA", "geneB"]

        
         #gene to geneinfo
         kegg.geneInfoDictionary["geneA"] = {'common_name': 'Apple',
                                            'kegg': 'NA',
                                            'Ensembl': 'ENSG00000127481', 
                                            'HGNC': '30313', 
                                            'HPRD': 'NA', 
                                            'NCBI-GeneID': '23352', 
                                            'NCBI-ProteinID': 'NP_065816', 
                                            'OMIM': '609890', 
                                            'UniProt': 'Q5T4S7', 
                                            'Vega': 'OTTHUMG00000002498', 
                                            'miRBase': 'NA', 
                                            'HMDB_protien_accession': 'NA',
                                            'Entrez': 'NA',
                                            'Enzyme Nomenclature': 'NA'}
         kegg.geneInfoDictionary["geneB"] = {'common_name': 'Banana',
                                            'kegg': 'NA',
                                            'Ensembl': 'ENSG00000100320', 
                                            'HGNC': '9906', 
                                            'HPRD': 'NA', 
                                            'NCBI-GeneID': '23543', 
                                            'NCBI-ProteinID': 'NP_065816', 
                                            'OMIM': '612149', 
                                            'UniProt': 'O43251', 
                                            'Vega': 'OTTHUMG00000150585', 
                                            'miRBase': 'NA', 
                                            'HMDB_protien_accession': 'NA',
                                            'Entrez': 'NA',
                                            'Enzyme Nomenclature': 'NA'}
         
         sql.write(kegg.pathwayDictionary, 
                 kegg.pathwayCategory,
                 kegg.metabolitesWithPathwaysDictionary,
                 kegg.metabolitesWithSynonymsDictionary,
                 kegg.metaboliteIDDictionary,
                 kegg.pathwaysWithGenesDictionary,
                 kegg.geneInfoDictionary,
                 kegg.biofluidLocation,
                 kegg.biofluid,
                 kegg.cellularLocation,
                 kegg.cellular,
                 kegg.pathwayOntology,
                 "kegg",
                 0,0,0,0)

         
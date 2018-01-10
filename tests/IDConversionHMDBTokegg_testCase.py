from hmdbData import hmdbData
from KeggData import KeggData
from IDconversion import IDconversion
import unittest
from writeToSQL import writeToSQL
from getStatistics import getStatistics


class TestIDconvsion(unittest.TestCase):
     '''
     The purpose of this test case is to check *METABOLITE* ID conversion from kegg to chebi and then from chebi to HMDB. 
     The section that is most important for this is noted in the test case. It is the portion holding the metabolite ID
     information. The rest of the code is essential for testing since it is required to write information to sql files. 
     
     '''
     
     
     def testKeggToHMDB(self):
         
         ###############################################################################################
         #IMPORTANT PART START
         
         hmdb = hmdbData()
         kegg = KeggData()
         sql = writeToSQL()
         idconvert = IDconversion()
         stat = getStatistics()
         
         #metabolite mapping for hmdb
         hmdb.metaboliteIDDictionary["HMDB00001"] = {"chebi_id": "NA", 
                           "drugbank_id": "NA", 
                           "drugbank_metabolite_id": "NA", 
                           "phenol_explorer_compound_id": "NA", 
                           "phenol_explorer_metabolite_id": "NA", 
                           "foodb_id": "FDB012119", 
                           "knapsack_id": "NA", 
                           "chemspider_id": "83153",
                           "kegg_id": "C14814",
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
                    "kegg_id": "C14814",
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
         
         idconvert.MetaboliteKeggIDToChebi(kegg.metaboliteIDDictionary, hmdb.metaboliteIDDictionary, "hmdb")
         idconvert.MetaboliteChebiToHMDB(kegg.metaboliteIDDictionary, hmdb.metaboliteIDDictionary, "kegg")
         
         #IMPORTANT PART END
         ######################################################################################################
         
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
         kegg.metabolitesWithPathwaysDictionary["C14814"] = ["00020"]

            
            
         #metabolites linkes with synonyms
         kegg.metabolitesWithSynonymsDictionary["C14814"] = ["MetaboliteSynonym1"]
     
        
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
         
         hmdb.metabolitesWithSynonymsDictionary["HMDB00001"] = ["1 Methylhistidine", "1-Methyl-L-histidine", "Pi-methylhistidine"]
         
         hmdb.metabolitesWithPathwaysDictionary["HMDB00001"] = ["SMP00716", "SMP00006"]

         
         hmdb.pathwayDictionary["SMP00716"] = "Thyroid hormone synthesis"
         hmdb.pathwayDictionary["SMP00006"] = "Tyrosine Metabolism"
         hmdb.pathwayDictionary["SMP00001"] = "Pathway1"
         hmdb.pathwayDictionary["SMP00002"] = "Pathway2"
         hmdb.pathwayDictionary["SMP00816"] = "Pathway3"
         
         hmdb.pathwayCategory["SMP00716"] = "NA"
         hmdb.pathwayCategory["SMP00006"] = "NA"
         hmdb.pathwayCategory["SMP00001"] = "NA"
         hmdb.pathwayCategory["SMP00002"] = "NA"
         hmdb.pathwayCategory["SMP00816"] = "NA"
         
         hmdb.pathwaysWithGenesDictionary["SMP00716"] = ["Q96KN2", "uniprot1"]
         hmdb.pathwaysWithGenesDictionary["SMP00006"] = ["Q96KN2", "uniprot1"]
         hmdb.pathwaysWithGenesDictionary["SMP00001"] = ["Q96KN2", "uniprot1"]
         hmdb.pathwaysWithGenesDictionary["SMP00002"] = ["Q96KN2", "uniprot1"]
         hmdb.pathwaysWithGenesDictionary["SMP00816"] = ["Q96KN2", "uniprot1"]
         
         hmdb.geneInfoDictionary["Q96KN2"] = { 'common_name': 'CNDP1',
                                           'kegg': 'NA',
                                           'Ensembl': 'NA', 
                                           'HGNC': 'NA', 
                                           'HPRD': 'NA', 
                                           'NCBI-GeneID': 'NA', 
                                           'NCBI-ProteinID': 'NA', 
                                           'OMIM': 'NA', 
                                           'UniProt': 'Q96KN2', 
                                           'Vega': 'NA', 
                                           'miRBase': 'NA', 
                                           'HMDB_protien_accession': 'HMDBP00473',
                                           'Entrez': 'NA',
                                           'Enzyme Nomenclature': 'NA'}
         
         hmdb.geneInfoDictionary["uniprot1"] = { 'common_name': 'genename1',
                                           'kegg': 'NA',
                                           'Ensembl': 'NA', 
                                           'HGNC': 'NA', 
                                           'HPRD': 'NA', 
                                           'NCBI-GeneID': 'NA', 
                                           'NCBI-ProteinID': 'NA', 
                                           'OMIM': 'NA', 
                                           'UniProt': 'uniprot1', 
                                           'Vega': 'NA', 
                                           'miRBase': 'NA', 
                                           'HMDB_protien_accession': 'HMDBP00321',
                                           'Entrez': 'NA',
                                           'Enzyme Nomenclature': 'NA'}
         
         
         hmdb.biofluidLocation["HMDB00001"] = ["Blood", "Cerebrospinal Fluid (CSF)", "Feces", "Saliva", "Urine"]
         
         hmdb.exoEndoDictionary["HMDB00001"] = ["Food"]
         
         
         hmdb.biofluid["Blood"] = "placeholder"
         hmdb.biofluid["Cerebrospinal Fluid (CSF)"] = "placeholder"
         hmdb.biofluid["Feces"] = "placeholder"
         hmdb.biofluid["Saliva"] = "placeholder"
         hmdb.biofluid["Urine"] = "placeholder"
         
         hmdb.cellularLocation["HMDB00001"] = ["Cytoplasm", "Location1"]

         
         hmdb.cellular["Cytoplasm"] = "placeholder"
         hmdb.cellular["Location1"] = "placeholder"
         hmdb.cellular["Location2"] = "placeholder"
         hmdb.cellular["Location3"] = "placeholder"
         
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
                 kegg.exoEndoDictionary,
                 "kegg",
                 0,0)
         
         sql.write(hmdb.pathwayDictionary, 
                 hmdb.pathwayCategory,
                 hmdb.metabolitesWithPathwaysDictionary,
                 hmdb.metabolitesWithSynonymsDictionary,
                 hmdb.metaboliteIDDictionary,
                 hmdb.pathwaysWithGenesDictionary,
                 hmdb.geneInfoDictionary,
                 hmdb.biofluidLocation,
                 hmdb.biofluid,
                 hmdb.cellularLocation,
                 hmdb.cellular,
                 hmdb.pathwayOntology,
                 hmdb.exoEndoDictionary,
                 "hmdb",
                 0,0)
         
         stat.analyteOverlaps(sql.rampCompoundIdInWhichDatabases)
         
         stat.analyteOverlaps(sql.rampGeneIdInWhichDatabases)
         


         
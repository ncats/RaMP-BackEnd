from hmdbData import hmdbData
from reactomeData import reactomeData
from IDconversion import IDconversion
import unittest
from getStatistics import getStatistics
from writeToSQL import writeToSQL


class TestIDconvsion(unittest.TestCase):
     
     def testreactomeToHMDB(self):  
         hmdb = hmdbData()
         reactome = reactomeData()
         sql = writeToSQL()
         idconvert = IDconversion()
         stat = getStatistics()
         
         #metabolite mapping for hmdb
         hmdb.metaboliteIDDictionary["HMDB00001"] = {"chebi_id": ["C14814"], 
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
         reactome.metaboliteIDDictionary["C14814"] = {"chebi_id": ["C14814"], 
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


         hmdb.geneInfoDictionary["Q96KN2"] = { 'common_name': 'genename1',
                                           'kegg': 'NA',
                                           'Ensembl': 'NA', 
                                           'HGNC': 'NA', 
                                           'HPRD': 'NA', 
                                           'NCBI-GeneID': 'NA', 
                                           'NCBI-ProteinID': 'NA', 
                                           'OMIM': 'NA', 
                                           'UniProt': ['Q96KN2'], 
                                           'Vega': 'NA', 
                                           'miRBase': 'NA', 
                                           'HMDB_protien_accession': 'HMDBP00321',
                                           'Entrez': 'NA',
                                           'Enzyme Nomenclature': 'NA'}
         reactome.geneInfoDictionary["Q96KN2"] = {'common_name': 'NA',
                                            'kegg': 'NA',
                                            'Ensembl': 'NA', 
                                            'HGNC': 'NA', 
                                            'HPRD': 'NA', 
                                            'NCBI-GeneID': 'NA', 
                                            'NCBI-ProteinID': 'NA', 
                                            'OMIM': 'NA', 
                                            'UniProt': ['Q96KN2'], 
                                            'Vega': 'NA', 
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
         
         hmdb.pathwaysWithGenesDictionary["SMP00716"] = ["Q96KN2"]
         hmdb.pathwaysWithGenesDictionary["SMP00006"] = ["Q96KN2"]
         hmdb.pathwaysWithGenesDictionary["SMP00001"] = ["Q96KN2"]
         hmdb.pathwaysWithGenesDictionary["SMP00002"] = ["Q96KN2"]
         hmdb.pathwaysWithGenesDictionary["SMP00816"] = ["Q96KN2"]
         
         hmdb.geneInfoDictionary["Q96KN2"] = { 'common_name': 'CNDP1',
                                           'kegg': 'NA',
                                           'Ensembl': 'NA', 
                                           'HGNC': 'NA', 
                                           'HPRD': 'NA', 
                                           'NCBI-GeneID': 'NA', 
                                           'NCBI-ProteinID': 'NA', 
                                           'OMIM': 'NA', 
                                           'UniProt': ['Q96KN2'], 
                                           'Vega': 'NA', 
                                           'miRBase': 'NA', 
                                           'HMDB_protien_accession': 'HMDBP00473',
                                           'Entrez': 'NA',
                                           'Enzyme Nomenclature': 'NA'}
         
         hmdb.biofluidLocation["HMDB00001"] = ["Blood", "Cerebrospinal Fluid (CSF)", "Feces", "Saliva", "Urine"]
         
         
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
         
         hmdb.exoEndoDictionary["HMDB00001"] = "Food"
         

         reactome.pathwayDictionary["R-HSA-210745"] = "Citrate cycle (TCA cycle)"
            
         #Pathway categories
         reactome.pathwayCategory["R-HSA-210745"] = "Human Diseases"
        
         #metabolites linked with pathways 
         reactome.metabolitesWithPathwaysDictionary["C14814"] = ["R-HSA-210745"]
         #metabolites linkes with synonyms
         reactome.metabolitesWithSynonymsDictionary["C14814"] = ["MetaboliteSynonym1"]
     
         #pathway to gene id
         reactome.pathwaysWithGenesDictionary["R-HSA-210745"] = ["Q96KN2"]
         
         idconvert.MetaboliteChebiToHMDB(reactome.metaboliteIDDictionary, hmdb.metaboliteIDDictionary, "reactome")
         idconvert.GeneUniprotToHMDBP(reactome.geneInfoDictionary, hmdb.geneInfoDictionary, "hmdb")
         
                  
         hmdbcompoundnum = sql.createRampCompoundID(hmdb.metaboliteIDDictionary, "hmdb", 0)
         reactomecompoundnum = sql.createRampCompoundID(reactome.metaboliteIDDictionary, "reactome", hmdbcompoundnum)
        
         hmdbgenenum = sql.createRampGeneID(hmdb.geneInfoDictionary, "hmdb", 0)
         reactomegenenum = sql.createRampGeneID(reactome.geneInfoDictionary, "reactome", hmdbgenenum)

         
         sql.write(reactome.pathwayDictionary, 
                 reactome.pathwayCategory,
                 reactome.metabolitesWithPathwaysDictionary,
                 reactome.metabolitesWithSynonymsDictionary,
                 reactome.metaboliteIDDictionary,
                 reactome.pathwaysWithGenesDictionary,
                 reactome.geneInfoDictionary,
                 reactome.biofluidLocation,
                 reactome.biofluid,
                 reactome.cellularLocation,
                 reactome.cellular,
                 reactome.pathwayOntology,
                 reactome.exoEndoDictionary,
                 "reactome",
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
         
       
         print("Compound:") 
         stat.analyteOverlaps(sql.rampCompoundIdInWhichDatabases, sql.rampCompoundIDdictionary, "Compound")
         print("\n")
         print("Gene:") 
         stat.analyteOverlaps(sql.rampGeneIdInWhichDatabases, sql.rampGeneIDdictionary, "Gene")
         
         stat.databaseContent(hmdb.pathwayDictionary, 
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
                 "hmdb")
         

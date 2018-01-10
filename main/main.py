import sys
sys.path.append("../")
from hmdbData import hmdbData
from KeggData import KeggData
from reactomeData import reactomeData
from wikipathwaysData import wikipathwaysData
from IDconversion import IDconversion
from getStatistics import getStatistics
from writeToSQL import writeToSQL
import os
os.chdir("../RaMP/main/")

class Main():

    def runEverything(self, getDatabaseFiles = False):
        sql = writeToSQL()
        idconvert = IDconversion()
        stat = getStatistics()
        hmdb = hmdbData()
        wikipathways = wikipathwaysData()
        reactome = reactomeData()
        kegg = KeggData()
        print(os.getcwd())
        #pulls needed files from each database if true. Otherwise, assumes files already present. Default false.
        if getDatabaseFiles:
            kegg.getDatabaseFiles()
            wikipathways.getDatabaseFiles()
            reactome.getDatabaseFiles()
            hmdb.getDatabaseFiles()
            

        print("Getting HMDB Metabolites...")
        hmdb.getMetaboliteOtherIDs()
        print("Getting HMDB pathways and synonyms...")
        hmdb.getPathwaysandSynonyms()
        print("Getting HMDB genes...")
        hmdb.getGenes()
        print("Getting HMDB biofluid and cellular locations...")
        hmdb.getBiofluidCellularLocationDisease()
        print("Getting HMDB pathways links to genes ...")
        hmdb.getPathwaysLinkedToGene()
        
        print("Getting wikipathways...")
        wikipathways.getEverything()
        wikipathways.getCommonNameForChebi()
        
        

        
        print("Getting reactome genes...")
        reactome.getGenes()
        print("Getting reactome metabolites...")
        reactome.getMetabolites()
        reactome.getCommonNameForChebi()
        reactome.getCommonNameForGenes()
        print("Getting kegg pathways...")
        kegg.getPathways()
        print("Getting kegg genes and metabolites...")
        kegg.getMetabolites()
        kegg.getSynonymsAndCHEBI()
        kegg.getGenes()
        kegg.getGeneInfo()
        
        print("Converting gene ids...")
        #Here are the identifiers that are present for each gene:
        #kegg: keggid (mainID), 'Ensembl', 'HGNC', 'HPRD', 'NCBI-GeneID', 'NCBI-ProteinID', 'OMIM', 'UniProt', 'Vega', 'miRBase'
        #wikipathways: (no mainID), 'Entrez', 'Enzyme Nomenclature', 'Uniprot (Uniprot-TrEMBL)
        #hmdb: HMDB-protien-accession (mainID), 'Uniprot'
        #reactome:Uniprot (mainID)
        idconvert.GeneConvert(wikipathways.geneInfoDictionary, "wikipathways")
        idconvert.GeneConvert(hmdb.geneInfoDictionary, "hmdb")
        idconvert.GeneConvert(reactome.geneInfoDictionary, "reactome")
        idconvert.GeneConvert(kegg.geneInfoDictionary, "kegg")
        
        idconvert.GeneUniprotToHMDBP(wikipathways.geneInfoDictionary, hmdb.geneInfoDictionary, "wikipathways")
        idconvert.GeneUniprotToHMDBP(reactome.geneInfoDictionary, hmdb.geneInfoDictionary, "reactome")
        idconvert.GeneUniprotToHMDBP(kegg.geneInfoDictionary, hmdb.geneInfoDictionary, "kegg")
        
        print("Converting metabolite ids...")
        idconvert.MetaboliteKeggIDToChebi(kegg.metaboliteIDDictionary, hmdb.metaboliteIDDictionary, "hmdb") 
        idconvert.MetaboliteChebiToHMDB(wikipathways.metaboliteIDDictionary, hmdb.metaboliteIDDictionary, "wikipathways")
        idconvert.MetaboliteChebiToHMDB(reactome.metaboliteIDDictionary, hmdb.metaboliteIDDictionary, "reactome")
        idconvert.MetaboliteChebiToHMDB(kegg.metaboliteIDDictionary, hmdb.metaboliteIDDictionary, "kegg")
        
        #check for dups
        print("Wikipathways compounds...")
        sql.checkForWithinDatabaseDuplicatesCompound(wikipathways.metaboliteIDDictionary, "wikipathways")
        print("Wikipathways genes...")
        sql.checkForWithinDatabaseDuplicatesGene(wikipathways.geneInfoDictionary, "wikipathways")
        print("Kegg compounds...")
        sql.checkForWithinDatabaseDuplicatesCompound(kegg.metaboliteIDDictionary, "kegg")
        print("kegg genes...")
        sql.checkForWithinDatabaseDuplicatesGene(kegg.geneInfoDictionary, "kegg")
        print("reactome compounds...")
        sql.checkForWithinDatabaseDuplicatesCompound(reactome.metaboliteIDDictionary, "reactome")
        print("reactome genes...")
        sql.checkForWithinDatabaseDuplicatesGene(reactome.geneInfoDictionary, "reactome")
        print("hmdb compounds...")
        sql.checkForWithinDatabaseDuplicatesCompound(hmdb.metaboliteIDDictionary, "hmdb")
        print("hmdb genes...")
        sql.checkForWithinDatabaseDuplicatesGene(hmdb.geneInfoDictionary, "hmdb")
        
        print('Generate compound id')
        hmdbcompoundnum = sql.createRampCompoundID(hmdb.metaboliteIDDictionary, "hmdb", 0)
        wikicompoundnum = sql.createRampCompoundID(wikipathways.metaboliteIDDictionary, "wiki", hmdbcompoundnum)
        reactomecompoundnum = sql.createRampCompoundID(reactome.metaboliteIDDictionary, "reactome", wikicompoundnum)
        keggcompoundnum = sql.createRampCompoundID(kegg.metaboliteIDDictionary, "kegg", reactomecompoundnum)
        print('Generate gene id ...')
        hmdbgenenum = sql.createRampGeneID(hmdb.geneInfoDictionary, "hmdb", 0)
        wikigenenum = sql.createRampGeneID(wikipathways.geneInfoDictionary, "wiki", hmdbgenenum)
        reactomegenenum = sql.createRampGeneID(reactome.geneInfoDictionary, "reactome", wikigenenum)
        kegggenenum = sql.createRampGeneID(kegg.geneInfoDictionary, "kegg", reactomegenenum)
        print('Write to sql file...')
        hmdbnumbers = sql.write(
                hmdb.metaboliteCommonName,
                hmdb.pathwayDictionary, 
                 hmdb.pathwayCategory,
                 hmdb.metabolitesWithPathwaysDictionary,
                 hmdb.metabolitesWithSynonymsDictionary,
                 hmdb.metaboliteIDDictionary,
                 hmdb.pathwaysWithGenesDictionary,
                 hmdb.metabolitesLinkedToGenes,
                 hmdb.geneInfoDictionary,
                 hmdb.biofluidLocation,
                 hmdb.biofluid,
                 hmdb.cellularLocation,
                 hmdb.cellular,
                 hmdb.pathwayOntology,
                 hmdb.exoEndoDictionary,
                 hmdb.exoEndo,
                 hmdb.tissueLocation,
                 hmdb.tissue,
                 "hmdb",
                 0,0)
        
        wikipathwaysnumbers = sql.write(
                wikipathways.metaboliteCommonName,
                wikipathways.pathwayDictionary, 
                 wikipathways.pathwayCategory,
                 wikipathways.metabolitesWithPathwaysDictionary,
                 wikipathways.metabolitesWithSynonymsDictionary,
                 wikipathways.metaboliteIDDictionary,
                 wikipathways.pathwaysWithGenesDictionary,
                 wikipathways.metabolitesLinkedToGenes,
                 wikipathways.geneInfoDictionary,
                 wikipathways.biofluidLocation,
                 wikipathways.biofluid,
                 wikipathways.cellularLocation,
                 wikipathways.cellular,
                 wikipathways.pathwayOntology,
                 wikipathways.exoEndoDictionary,
                 wikipathways.exoEndo,
                 wikipathways.tissueLocation,
                 wikipathways.tissue,
                 "wiki",
                 hmdbnumbers[0],hmdbnumbers[1])
        
        reactomenumbers = sql.write(
                reactome.metaboliteCommonName,
                reactome.pathwayDictionary, 
                 reactome.pathwayCategory,
                 reactome.metabolitesWithPathwaysDictionary,
                 reactome.metabolitesWithSynonymsDictionary,
                 reactome.metaboliteIDDictionary,
                 reactome.pathwaysWithGenesDictionary,
                 reactome.metabolitesLinkedToGenes,
                 reactome.geneInfoDictionary,
                 reactome.biofluidLocation,
                 reactome.biofluid,
                 reactome.cellularLocation,
                 reactome.cellular,
                 reactome.pathwayOntology,
                 reactome.exoEndoDictionary,
                 reactome.exoEndo,
                 reactome.tissueLocation,
                 reactome.tissue,
                 "reactome",
                 wikipathwaysnumbers[0],wikipathwaysnumbers[1])
        
        keggnumbers = sql.write(
                kegg.metaboliteCommonName,
                kegg.pathwayDictionary, 
                 kegg.pathwayCategory,
                 kegg.metabolitesWithPathwaysDictionary,
                 kegg.metabolitesWithSynonymsDictionary,
                 kegg.metaboliteIDDictionary,
                 kegg.pathwaysWithGenesDictionary,
                 kegg.metabolitesLinkedToGenes,
                 kegg.geneInfoDictionary,
                 kegg.biofluidLocation,
                 kegg.biofluid,
                 kegg.cellularLocation,
                 kegg.cellular,
                 kegg.pathwayOntology,
                 kegg.exoEndoDictionary,
                 kegg.exoEndo,
                 kegg.tissueLocation,
                 kegg.tissue,
                 "kegg",
                 reactomenumbers[0],reactomenumbers[1])
        print("Done ... for importing database")
        
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
         
        stat.databaseContent(kegg.pathwayDictionary, 
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
                 "kegg")
          
        stat.databaseContent(reactome.pathwayDictionary, 
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
                 "reactome")
           
        stat.databaseContent(wikipathways.pathwayDictionary, 
                 wikipathways.pathwayCategory,
                 wikipathways.metabolitesWithPathwaysDictionary,
                 wikipathways.metabolitesWithSynonymsDictionary,
                 wikipathways.metaboliteIDDictionary,
                 wikipathways.pathwaysWithGenesDictionary,
                 wikipathways.geneInfoDictionary,
                 wikipathways.biofluidLocation,
                 wikipathways.biofluid,
                 wikipathways.cellularLocation,
                 wikipathways.cellular,
                 wikipathways.pathwayOntology,
                 wikipathways.exoEndoDictionary,
                 "wiki")
        
        stat.Apoptosis(sql.rampGeneIDdictionary,
                        wikipathways.pathwaysWithGenesDictionary, 
                        kegg.pathwaysWithGenesDictionary, 
                        reactome.pathwaysWithGenesDictionary)

        
        
    
        
        
main = Main()
main.runEverything()

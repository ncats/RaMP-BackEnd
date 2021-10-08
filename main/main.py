import sys
sys.path.append('../src')
from parse.wikipathwayRDF import WikipathwaysRDF
from parse.hmdbData import hmdbData
from parse.reactomeData import reactomeData
from parse.lipidmapsChemData import lipidmapsChemData
from parse.KeggData import KeggData
from util.EntityBuilder import EntityBuilder
from getStatistics import getStatistics
from writeToSQL import writeToSQL
import os
import time

class Main():

    def runEverything(self, getDatabaseFiles = True, resourceConfigFile):
        sql = writeToSQL()
        
        # build the ramp resource config
        resourceConf = RampConfig()
        resourceConf.loadConfig(resourceConfigFile)
        
        stat = getStatistics()
        hmdb = hmdbData()
        wikipathways = WikipathwaysRDF()
        reactome = reactomeData()
        kegg = KeggData()
        lipidmaps = lipidmapsChemData()
        
        # works based on your computer, setup working directory
        os.chdir('../main/')

        #kegg.getEverything(False)
        #print("KEGG Wonder")
        print("Getting hmdb...")
        hmdb.getEverything(True)
        print("Getting wiki...")
        wikipathways.getEverything(True)
        print("Getting reactome...")
        reactome.getEverything(True)
        
        # This parses and writes lipid maps
        # sql write will be handled by EntityBuilder
        print("Getting LipidMaps...")
        lipidmaps.getEverything(True)

        #Here are the identifiers that are present for each gene:
        #kegg: keggid (mainID), 'Ensembl', 'HGNC', 'HPRD', 'NCBI-GeneID', 'NCBI-ProteinID', 'OMIM', 'UniProt', 'Vega', 'miRBase'
        #wikipathways: (no mainID), 'Entrez', 'Enzyme Nomenclature', 'Uniprot (Uniprot-TrEMBL)
        #hmdb: HMDB-protien-accession (mainID), 'Uniprot'
        #reactome:Uniprot (mainID)
        
        print('Generate compound id')
        hmdbcompoundnum = sql.createRampCompoundID(hmdb.metaboliteIDDictionary, "hmdb", 0)
        print("hmdbcompoundnum:   ", hmdbcompoundnum)
        keggcompoundnum = sql.createRampCompoundID(kegg.metaboliteIDDictionary, "kegg", hmdbcompoundnum)
        wikicompoundnum = sql.createRampCompoundID(wikipathways.metaboliteIDDictionary, "wiki", keggcompoundnum)
        print("wikicompoundnum:   ", wikicompoundnum)
        reactomecompoundnum = sql.createRampCompoundID(reactome.metaboliteIDDictionary, "reactome", wikicompoundnum)
        
        print('Generate gene id ...')
        hmdbgenenum = sql.createRampGeneID(hmdb.geneInfoDictionary, "hmdb", 0)
        kegggenenum = sql.createRampGeneID(kegg.geneInfoDictionary, "kegg", hmdbgenenum)
        wikigenenum = sql.createRampGeneID(wikipathways.geneInfoDictionary, "wiki", kegggenenum)
        reactomegenenum = sql.createRampGeneID(reactome.geneInfoDictionary, "reactome", wikigenenum)
        print(" hmdbgenenum ", hmdbgenenum, " kegggenenum ", kegggenenum, " wikigenenum ", wikigenenum, " reactomegenenum ", reactomegenenum)
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
                 hmdb.metaInchi,
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
                 dict(),
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
                dict(),
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
                 dict(),
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
        
        
        # constructs the entity builder
        builder = EntityBuilder()
        
        # performs a full build of entities for loading
        # the input are files in /misc/output
        # the result are files for DB loading in /misc/sql
          
        builder.fullBuild()
        
        # Database loading is handled as a separate, un-coupled step.
            

resourceConfFile = "../../config/external_resource_config.txt"                
main = Main(resoureceConfFile)
main.runEverything()







import urllib.request
import os
import zipfile
from parse.MetabolomicsData import MetabolomicsData
from chemprop.ChemWrangler import ChemWrangler
from rampConfig.RampConfig import RampConfig

class lipidmapsChemData(MetabolomicsData):
    '''This class parses lipidmaps SDF file to capture lipidmaps structures specifically for the following tables
    source, met_classes and chem_props. The goal here is to increase chemical knowledge and expand ids.
    This data source does not conform to typical MetabolimicsData but is better suited specifically as a collection of Molecule entries.
    '''
    
    def __init__(self, resConfig):
        
        super().__init__()

        self.resourceConfig = resConfig
        
        self.moleculeList = list()
        
        self.sourceName = "LIPIDMAPS"

#         self.chemPropsDir = "../misc/data/chemprops/"
# 
        self.lipidMapsOutputDir = "../misc/output/lipidmaps/"
#         
        
        ####DICTIONARIES IN COMMON WITH OTHER CLASSES######################################
        # common name dictionary Key: HMDB ID Value: Common Name
        self.metaboliteCommonName = dict()
        #pathway dictionary. Key: hsaID for pathway, Value: pathway name
        self.pathwayDictionary = dict()
        # pathway id mapping Key SMP ID Value: Kegg id
        self.SMPToKegg = dict()
        #hsaID for pathway, value: pathway category (all will be "NA" for hmdb)
        self.pathwayCategory = dict()
        
        #key: metabolite id , value: list of pathway id
        self.metabolitesWithPathwaysDictionary = dict()
        
        #key: metabollite id, value: list of synonyms 
        self.metabolitesWithSynonymsDictionary = dict()
        
        #key: metabolite id, value: mapping of other ids
        self.metaboliteIDDictionary = dict()
        self.metaInchi = dict()
        # key: pathway SMP id, value: list of gene HMDBP id
        self.pathwaysWithGenesDictionary = dict()      
        # key: pathway id, value: list of metabolites id with this pathways.
        self.pathwaysWithMetabolitesDictionary = dict()
        #key: gene id, value: list of FOUR gene identifiers 
        self.geneInfoDictionary = dict()
        
        #only not empty when a catalyzed class exists 
        #key: matabole, value: list of genes
        self.metabolitesLinkedToGenes = dict()
        #self.inchiDict = dict()
        ###################################################################
        
        #stays empty for this class
        self.pathwayOntology = dict()
        
        #key: metabolite id, value: list of metabolite locations
        self.biofluidLocation = dict()
        
        #key: biofluid location, value: the string "placeholder"
        self.biofluid = dict()
        
        #key: metabolite id, value: list of cellular locations
        self.cellularLocation = dict()
        
        #key: cellular location, value: the string "placeholder"
        self.cellular = dict()
        
        #key: metaboliteID, value: exo/endo/Drug metabolite
        self.exoEndoDictionary = dict()
        
        #key: origin, value: exo/endo/ Drug ,etc.
        self.exoEndo = dict()
        #key: metaboliteID, value: tissue_locations
        self.tissueLocation = dict()
        
        #key: tissue location, value : "placeholder"
        self.tissue = dict()
        self.idDictForMetabolite = dict()
        # key: source ID e.g. HMDBID, value: dictionary that has key of sub,class,super class
        # value as the class name
        self.metaboliteClass = dict()


    def parseLipidMaps(self, writeToFile=False):
        chemist = ChemWrangler(self.resourceConfig)

        try:
            os.makedirs(self.lipidMapsOutputDir)
        except FileExistsError:
            # directory already exists
            pass

        # handled in get database files
        # chemist.fetchFile(self.sourceName, self.lipidMapsDir, "https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip", "LMSD.sdf.zip", "zip")

        metConfig = self.resourceConfig.getConfig("lipidmaps_met")
        localDir = metConfig.localDir
        metFile = metConfig.extractFileName

        chemist.readLipidMapsSDF(self.sourceName, localDir + metFile)
        lipidMapMolecules = chemist.chemLibDict[self.sourceName]

        if writeToFile:
            self.write_myself_files('lipidmaps')
            self.writeFiles(lipidMapMolecules, self.lipidMapsOutputDir)


    def writeFiles(self, molDict, lipidMapsOutputDir):
        
        print("Writing LipidMaps mol count=" + str(len(molDict)))
        
        classFile = "lipidmapsmetaboliteClass.txt"
        metIdFile = "lipidmapsmetaboliteIDDictionary.txt"
        commonNameFile = "lipidmapsmetaboliteCommonName.txt"
        synonymsFile = "lipidmapsmetaboliteSynonyms.txt"
      
        try:
            os.makedirs(lipidMapsOutputDir)
        except FileExistsError:
            # directory already exists
            pass
        
        classFileHandle  = open(lipidMapsOutputDir+classFile, "w+", encoding='utf-8')

        for id in molDict :
            classFileHandle.write(molDict[id].toClassString())
        
        classFileHandle.close()

        idFileHandle  = open(lipidMapsOutputDir+metIdFile, "w+", encoding='utf-8')

        for id in molDict :
            idFileHandle.write(molDict[id].toSourceString())
        
        idFileHandle.close()

        commonNameFileHandle  = open(lipidMapsOutputDir+commonNameFile, "w+", encoding='utf-8')

        for id in molDict :
            commonNameFileHandle.write(molDict[id].toCommonNameString())
        
        commonNameFileHandle.close()
        
        synonymsFileHandle  = open(lipidMapsOutputDir+synonymsFile, "w+", encoding='utf-8')

        for id in molDict :
            synonymsFileHandle.write(molDict[id].toSynonymsString())
        
        synonymsFileHandle.close()
        
    def getEverything(self, writeToFile = False):
        # get file resources
        self.getDatabaseFiles()
        self.parseLipidMaps(writeToFile)    
    
    def getDatabaseFiles(self):
        
        metConfig = self.resourceConfig.getConfig("lipidmaps_met")
                
        metFile = metConfig.sourceFileName
        metUrl = metConfig.sourceURL

        localDir = metConfig.localDir
    
        # check if this path exists
        self.check_path(localDir)
    
        self.download_files(metUrl, localDir + metFile)
        
        with zipfile.ZipFile(localDir+metFile,"r") as zip_ref:
            zip_ref.extractall(localDir)
    
# Test    
# resourceConfFile = "../../config/external_resource_config.txt" 
# resourceConf = RampConfig()
# resourceConf.loadConfig(resourceConfFile)
# lpData = lipidmapsChemData(resourceConf)
# lpData.getEverything(True)       
#         
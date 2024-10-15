import sys
import shutil
sys.path.append('../src')
from util.SQLiteDBBulkLoader import SQLiteDBBulkLoader



class mainSQLiteDBLoad(object):
    
    def __init__(self):
        # config for tables to load
        # a tab delimited file indicating which tables to load.
        self.dbConfigFilePath = "../config/db_load_resource_config.txt"
        
    
    
    def loadDBAfterTruncatingTables(self, schema_file='../schema/RaMP_SQLite_BASE.sqlite', incrementLevel = 'increment_patch_release', optionalVersionOveride = None, optionalVersionNote = None, truncateTables = False, tablesToKeep=['db_version', 'version_info']):

        if optionalVersionOveride:
            sqliteFile = schema_file.replace('BASE', f'v{optionalVersionOveride}')
            shutil.copy(schema_file, sqliteFile)
        ################# DB Loading Instructions
        # Sets logging level
        # config file holds login credentials in this format:
        # pass the credentials object to the constructed rampDBBulLoader
        
        loader = SQLiteDBBulkLoader(sqliteFileName=sqliteFile)


        # truncate tables
        if truncateTables:
            loader.truncateTables(tablesToSkip=tablesToKeep)
        
        
        # update methods
        # the sql_resource_config.txt is a tab delimited file indicating which resources to load
        # those marked as 'ready' will be updated. Usually all database tables are updated in one run.
        # this method loads the intermediate parsing results from the ../../misc/sql/ directory.
        loader.load(self.dbConfigFilePath)     
        
        # update Ontology Metabolite counts
        loader.updateOntologyMetaboliteCounts()
        
        # update Source pathwayCount
        loader.updateSourcePathwayCount()
        
        # sets the new updated version
        loader.updateDBVersion(incrementLevel = incrementLevel, optionalVersion = optionalVersionOveride, optionalNote = optionalVersionNote)
        
        # sets the analyte intercept json in the version table.
        # precondition: the updateDBVersion must have been set so that the
        # intersections can be attached to the current version
        loader.updateEntityIntersects()
        
        # this optional method tracks database version information supplied in this file.
        loader.updateVersionInfo("../config/ramp_resource_version_update.txt")
        
        # this method populates a table that reflects the current status of the database.
        # metrics such as gene and metabolite counts for reach data sets are tallied.
        loader.updateDataStatusSummary()

        # generate pathway similarity matrices, analyte lists and whatnot
        # this process replaced the old system of having Rdata in the package
        loader.generateAndLoadRampSupplementalData()

loader = mainSQLiteDBLoad()

# increment level 'increment_patch_release', 'increment_minor_release', 
# or 'specified' (new version, perhaps major release)
loader.loadDBAfterTruncatingTables(incrementLevel = 'specified',
       optionalVersionOveride = "2.6.4",
       optionalVersionNote = "20240822 data update, new datasource for pathways from PFOCR, updated MW check for mismerged metabolites, added field for best analyte name",
       truncateTables=True)

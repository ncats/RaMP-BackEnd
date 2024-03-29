import sys
sys.path.append('../src')
from util.SQLiteDBBulkLoader import SQLiteDBBulkLoader



class mainSQLiteDBLoad(object):
    
    def __init__(self):
        
        
        # db login credentials and host info
        # this is a private text file for login credentials
        # Format:
            # host=<host_uri>
            # dbname=<db_name_usually_ramp>
            # username=<db_user_name_often_root>
            # conpass=<db_connection_password>
        self.dbPropsFile = "../config/ramp_db_props.txt"
        
        # config for tables to load
        # a tab delimited file indicating which tables to load.
        self.dbConfigFilePath = "../config/db_load_resource_config.txt"
        
    
    
    def loadDBAfterTruncatingTables(self, sqliteFile, incrementLevel = 'increment_patch_release', optionalVersionOveride = None, optionalVersionNote = None, truncateTables = False, tablesToKeep=['db_version', 'version_info']):
        
        
        
    ################# DB Loading Instructions
        
        # Sets logging level
   
        # config file holds login credentials in this format:
 
        # pass the credentials object to the constructed rampDBBulLoader
        
        loader = SQLiteDBBulkLoader(dbPropsFile=self.dbPropsFile, sqliteFileName=sqliteFile)


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
loader.loadDBAfterTruncatingTables(sqliteFile = '/mnt/ncatsprod/braistedjc/tmp_work/RaMP_SQLite_BASE.sqlite', incrementLevel = 'specified',  
                                   optionalVersionOveride = "2.4.2", 
                                   optionalVersionNote = "20231107 data update, Rhea reaction to EC reaction class. Reactome Genes Patch.", 
                                   truncateTables=True)


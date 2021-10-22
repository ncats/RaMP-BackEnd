from util.rampDBBulkLoader import rampDBBulkLoader



class mainDBLoad():
    
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
        
    
    
    def loadDBAfterTruncatingTables(self, incrementLevel = 'increment_patch_release', optionalVersionOveride = None, optionalVersionNote = None):

    ################# DB Loading Instructions
        
        # Sets logging level
   
        # config file holds login credentials in this format:
 
        # pass the credentials object to the constructed rampDBBulLoader
        loader = rampDBBulkLoader(self.dbPropsFile)
        
        # update methods
        # the sql_resource_config.txt is a tab delimited file indicating which resources to load
        # those marked as 'ready' will be updated. Usually all database tables are updated in one run.
        # this method loads the intermediate parsing results from the ../../misc/sql/ directory.
        loader(dbConf, self.dbConfigFilePath)     
        
        loader.updateDBVersion(incrementLevel, optionalVersionOverride, optionalNote)
        
        
        # this optional method tracks database version information supplied in this file.
        #loader.updateVersionInfo("../../misc/resourceConfig/ramp_version_update_info.txt")
        
        # this method populates a table that reflects teh current status of the database.
        # metrics such as gene and metabolite counts for reach data sets are tallied.
        loader.updateDataStatusSummary()



loader = mainDBLoad()
loader.loadDBAfterTruncatingTables(incrementLevel = 'increment_patch_release', optionalVersionOveride = None, optionalVersionNote = None):


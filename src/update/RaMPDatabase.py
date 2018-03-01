import pymysql.cursors
from fileinput import close
import pandas as pd


class RaMPDatabase():
    '''
    This class is the super class of all checker, updater class 
    It contains general functions to aid genral functionality of other classes
    
    attribute str dbname database name for the mysql
     
    '''
    def __init__(self):
        
        self.table_names = [
                       "analyte",
                       "analytehasontology",
                       "analytehaspathway",
                       "analytesynonym",
                       "catalyzed",
                       "ontology",
                       "pathway",
                       "source"]
        
        
        
    def check_path(self,dir):
        '''
        This fucntion check if this directory exists, otherwise it will create one
        - param dict dir: The directory to check or created.
        - return: True if the path has been created successfully
        '''
        if not os.path.exists(dir):
            try:
                os.makedirs(dir) # check if the directory exists, create one if not
                return True
            except OSError as e: # Trap the OS error and show the embedded error code
                if e.errno != errno.EEXIST:
                    raise
        
    
    def connectToRaMP(self, host= "localhost", user = "root"
                      ,password = "Ehe131224",dbname = "mathelabramp"):
        '''
        Connect to local RaMP database by MySQL
        - param str host host name for the mysql connection
        - param str user username for the mysql conncection 
        - param str dbname database name for connection if None: connect to the database page
        instead of table page
        - param str password the password you used for you computer's mysql database 
        '''
        if dbname is not None:
            conn = pymysql.connect(host = host,
                                   user= user,
                                   password = password,
                                   db = dbname,
                                   charset = "utf8mb4",
                                   cursorclass = pymysql.cursors.DictCursor)
        else:
            conn = pymysql.connect(host = host,
                                   user= user,
                                   password = password,
                                   charset = "utf8mb4",
                                   cursorclass = pymysql.cursors.DictCursor)
        return conn     
        

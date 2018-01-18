import pymysql.cursors
from fileinput import close
import pandas as pd


class RaMPDatabase():
    '''
    This class is the super class of all checker, updater class 
    It contains general functions to aid genral functionality of other classes
    '''
    def __init__(self):
        self.name = "Hello World"
        self.tables = ["analyte",
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
        - dir the directory to check or created.
        - return True if the path has been created successfully
        '''
        if not os.path.exists(dir):
            try:
                os.makedirs(dir) # check if the directory exists, create one if not
                return True
            except OSError as e: # Trap the OS error and show the embedded error code
                if e.errno != errno.EEXIST:
                    raise
        
    
    def connectToRaMP(self):
        conn = pymysql.connect(host = "localhost",user="root",
                               password = "Ehe131224",
                               db = "mathelabramp",
                               charset = "utf8mb4",
                               cursorclass = pymysql.cursors.DictCursor)
        
        try: 
            with conn.cursor() as cur:
                df = pd.read_sql("select * from source;",conn)
                return df
                
        
        finally:
            conn.close()
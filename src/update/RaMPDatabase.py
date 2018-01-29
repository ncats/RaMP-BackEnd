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
        
        self.table_names = ["analyte",
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
        def create_new_db(self,database):
            '''
            This function will create a new database on MySQL
            if the database already existed, it will drop the previous one 
            '''
            conn = self.connectToRaMP(dbname = None)
            
            with conn.cursor() as cur:
                query = "create database "+database + ";"
                try:
                    cur.execute(query)
                    conn.commit()
                except pymysql.err.ProgrammingError:
                    query2 = "drop database mathelabramp2;"
                    cur.execute(query2)
                    conn.commit()
                    cur.execute(query)
                    conn.commit()  
                    
                finally:
                    conn.close()           
        
    def create_tbs(self,database):
        '''
        This function create tables based on the super class table_name attributes
        This part of sql code is exactly same as the C# part. If the table already existes,
        the table will be dropped.
        param str database the database name we want to connect
        '''        
        conn = self.connectToRaMP(dbname = "mathelabramp2")
        
        
        CreateTable = ["create table source(sourceId VARCHAR(30), rampId VARCHAR(30), IDtype VARCHAR(30), geneOrCompound VARCHAR(30),commonName VARCHAR(30), PRIMARY KEY (sourceId)) engine = InnoDB;",
                       "create table analyte (rampId VARCHAR(30), type VARCHAR(30), PRIMARY KEY (rampId)) engine = InnoDB;",
                       "create table pathway (pathwayRampId VARCHAR(30), sourceId VARCHAR(30),type VARCHAR(30), pathwayCategory VARCHAR(30), pathwayName VARCHAR(250), PRIMARY KEY (pathwayRampId)) engine = InnoDB;",
                       "create table analyteHasPathway (rampId VARCHAR(30), pathwayRampId VARCHAR(30), type VARCHAR(30)) engine = InnoDB;",
                       "create table analyteSynonym (Synonym VARCHAR(500), rampId VARCHAR(30), geneOrCompound VARCHAR(30), source VARCHAR(30)) engine = InnoDB;",
                       "create table catalyzed (rampCompoundId VARCHAR(30), rampGeneId VARCHAR(30)) engine = InnoDB;",
                       "create table analyteHasOntology (rampCompoundId VARCHAR(30), rampOntologyIdLocation VARCHAR(30)) engine = InnoDB;",
                       "create table ontology (rampOntologyIdLocation VARCHAR(30), commonName VARCHAR(30), biofluidORcellular VARCHAR(30)) engine = InnoDB;"]
        with conn.cursor() as cur:
            for sql in CreateTable:
                cur.execute(sql)
                conn.commit()
        
        conn.close()
        
    def drop_database(self,database = None):
        '''
        Drop the given database, call this function before you want to create new database
        param str database the database name you want to drop
        '''
        assert database != None, "No input! You should provide database name"
        assert type(database) is str, "The database name should be string"
        conn = self.connectToRaMP(dbname = None)
        with conn.cursor() as cur:
            sql = "drop database " + database +";"
            try:
                cur.execute(sql)
                conn.commit()
            finally:
                conn.close()
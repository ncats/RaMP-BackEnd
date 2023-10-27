'''
Created on Aug 2, 2023

@author: braistedjc
'''
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sklearn.metrics.pairwise import pairwise_distances
#from .rampDBBulkLoader import dbConfig

class RampSupplementalDataBuilder(object):
    '''
    classdocs
    '''


    def __init__(self, dbType, sqliteCreds=None, dbConf=None):
        '''
        Constructor
        '''
        # the type of DB, MySQL or SQLite
        self.dbType = dbType
        
        # a MySQL RaMP db_properties file, or an SQLite DB file 
        self.credInfo = sqliteCreds

        # sqlalchemy engine to provide connections to DB        
        self.engine = None
        
        if self.dbType == 'sqlite':
            self.engine = self.createSQLiteEngine(self.credInfo)
        else:
            self.engine = self.createMySQLEngine(dbConf)
        
        # all analyte pathway similarity matrix
        self.analyteResult = None

        # all analyte pathway similarity matrix        
        self.metsResult = None
        
        # all analyte pathway similarity matrix
        self.genesResult = None
        
    
    
    def createSQLiteEngine(self, sqliteFile=None):
        engine = create_engine('sqlite:///'+sqliteFile, echo=False)
        return engine

    def createMySQLEngine(self, dbConf = None):
        engine = create_engine((("mysql+pymysql://{username}:{conpass}@{host_url}/{dbname}").format(username=dbConf.username, conpass=dbConf.conpass, host_url=dbConf.host,dbname=dbConf.dbname)), echo=False)
        return engine
    
    
    def listTables(self):
        if self.dbType == 'mysql':
            sql = 'show tables'
        elif self.dbType == 'sqlite':
            sql = "SELECT name FROM sqlite_master WHERE type ='table' AND name NOT LIKE 'sqlite_%%'";
        else:
            print("Unsupported DB Type: " + self.dbType)
            return
        
        with self.engine.connect() as conn:
            tables = conn.execute(sql).all()            
            tables = pd.DataFrame(tables)
            print("tables shape:" + str(tables.shape))
            print(tables)
            conn.close()
            
    def buildPathwaySimilarityMatrices(self):
        x = None    

    def buildAnalyteSetStats(self):
        x = None

    def buildSimilarityMatrix(self, matrixType):
        df = None
        
        analyteKey = 'RAMP_%%'
        minPathwaySize = 10
        
        if matrixType == 'mets':
            analyteKey = 'RAMP_C%%'
            minPathwaySize = 5
        elif matrixType == 'genes':
            analyteKey = 'RAMP_G%%'
            minPathwaySize = 5
        
        sql = "select ap.pathwayRampId, ap.rampID from analytehaspathway ap, pathway p "\
        "where p.type != 'hmdb' and ap.pathwayRampId = p.pathwayRampId and ap.rampId like '" + analyteKey + "'"
        
        with self.engine.connect() as conn:
            df = conn.execute(sql).all()
            df = pd.DataFrame(df)
            df.columns = ['pathwayRampId', 'rampId']
            print(df.shape)
            print(list(df.columns))
            
            crossTab = pd.crosstab(df['rampId'], df['pathwayRampId'])
            ctSums = crossTab.sum(axis=0)
            pwSubset = ctSums[ctSums >= minPathwaySize]
            
            pwNames = pwSubset.index.values.tolist()
            crossTab = crossTab.loc[:,pwNames]
                        
            dm = 1.0 - pairwise_distances(crossTab.T.to_numpy(), metric='jaccard')


            dm = pd.DataFrame(dm)
                 
            dm.columns = crossTab.columns       
            dm.index = crossTab.columns

            conn.close()

        return dm



    def buildAnalyteSet(self, dataSource, geneOrMet):
        
        print("building analyte stat set")
        
        # NOTE the % has to be escaped for mysql, also works for sqlite, but is optional for sqlite.
        rampIdPrefix = "RAMP_C%%"
        if geneOrMet == 'gene':
            rampIdPrefix = "RAMP_G%%"
        
        sql = "select ap.pathwayRampId as pathwayRampId, count(distinct(ap.rampId)) as Freq, p.type as pathwaySource "\
        "from analytehaspathway ap, pathway p "\
        "where p.type = '" + dataSource + "' and ap.pathwayRampId = p.pathwayRampId and ap.rampId like '" + rampIdPrefix + "' group by ap.pathwayRampId"

        df = None

        with self.engine.connect() as conn:
            
            print(sql)
            
            df = conn.execute(sql).all()
            df = pd.DataFrame(df)
            
            print("Stats shape")
            print(df.shape)
            print("Stats header")
            print(df.columns)
            print(type(df))
            
            print(df.head(5))
            
            conn.close()

        return df

        
pwob = RampSupplementalDataBuilder(dbType = "sqlite", sqliteCreds = "X:\\braistedjc\\tmp_work\\RaMP_SQLite_v2.3.1b.sqlite")
#pwob.listTables()
dm = pwob.buildSimilarityMatrix(matrixType = "analytes")
print(dm.values.sum())

# pwob.buildSimilarityMatrix(matrixType = "genes")

#pwob.buildAnalyteSet("wiki", "met")
#pwob.buildAnalyteSet("wiki", "gene")

#pwob.buildAnalyteSet("reactome", "met")
#pwob.buildAnalyteSet("reactome", "gene")

#pwob.buildAnalyteSet("hmdb", "met")
# pwob.buildAnalyteSet("hmdb", "gene")
#pwob.buildBaseMatrix(matrixType = "genes")

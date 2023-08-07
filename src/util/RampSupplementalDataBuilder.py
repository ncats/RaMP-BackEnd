'''
Created on Aug 2, 2023

@author: braistedjc
'''
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sklearn.metrics.pairwise import pairwise_distances

class RampSupplementalDataBuilder(object):
    '''
    classdocs
    '''


    def __init__(self, dbType, credInfo):
        '''
        Constructor
        '''
        # the type of DB, MySQL or SQLite
        self.dbType = dbType
        
        # a MySQL RaMP db_properties file, or an SQLite DB file 
        self.credInfo = credInfo

        # sqlalchemy engine to provide connections to DB        
        self.engine = None
        
        if self.dbType == 'sqlite':
            self.engine = self.createSQLiteEngine(self.credInfo)
        
        # all analyte pathway similarity matrix
        self.analyteResult = None

        # all analyte pathway similarity matrix        
        self.metsResult = None
        
        # all analyte pathway similarity matrix
        self.genesResult = None
        
    
    
    def createSQLiteEngine(self, sqliteFile=None):
        engine = create_engine('sqlite:///'+sqliteFile, echo=False)
        return engine
    
    def listTables(self):
        if self.dbType == 'mysql':
            sql = 'show tables'
        elif self.dbType == 'sqlite':
            sql = "SELECT name FROM sqlite_master WHERE type ='table' AND name NOT LIKE 'sqlite_%'";
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
        
        analyteKey = 'RAMP_%'
        minPathwaySize = 10
        
        if matrixType == 'mets':
            analyteKey = 'RAMP_C%'
            minPathwaySize = 5
        elif matrixType == 'genes':
            analyteKey = 'RAMP_G%'
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
        
        rampIdPrefix = "RAMP_C%"
        if geneOrMet == 'genes':
            rampIdPrefix = "RAMP_G%"
        
        sql = "select ap.pathwayRampId, count(distinct(ap.rampId)) as Freq, p.type as pathwaySource "\
        "from analytehaspathway ap, pathway p "\
        "where p.type = '" + dataSource + "' and ap.pathwayRampId = p.pathwayRampId and ap.rampId like '" + rampIdPrefix + "' group by ap.pathwayRampId"

        df = None

        with self.engine.connect() as conn:
            df = conn.execute(sql).all()
            df = pd.DataFrame(df)
            
            print("Stats shape")
            print(df.shape)
            print("Stats header")
            print(df.columns)
            
            conn.close()

        return df

        
#pwob = PathwayOverlapBuilder(dbType = "sqlite", credInfo = "X:\\braistedjc\\tmp_work\\RaMP_SQLite_v2.3.0_Structure.sqlite")
#pwob.listTables()
#pwob.buildBaseMatrix(matrixType = "analytes")
# pwob.buildSimilarityMatrix(matrixType = "genes")

#pwob.buildAnalyteSet("wiki", "met")
#pwob.buildAnalyteSet("wiki", "gene")

#pwob.buildAnalyteSet("reactome", "met")
#pwob.buildAnalyteSet("reactome", "gene")

#pwob.buildAnalyteSet("hmdb", "met")
# pwob.buildAnalyteSet("hmdb", "gene")
#pwob.buildBaseMatrix(matrixType = "genes")

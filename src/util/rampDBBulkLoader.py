'''
Created on Aug 25, 2020

@author: braistedjc
'''
# from odo import odo
import mysql.connector
from mysql.connector import Error
import pandas as pd
from pandas.io.common import file_path_to_url
from pandas.api.types import is_string_dtype

import os.path
from os import path
from sqlalchemy import MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import insert
import logging

import pymysql.cursors


#from distutils.command.config import config
#from tornado.routing import _unquote_or_none
#from builtins import None

#from schema import RaMP_schema
from collections import namedtuple
from pprint import pprint
from pymysql import charset


#import pymysql
class dbConfig(object):
    def __init__(self):
        password = 'ramptest'
        username = 'ramp'
        host = 'ramp-db.ncats.io'
        dbname = 'ramp2'
        
        
class rampFileResource(object):

    def __init__(self):
        self.loadStatus = ""
        self.destTable = ""
        self.loadType = ""
        self.stagingFile = ""
        self.primaryKey = ""
        self.columnNames = []
        
    def initFileResource(self, resource):
        self.loadStatus = resource.status
        self.stagingFile = resource.file
        self.loadType = resource.loadType
        self.destTable = resource.table
        self.primaryKey = resource.primaryKey
        self.columnNames = resource.colNames.split(",")
        
    
    def printResource(self):
        pprint(vars(self))
                        
                        
class rampDBBulkLoader(object):
    
    def __init__(self):
        print("rampDBBulkLoaer__init__")
    
    # def create_connection(host_name, user_name, user_password, myDatabase):
    # 
    #     connection = None
    # 
    #     try:
    #         connection = mysql.connector.connect(
    #             host=host_name,
    #             user=user_name,
    #             passwd=user_password,
    #             database = myDatabase
    #         )
    #         
    # 
    #         print("Connection to MySQL DB successful")
    # 
    #     except Error as e:
    #         print(f"The error '{e}' occurred")
    #     return connection
    # 
    #     
    
    def remove_whitespace(self, dF):     
        for colName in dF.columns:
            if is_string_dtype(dF[colName]):
                dF[colName] = dF[colName].str.strip()
                print("fixing column...")
        return dF
                
            
        
    
    # lets load one file only....
    
    # def loadFile(conn, file_path, table, stmt, colIndices):
    #     print(os.path.abspath(os.getcwd()))
    #     print(str(path.exists(file_path)))
    #     data = pd.read_csv (file_path, sep='\t', header=None)
    #     df = pd.DataFrame(data)
    #     cursor = conn.cursor()
    #     cursor.executemany(stmt, df.iterrows())
    #     cursor.commit()
            
    def loadFile2(self, engine, file_path, table, stmt, colNames, colIndices, uniqueIndex = False, uniqueCol = None):
        print(os.path.abspath(os.getcwd()))
        print(str(path.exists(file_path)))
   
        data = pd.read_csv(file_path, delimiter=r'\t+', header=None, index_col=None, usecols=colIndices)
        df = pd.DataFrame(data)
    
        # issue with whitespace    
        df = self.remove_whitespace(df)
        
        # grab key columns append column names
        df = df[colIndices]     
        df.columns = colNames
        
        # should there be a unique key column? If so, keep the first instance
        # this came up with a redundant recorde in the intermediate file of a 1:1. 
        # final catch of a situation that can be avoided   
        if(uniqueCol):
            df = df.drop_duplicates(subset=[uniqueCol], inplace=False)
    
        # this loads the data frame into the table.
        try:
            df.to_sql(table, engine, if_exists = 'append', index=False)
        except Exception as ex:   
            print(ex)
            
        
        
#     def odoLoadFile(self, resource):
#         table = resource.destTable
#         file_path = "../../misc/sql/"+resource.stagingFile
#         url = ("mysql+pymysql://{user}:{pw}@ramp-db.ncats.io/{db}::{table}").format(user="ramp", pw="ramptest", db="ramp2", table="")
#         res = odo(file_path, url)
#         print(str(res))
        
    def loadFile(self, resource, engine):
        print(os.path.abspath(os.getcwd()))

        file_path = "../../misc/sql/"+resource.stagingFile
        colNames = resource.columnNames
        
        #data = pd.read_csv(file_path, sep="\t+", header=None, index_col=None, engine="python")
        data = pd.read_table(file_path, sep="\t", header=None, names=colNames, index_col=None, engine="python")     
        df = pd.DataFrame(data)
    
        # issue with whitespace    
        df = self.remove_whitespace(df)
        
        # grab key columns append column names
        
        print(df.shape)
        if(len(df.columns) > len(colNames)):
            print("too wide... cut down")           
            df = df.iloc[ : , 0:len(colNames)]
            print("new shape = "+str(df.shape))                   
        
        df.columns = colNames
        
        # should there be a unique key column? If so, keep the first instance
        # this came up with a redundant recorde in the intermediate file of a 1:1. 
        # final catch of a situation that can be avoided   
        
        print(str(df.shape))
        primaryKey = resource.primaryKey
        if(primaryKey != "None"):
            df = df.drop_duplicates(subset=[primaryKey], inplace=False, keep='first')
        else:
            # rationale.... final check if we have a completely duplicated record over all columns
            print("Droping full duplicates!")
            df = df.drop_duplicates(ignore_index=False, inplace=False, keep='first')
            print(str(df.shape))

    
            
    
        print(df.head(n=5))
        table = resource.destTable
        # this loads the data frame into the table.
        try:
            df.to_sql(table, engine, if_exists = 'append', index=False)
        except Exception as ex:   
            print(ex)
            pass
    
    
#         metaData = MetaData()
#         metaData.reflect(bind=engine)
#         myTable = metaData.tables[table]
#         __table_args__ = {'quote':False}
#         
#         for dataRow in df.to_dict(orient="records"):
#             #print(dataRow)
#             resource.columnNames;
#             
#             insert_stmt = insert(myTable).values(dataRow)
#            
# #             on_duplicate_key_stmt = insert_stmt.on_duplicate_key_update(
# #                 data=dataRow,
# #                 status='U')
#             
#             # print(str(insert_stmt.compile(engine)))
#             engine.execute(insert_stmt)
    

    def loadIgnore(self, engine, resource):
#         connection = pymysql.connect(host='ramp-db.ncats.io',
#                          user='ramp',
#                          password='ramptest',
#                          db='ramp2',
#                          charset='utf8mb4',
#                          cursorclass=pymysql.cursors.DictCursor)

        conn = mysql.connector.connect(host='ramp-db.ncats.io',
                         user='ramp',
                         password='ramptest',
                         db='ramp2',
                         charset = 'utf8',
                         use_unicode=True)
        #conn.set_charset_collation('utf16')
        cursor = conn.cursor()
        

        table = resource.destTable
        colNames = resource.columnNames
        fileName = resource.stagingFile
        
        cols = "`,`".join([str(i) for i in colNames])
        sql = "INSERT IGNORE INTO `"+table+"` (`" +cols+ "`) VALUES (" + "%s,"*(len(colNames)-1) + "%s)"
        #sql = "INSERT IGNORE INTO '"+table+" VALUES (" + "%s,"*(len(colNames)-1) + "%s)"
        #sql = "INSERT IGNORE INTO '"+table+" VALUES (%s)"
        valueString = "VALUES ("
        for col in colNames:
            valueString = valueString + "%("+col+")s,"
        valueString = valueString[:-1]    
        valueString = valueString + ")"
        #sql = "INSERT IGNORE INTO '"+table+"' (`" +cols+ "`) " + valueString
        
        print(sql)
        stmt = "insert into {} ({})".format(table, colNames)

        file_path = "../../misc/sql/"+fileName
        #data = pd.read_csv(file_path, sep="\t+", header=None, index_col=None, engine="python")
        data = pd.read_table(file_path, sep="\t+", header=None, index_col=None, engine="python", keep_default_na=False)     

        print(list(data.columns))

        df = pd.DataFrame(data)

        print(list(df.columns))

    
        # issue with whitespace    
        df = self.remove_whitespace(df)
        
        # grab key columns append column names
        colNames = resource.columnNames
        print(df.shape)
        if(len(df.columns) > len(colNames)):
            print("too wide... cut down")           
            df = df.iloc[ : , 0:len(colNames)]
            print("new shape = "+str(df.shape))       
        
        df.columns = colNames
        print("ColNames: "+str(colNames))
        # should there be a unique key column? If so, keep the first instance
        # this came up with a redundant recorde in the intermediate file of a 1:1. 
        # final catch of a situation that can be avoided   
        primaryKey = resource.primaryKey
        if(primaryKey != "None"):
            df = df.drop_duplicates(subset=[primaryKey], inplace=False, keep='first')
        else:
            # rationale.... final check if we have a completely duplicated record over all columns
            df = df.drop_duplicates(ignore_index=False, inplace=False, keep="first")
    # loadFile2(engine = sqlengine, file_path = "../../misc/sql/hmdbsource.sql", table = "source", 
    #          stmt= "INSERT INTO source (sourceId, rampId, IDtype, geneOrCompound, commonName) VALUES ( %s, %s, %s, %s, %s )",
    #          colNames = ['sourceId', 'rampId', 'IDtype', 'geneOrCompound', 'commonName'], colIndices = [0,1,2,3,4], uniqueIndex = True, uniqueCol='sourceId')
        
        records = list(df.itertuples(index=False))
        print(records[1])
      #  with connection.cursor() as cursor:        
        try:
            
            # cursor.executemany(sql, records)
            for i,row in df.iterrows():
#                 #print(tuple(row))
# #                 print(list(row.values()))
# #                 print(row.values())
# #                 values = tuple(row)
# #                 valuesStmt = " values ({})".format(df.iloc[i:i+1])                            
# #                 insertStmt = stmt + valuesStmt
#                 #print(insertStmt)
# #                 if(i > 1):
# #                     break
                cursor.execute(sql, tuple(row))
                if(i % 1000 == 0):
                    print(i)
                    conn.commit()

                
#            df.iloc[i:i+1].to_sql(table, engine, if_exists = 'append', index=False)
        except Exception as ex:   
            print(ex)
            pass
        
        conn.commit()
        #cursor.commit()
        #conn.commit()
    
    def loadConfig(self):
        print("nothing")

        
    def load(self, dbConf, rampResourceConfigFile):
    
        resourceConfig = pd.read_csv(rampResourceConfigFile, sep='\t', index_col=None)
        resourceConfig = pd.DataFrame(resourceConfig)
        resources = []
        fileResource = rampFileResource()
                
        for config in resourceConfig.itertuples():
            print(config.file)                    
            fileResource = rampFileResource()
            fileResource.initFileResource(config)            
            resources.append(fileResource)
            #fileResource.printResource(fileResource)
        
         
        engine = create_engine(("mysql+pymysql://{user}:{pw}@ramp-db.ncats.io/{db}").format(user="ramp", pw="ramptest", db="ramp2"), echo=False)
        # engine.dialect.identifier_preparer.initial_quote = "'"
        # engine.dialect.identifier_preparer.final_quote = "'"         

        print("Hey in loading loop now") 
        for resource in resources:
            print(resource.stagingFile)
            if(resource.loadStatus == "ready"):
                #self.loadFile(resource, engine)
                #self.odoLoadFile(resource)
                if(resource.loadType != "bulk"):
                    print("\n\nloadIgnore: "+resource.stagingFile+"\n\n")
                    self.loadIgnore(engine, resource)
                else:
                    self.loadFile(resource,engine)
                    print("\n\nbulkLoadFile: "+resource.stagingFile+"\n\n")


    def updateVersionInfo(self, infoFile):
        print("Updating Version Info")

        engine = create_engine(("mysql+pymysql://{user}:{pw}@ramp-db.ncats.io/{db}").format(user="ramp", pw="ramptest", db="ramp2"), echo=False)
        
        versionInfo = pd.read_csv(infoFile, sep='\t', index_col=None)
    
        #print(versionInfo)
    
        # change current status to archive
        sql = "update version_info set status = 'archive' where status = 'current'"

        with engine.connect() as conn:
            conn.execute(sql)
            
        versionInfo.to_sql("version_info", engine, if_exists = 'append', index=False)

        
    def updateDataStatusSummary(self):
        
        engine = create_engine(("mysql+pymysql://{user}:{pw}@ramp-db.ncats.io/{db}").format(user="ramp", pw="ramptest", db="ramp2"), echo=False)
        
        sqlMets = "select dataSource, count(distinct(rampId)) from source where geneOrCompound = 'compound' and dataSource not like '%%kegg' group by dataSource"
        sqlKeggMets = "select count(distinct(rampId)) from source where geneOrCompound = 'compound' and dataSource like '%%_kegg'"
        
        sqlGenes = "select dataSource, count(distinct(rampId)) from source where geneOrCompound = 'gene' and dataSource not like '%%kegg'group by dataSource"
        sqlKeggGenes = "select count(distinct(rampId)) from source where geneOrCompound = 'gene' and dataSource like '%%_kegg'"
        
        sqlPathways = "select type, count(distinct(pathwayRampId)) from pathway where sourceId not like 'map%%' group by type"
        sqlKeggPathways = "select count(distinct(pathwayRampId)) from pathway where sourceId like 'map%%'"
        
        sqlPathAssocMets = "select p.type, count(1) from pathway p, analytehaspathway a where a.rampId like 'RAMP_C%%' and p.pathwayRampId = a.pathwayRampId and p.sourceId not like 'map%%' group by p.type"
        sqlPathAssocGenes = "select p.type, count(1) from pathway p, analytehaspathway a where a.rampId like 'RAMP_G%%' and p.pathwayRampId = a.pathwayRampId and p.sourceId not like 'map%%' group by p.type"
        
        sqlKeggAssocMets = "select p.type, count(1) from pathway p, analytehaspathway a where a.rampId like 'RAMP_C%%' and p.pathwayRampId = a.pathwayRampId and p.sourceId like 'map%%' group by p.type"
        sqlKeggAssocGenes = "select p.type, count(1) from pathway p, analytehaspathway a where a.rampId like 'RAMP_G%%' and p.pathwayRampId = a.pathwayRampId and p.sourceId like 'map%%' group by p.type"

        sqlChemProps = "select chem_data_source, count(distinct(chem_source_id)) from chem_props group by chem_data_source"

        statusTable = dict()
        
        sourceNameDict = {'hmdb':'HMDB', 'kegg':'KEGG', 'lipidmaps':'Lipid Maps', 'reactome':'Reactome', 'wiki':'WikiPathways', 'chebi':'ChEBI'}
        
        with engine.connect() as conn:

            conn.execute("delete from entity_status_info")

            rs = conn.execute(sqlMets)
            statusTable["Metabolites"] = dict()
            for row in rs:
                statusTable["Metabolites"][row[0]] = row[1]
                
            rs = conn.execute(sqlKeggMets)
            for row in rs:
                statusTable["Metabolites"]['kegg'] = row[0]        

            rs = conn.execute(sqlGenes)
            statusTable["Genes"] = dict()
            for row in rs:
                statusTable["Genes"][row[0]] = row[1]
                
            rs = conn.execute(sqlKeggGenes)
            for row in rs:
                statusTable["Genes"]['kegg'] = row[0]     
            
            rs = conn.execute(sqlPathways)
            statusTable["Pathways"] = dict()
            for row in rs:
                statusTable["Pathways"][row[0]] = row[1]

            rs = conn.execute(sqlKeggPathways)
            for row in rs:
                statusTable["Pathways"]['kegg'] = row[0]            

            rs = conn.execute(sqlPathAssocMets)
            statusTable["Metabolite-Pathway Associations"] = dict()
            for row in rs:
                statusTable["Metabolite-Pathway Associations"][row[0]] = row[1]

            rs = conn.execute(sqlPathAssocGenes)
            statusTable["Gene-Pathway Associations"] = dict()
            for row in rs:
                statusTable["Gene-Pathway Associations"][row[0]] = row[1]

            rs = conn.execute(sqlKeggAssocMets)
            for row in rs:
                statusTable["Metabolite-Pathway Associations"]['kegg'] = row[1]    

            rs = conn.execute(sqlKeggAssocGenes)
            for row in rs:
                statusTable["Gene-Pathway Associations"]['kegg'] = row[1]    

            rs = conn.execute(sqlChemProps)
            statusTable["Chemical Property Records"] = dict()
            for row in rs:
                statusTable["Chemical Property Records"][row[0]] = row[1]  

            cols=('status_category','entity_source_id','entity_source_name','entity_count')
            df = pd.DataFrame(columns=cols)
            dataList = list()
            for cat in statusTable:
                for source in statusTable[cat]:
                    print(cat + " " + source + " " + str(statusTable[cat][source]))
                    row = [cat,source,sourceNameDict[source],statusTable[cat][source]]
                    print(row)
                    dataList.append(row)
             
            df = pd.DataFrame(dataList, columns=cols)        
            print(df)   
            
            df.to_sql("entity_status_info", engine, if_exists = 'append', index=False)
    

logging.basicConfig()
logging.getLogger('sqlalchemy').setLevel(logging.ERROR)
pd.set_option('display.max_columns', None)   
dbConf = dbConfig()
loader = rampDBBulkLoader()
#rampResourceConfigFile = "../../misc/resourceConfig/sql_resource_config.txt"
#loader.load(dbConf, rampResourceConfigFile)     

# loader.updateVersionInfo("../../misc/resourceConfig/ramp_version_update_info.txt")

loader.updateDataStatusSummary()

        




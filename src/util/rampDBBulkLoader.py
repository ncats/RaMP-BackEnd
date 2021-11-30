'''
Created on Aug 25, 2020

@author: braistedjc
'''
import mysql.connector
import pandas as pd
from pandas.api.types import is_string_dtype
import os.path
from os import path
from sqlalchemy import create_engine
from sqlalchemy import MetaData
import logging
from jproperties import Properties

class rampDBBulkLoader(object):

    def __init__(self, dbPropsFile):
        print("rampDBBulkLoaer__init__")
        
        # holds db credentials
        self.dbConf = dbConfig(dbPropsFile)
        
        logging.basicConfig()
        logging.getLogger('sqlalchemy').setLevel(logging.ERROR)
         
        pd.set_option('display.max_columns', None)   
     
   
    
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
            
        
    def loadFile(self, resource, engine):
        print(os.path.abspath(os.getcwd()))

        file_path = "../misc/sql/"+resource.stagingFile
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
    
    

    def loadIgnore(self, engine, resource):

        conn = mysql.connector.connect(host= self.dbConf.host,
                         user=self.dbConf.username,
                         password=self.dbConf.conpass,
                         db=self.dbConf.dbname,
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

        records = list(df.itertuples(index=False))
        print(records[1])
      #  with connection.cursor() as cursor:        
        try:
            
            # cursor.executemany(sql, records)
            for i,row in df.iterrows():
                cursor.execute(sql, tuple(row))
                if(i % 1000 == 0):
                    print(i)
                    conn.commit()

        except Exception as ex:   
            print(ex)
            pass
        
        conn.commit()

    
    def loadConfig(self):
        print("nothing")

        
    def load(self, rampResourceConfigFile):
    
        resourceConfig = pd.read_csv(rampResourceConfigFile, sep='\t', index_col=None)
        resourceConfig = pd.DataFrame(resourceConfig)
        resources = []
        fileResource = rampFileResource()
                
        for config in resourceConfig.itertuples():
            print(config.file)                    
            fileResource = rampFileResource()
            fileResource.initFileResource(config)            
            resources.append(fileResource)
            
        engine = create_engine(("mysql+pymysql://{username}:{conpass}@{host_url}/{dbname}").format(username=self.dbConf.username, conpass=self.dbConf.conpass, host_url=self.dbConf.host,dbname=self.dbConf.dbname), echo=False)

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

        engine = create_engine(("mysql+pymysql://{username}:{conpass}@{host_url}/{dbname}").format(username=self.dbConf.username, conpass=self.dbConf.conpass, host_url=self.dbConf.host,dbname=self.dbConf.dbname), echo=False)
       
        versionInfo = pd.read_csv(infoFile, sep='\t', index_col=None)

        # change current status to archive
        sql = "update version_info set status = 'archive' where status = 'current'"

        with engine.connect() as conn:
            conn.execute(sql)
            
        versionInfo.to_sql("version_info", engine, if_exists = 'append', index=False)

        
    def updateDataStatusSummary(self):
        
        print("starting update entity summary")
        
        engine = create_engine(("mysql+pymysql://{username}:{conpass}@{host_url}/{dbname}").format(username=self.dbConf.username, conpass=self.dbConf.conpass, host_url=self.dbConf.host,dbname=self.dbConf.dbname), echo=False)
        
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
                if row[0] > 0:
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
                    row = [cat,source,sourceNameDict[source],statusTable[cat][source]]
                    dataList.append(row)
             
            df = pd.DataFrame(dataList, columns=cols)        
            
            df.to_sql("entity_status_info", engine, if_exists = 'append', index=False)

            print("finished update entity summary")


    def updateDBVersion(self, incrementLevel = 'increment_patch_release', optionalVersion = None, optionalNote = None):
        
        self.dbConf.dumpConfig()
        
        engine = create_engine(("mysql+pymysql://{username}:{conpass}@{host_url}/{dbname}").format(username=self.dbConf.username, conpass=self.dbConf.conpass, host_url=self.dbConf.host,dbname=self.dbConf.dbname), echo=False)
        
        versionSQL = "select * from db_version where load_timestamp = (select max(load_timestamp) from db_version)"

        print("Updating DB Version")
                
        with engine.connect() as conn:
            
            meta_data = MetaData(bind=conn)
            meta_data.reflect()
            db_version = meta_data.tables['db_version']
            
            if incrementLevel != 'specified':
                                
                rp = conn.execute(versionSQL)
                res = rp.fetchone()
                                
                print("have current version of ramp: " + str(res['ramp_version']) + " " + str(res['load_timestamp']))

                vers = res['ramp_version']
                
                if incrementLevel == 'increment_patch_release':
                    end = vers.rfind(".")
                    start = end + 1                
                    newPatchLevel = (int)(vers[start:]) + 1                
                    newVersion = vers[0:start] + str(newPatchLevel)
                    print("increment patch release - new version = " + newVersion)
                elif incrementLevel == 'increment_minor_release':
                    end = vers.find(".") + 1
                    end2 = vers.rfind(".")
                    releaseVersion = (int)(vers[end:end2]) + 1
                    newVersion = vers[0:end] + str(releaseVersion) + ".0"
                    print("increment minor release - new version = " + newVersion)

            else:
                newVersion = optionalVersion
                print('set explicit db version = ' + newVersion)
            
            valDict = dict()
            valDict['ramp_version'] = newVersion
            
            if optionalNote is not None:
                valDict['version_notes'] = optionalNote

            vals = []
            vals.append(valDict)
            
            conn.execute(db_version.insert(), vals)            
            conn.close()
            
            

class dbConfig(object):
    
    def __init__(self, configFile):
                
        dbConfig = Properties()
        
        with open(configFile, 'rb') as config_file:
            dbConfig.load(config_file)
        
        self.conpass = dbConfig.get("conpass").data
        self.username = dbConfig.get("username").data
        self.host = dbConfig.get("host").data
        self.dbname = dbConfig.get("dbname").data
        
    def dumpConfig(self):        
        print(self.host)
        print(self.dbname)
        print(self.username)
        print(self.conpass)    
        
        
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
        
                        
        

# loader = rampDBBulkLoader("../../config/ramp_db_props.txt")
# 
# loader.updateDBVersion('increment_patch_release', None, "Testing the increment patch release")
# loader.updateDBVersion('increment_minor_release', None, "Testing the increment minor release")
# loader.updateDBVersion('specified', "v3.0.0", "Testing explicit version set")


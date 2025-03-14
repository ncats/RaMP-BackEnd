'''
Created on Aug 25, 2020

@author: braistedjc
'''
import pandas as pd
import numpy as np
from pandas.api.types import is_string_dtype
import os.path
from os import path
from sqlalchemy import create_engine, inspect, MetaData
import logging
from jproperties import Properties
import itertools
import time
from datetime import datetime
import json
from util.RampSupplementalDataBuilder import RampSupplementalDataBuilder

"""
from sqlalchemy import create_engine

engine = create_engine('sqlite://', echo=False)

I think building the engine is just a matter of supplying a file. This example above is in-memory.

The SQL file will need to have an existing DB.

"""
class SQLiteDBBulkLoader(object):

    def __init__(self, sqliteFileName):
        print("rampDBBulkLoaer__init__")
        
        self.sqliteFileName = sqliteFileName
        
        self.engine = self.createSQLiteEngine(sqliteFileName)
        
        self.currDBVersion = None

        self.sourceDisplayNames = {
            'kegg':'KEGG',
            'chebi':'ChEBI',
            'wikipathways_kegg':'KEGG',
            'hmdb_kegg':'KEGG',
            'hmdb':'HMDB',
            'reactome':'Reactome',
            'wiki':'WikiPathways',
            'lipidmaps':'LIPIDMAPS',
            'rhea':'Rhea',
            'pfocr': "PFOCR",
            'refmet': 'RefMet'
        }
        
        self.keggSubSources = ['hmdb_kegg', 'wikipathways_kegg']
        
        
        logging.basicConfig()
        logging.getLogger('sqlalchemy').setLevel(logging.ERROR)
         
        pd.set_option('display.max_columns', None)   
     
    
    def createSQLiteEngine(self, sqliteFile=None):
        engine = create_engine('sqlite:///'+sqliteFile, echo=False)
        return engine
    
    
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
    #     data = pd.read_csv (file_path, sep='\t', header=None, engine='python')
    #     df = pd.DataFrame(data)
    #     cursor = conn.cursor()
    #     cursor.executemany(stmt, df.iterrows())
    #     cursor.commit()
            
    def loadFile2(self, engine, file_path, table, stmt, colNames, colIndices, uniqueIndex = False, uniqueCol = None):
        print(os.path.abspath(os.getcwd()))
        print(str(path.exists(file_path)))
   
        data = pd.read_csv(file_path, delimiter=r'\t+', header=None, index_col=None, usecols=colIndices, engine='python')
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
        else:
            df = df.drop_duplicates(ignore_index=False, inplace=False, keep='first')
            
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
        if(primaryKey != "None" and not (isinstance(primaryKey, float) and np.isnan(primaryKey))):
            df = df.drop_duplicates(subset=[primaryKey], inplace=False, keep='first')
        else:
            # rationale.... final check if we have a completely duplicated record over all columns
            print("Dropping full duplicates!")
            df = df.drop_duplicates(ignore_index=False, inplace=False, keep='first')
            print(str(df.shape))

        table = resource.destTable
        # this loads the data frame into the table.
        try:
            df.to_sql(table, self.engine, if_exists = 'append', index=False)
        except Exception as ex:   
            print(ex)
            raise ex

    def loadConfig(self):
        print("nothing")

        
    def load(self, rampResourceConfigFile):
    
        resourceConfig = pd.read_csv(rampResourceConfigFile, sep='\t', index_col=None, engine='python')
        resourceConfig = pd.DataFrame(resourceConfig)
        resources = []
        fileResource = rampFileResource()
                
        for config in resourceConfig.itertuples():
            print(config.file)                    
            fileResource = rampFileResource()
            fileResource.initFileResource(config)            
            resources.append(fileResource)
            
        print("Hey in loading loop now") 
        for resource in resources:
            print(resource.stagingFile)
            if(resource.loadStatus == "ready"):
                #self.loadFile(resource, engine)
                #self.odoLoadFile(resource)
                self.loadFile(resource, self.engine)
                print("\n\nbulkLoadFile: "+resource.stagingFile+"\n\n")


    def updateVersionInfo(self, infoFile):
        print("Updating Version Info")

        sql = "select ramp_version, load_timestamp from db_version order by load_timestamp desc limit 1"
        
        dbVersion = None
        
        with self.engine.connect() as conn:
            dbVersionDF = conn.execute(sql).all()
            dbVersionDF = pd.DataFrame(dbVersionDF)
            conn.close()

        print(dbVersionDF.iloc[0,0])
        print(dbVersionDF.iloc[0,1])
        dbVersion = dbVersionDF.iloc[0,0]
        today = dbVersionDF.iloc[0,1]

        versionInfo = pd.read_csv(infoFile, sep='\t', index_col=None, header=0, engine='python')

        for c in versionInfo.columns:
            print("column: "+c)

        versionInfo['ramp_db_version'] = dbVersion
        versionInfo['db_mod_date'] = today

        print(versionInfo)

        # change current status to archive
        sql = "update version_info set status = 'archive' where status = 'current'"
  
        with self.engine.connect() as conn:
            conn.execute(sql)
            conn.close()
              
        versionInfo.to_sql("version_info", self.engine, if_exists = 'append', index=False)

        
    def updateDataStatusSummary(self):
        
        print("starting update entity summary")

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

        sqlReactAssocMets = "select count(*) from reaction2met"
        sqlReactAssocGenes = "select count(*) from reaction2protein"

        sqlHMDBGeneMetAssoc = "select count(*) from (select distinct rampCompoundId, rampGeneId from catalyzed)"
        sqlRheaGeneMetAssoc = "select count(*) from (select distinct rm.ramp_cmpd_id, rp.ramp_gene_id from reaction2met rm, reaction2protein rp where rm.ramp_rxn_id = rp.ramp_rxn_id)"

        statusTable = dict()

        inspector = inspect(self.engine)
        table_exists = "entity_status_info" in inspector.get_table_names()

        with self.engine.connect() as conn:
            if table_exists:
                conn.execute("delete from entity_status_info")

            rs = conn.execute(sqlReactAssocMets)
            statusTable["Metabolite-Reaction Associations"] = dict()
            for row in rs:
                statusTable["Metabolite-Reaction Associations"]['rhea'] = row[0]

            rs = conn.execute(sqlReactAssocGenes)
            statusTable["Gene-Reaction Associations"] = dict()
            for row in rs:
                statusTable["Gene-Reaction Associations"]['rhea'] = row[0]

            rs = conn.execute(sqlHMDBGeneMetAssoc)
            statusTable["Metabolite-Gene Associations"] = dict()
            for row in rs:
                statusTable["Metabolite-Gene Associations"]['hmdb'] = row[0]

            rs = conn.execute(sqlRheaGeneMetAssoc)
            for row in rs:
                statusTable["Metabolite-Gene Associations"]['rhea'] = row[0]


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
                    row = [cat,source,self.sourceDisplayNames[source],statusTable[cat][source]]
                    dataList.append(row)
             
            df = pd.DataFrame(dataList, columns=cols)        
            
            df.to_sql("entity_status_info", self.engine, if_exists = 'append', index=False)

            print("finished update entity summary")


    def updateDBVersion(self, incrementLevel = 'increment_patch_release', optionalVersion = None, optionalNote = None):
        
        #self.dbConf.dumpConfig()
        
        versionSQL = "select * from db_version where load_timestamp = (select max(load_timestamp) from db_version)"

        print("Updating DB Version")
                
        with self.engine.connect() as conn:
            
            meta_data = MetaData(bind=conn)
            meta_data.reflect()
            db_version = meta_data.tables['db_version']
            newVersion = ""
            
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

            valDict['load_timestamp'] = datetime.now()
            vals = []
            vals.append(valDict)
            
            conn.execute(db_version.insert(), vals)
            conn.close()
            
            self.currDBVersion = newVersion
            
    
    def updateEntityIntersects(self, filterComps=False):
        print("resolving analyte intersects")
        cmpdIntersects =  self.collectEntityIntersects(analyteType= 'compound', format='json', filterMets=filterComps)
        geneIntersects =  self.collectEntityIntersects(analyteType= 'gene', format='json', filterMets=filterComps)
        
        cmpdIntersectsInPW =  self.collectEntityIntersectsMappingToPathways(analyteType= 'compound', format='json', filterMets=filterComps)
        geneIntersectsInPW =  self.collectEntityIntersectsMappingToPathways(analyteType= 'gene', format='json', filterMets=filterComps)
                
        vals = list()
        vals.append({'met_intersects_json':cmpdIntersects, 'gene_intersects_json':geneIntersects, 'met_intersects_json_pw_mapped':cmpdIntersectsInPW, 'gene_intersects_json_pw_mapped':geneIntersectsInPW})
        
        if self.currDBVersion != None:
 
            with self.engine.connect() as conn:
                meta_data = MetaData(bind=conn)
                meta_data.reflect()
                db_version = meta_data.tables['db_version']
                conn.execute(db_version.update().where(db_version.c.ramp_version == self.currDBVersion).values(
                    met_intersects_json=cmpdIntersects, gene_intersects_json=geneIntersects, met_intersects_json_pw_mapped=cmpdIntersectsInPW, gene_intersects_json_pw_mapped=geneIntersectsInPW))
                #conn.execute(db_version.update().where(db_version.db_version == self.currDBVersion).values(vals))
            conn.close()
        print("updated DB analyte intersects")
    
            

    def collectEntityIntersectsMappingToPathways(self, analyteType='compound', format='json', filterMets=False, dropSMPD=False):
        sourceInfo = pd.read_table('../misc/sql/analytesource.txt', sep = '\t', header=None, dtype=str)
        sourceInfo = pd.DataFrame(sourceInfo)
        sourceInfo.columns = ['sourceId','rampId', 'idType', 'analyteType', 'commonName', 'status', 'dataSource']
        #sourceInfo.replace('hmdb_kegg', value='kegg', inplace=True)
        #sourceInfo.replace('wikipathways_kegg', value='kegg', inplace=True)
        print(sourceInfo.shape)
        mappingToPathways = pd.read_table('../misc/sql/analytetopathway.txt', sep = '\t', header=None, dtype=str)
        mappingToPathways = pd.DataFrame(mappingToPathways)
        mappingToPathways.columns = ['rampId', 'pathwayRampId', "dataSource"]
        print(mappingToPathways.shape)
#         mappingToPathways = pd.read_table('../../misc/sql/analytetopathway.txt', sep = '\t', header=None, dtype=str)
#         mappingToPathways = pd.DataFrame(mappingToPathways)
#         mappingToPathways.columns = ['rampId', 'pathwayRampId', "dataSource"]
#         print(mappingToPathways.shape)
        pathwayInfo = pd.read_table('../misc/sql/pathway.txt', sep = '\t', header=None, dtype=str)
        pathwayInfo = pd.DataFrame(pathwayInfo)
        pathwayInfo.columns = ['pathwayRampId','pathwayId','pathwaySource','pathwayCat', 'pathwayName']
        smpdbVersions = ['smpdb2', 'smpdb3']
        print(pathwayInfo.shape)
        #eliminate pathways mapping to smpdb
        pathwayInfo = pathwayInfo[~pathwayInfo['pathwayCat'].isin(smpdbVersions)]
        print(pathwayInfo.shape)
        #now restrict mappingToPathways to this subset
        mappingToPathways = mappingToPathways[mappingToPathways['pathwayRampId'].isin(pathwayInfo['pathwayRampId'])]
        print(mappingToPathways.shape)

        #here are the ramp ids that count
        rampIdsInPathways = mappingToPathways['rampId']
        print(rampIdsInPathways.shape)
        # limit ramp ids to pathway mapping, non-smpdb
        sourceInfo = sourceInfo[sourceInfo['rampId'].isin(rampIdsInPathways)]
        print(sourceInfo.shape)
        
        #  we used to collapse to display names here, but later moved this back in the process
        #  after we do more kegg accounting.
        
        #  for source in self.sourceDisplayNames:
        #      sourceInfo.replace(source, value=self.sourceDisplayNames[source], inplace=True)
        
        if filterMets:
            sourceInfo = sourceInfo[~sourceInfo['status'].isin(['predicted', 'expected'])]
        
        print(sourceInfo.shape)
        
        sourceInfo = sourceInfo[sourceInfo['analyteType'] == analyteType]
        sourceSet = set(sourceInfo['dataSource'])
        smSourceData = sourceInfo[['rampId', 'dataSource', 'analyteType']]            
        smSourceData = smSourceData.drop_duplicates()
        combos = []
        nodeList = list()
        print("source set")
        print(sourceSet)
        for r in range(1,len(sourceSet)+1):
            currCombos = itertools.combinations(sourceSet, r)
            #print("currCombos")
            #print(currCombos)
            combos += list(currCombos)
        
        intersectIndex = 0        
        
        #print("combos")
        #print(combos)
        
        for comb in combos:
            intersectIndex += 1
            combSet = set()
            combIndex = 0

            for s in comb:
                if combIndex == 0:
                    currCombSet = set(smSourceData.loc[smSourceData['dataSource'] == s]['rampId'])
                else:
                    currCombSet = currCombSet.intersection(set(smSourceData.loc[smSourceData['dataSource'] == s]['rampId']))
                    
                combIndex += 1
                                
            restData = smSourceData.loc[~smSourceData['dataSource'].isin(comb)]
            restSet = set(restData['rampId'])
            
            node = intersectNode()
            node.sets = list(comb)
            node.size = len(currCombSet-restSet)
            if(analyteType == 'compound'):
                node.id = "cmpd_src_set_" + str(intersectIndex)
            else:
                node.id = "gene_src_set_" + str(intersectIndex)
                                            
            if(node.size > 0):
                nodeList.append(node) 
        
        # need to deal with wikipathways_kegg and hmdb_kegg entities
        # if it's one of these, it should overlap hmdb or wiki
        # remove kegg subsource if only kegg subsource, add to the parent source tally
        print("node list length, initial: "+str(len(nodeList)))
        for node in nodeList:
            #print(node.id)
            #print(node.sets)
            #print(str(node.size))
            for s in node.sets:
                
                #### NOTE below here only happens if the current set has exactly 1 data source in the set and its either wikipathways_kegg or hmdb_kegg
                #### So this creates a node called 'mainSource' of the primary data source, *But currently isn't doing anything with 'mainSource' ???
                #### Then it seems to remove that single _kegg entry leaving an empty node.set, this is at a high level removing
                #### a node with one _kegg source in it's set. 
                #### then it traverses the entire nodeList until it finds a ndoe that has set size 2 and the two items are *_kegg and * where * is the primary source
                #### now we increment this particular joint node's count by the singleton _kegg node's size.
                #### This gives credit to the main source.
                
                #### This is ok because there should be no single-item set nodes that are *_kegg. So these are eliminated and the main source is incremented.                
                
                if len(node.sets) == 1 and s in self.keggSubSources:                                        
                    #mainSource = s.replace("_kegg", "")                
                    print("Node sets with kegg subsource: " + s)
                    if s == 'wikipathways_kegg':
                        nodeList.remove(node)
                        for n2 in nodeList:
                            if len(n2.sets) == 2 and 'wiki' in n2.sets and 'wikipathways_kegg' in n2.sets:
                                n2.size = n2.size + node.size
                    if s == 'hmdb_kegg':
                        nodeList.remove(node)
                        for n2 in nodeList:
                            if len(n2.sets) == 2 and 'hmdb' in n2.sets and 'hmdb_kegg' in n2.sets:
                                n2.size = n2.size + node.size
                                                            
        # drop in display names
        # 20220517 - correcting issue. We now have possibly more than one kegg node.
        for node in nodeList:
            sourceIndex = 0
            for source in node.sets: 
                node.sets[sourceIndex] = self.sourceDisplayNames[source] 
                sourceIndex = sourceIndex + 1  
        
        # Now we have to correct node sets that have two KEGG nodes. 
        # The node level accounting is the same, it's just that wikipathways_kegg and hmdb_kegg are both 
        # changed just to KEGG and now some nodes have two KEGG's listed.
        # So it's fair to drop a kegg entry within a node that has two keggs in it's set... but the set may or will become redundant, perhaps...
        for node in nodeList:
            keggCount = node.sets.count("KEGG")
            #print("keggCount = "+str(keggCount))
            if keggCount == 2:
                for s in node.sets:
                    if(s == "KEGG"):
                        node.sets.remove(s)
                        break
                #print("keggCountFixed = "+str(node.sets.count("KEGG")))
                #print(node.sets)
                
        # OK, redundant KEGG nodes have removed, BUT now we may have replicates sets
        # Hmmmm doesn't seem like we have redundant source nodes now. 
        nodesToRemove = list()
        touchedNodePairs = list()
        eqCount = 0
        for n in nodeList:
            for n2 in nodeList:
                if n.id != n2.id:
                    if (n.id+n2.id) not in touchedNodePairs and (n2.id+n.id) not in touchedNodePairs:
                        touchedNodePairs.append(n.id+n2.id)
                        touchedNodePairs.append(n2.id+n.id)
                        if n.sets == n2.sets:
                            eqCount = eqCount + 1
                            n.size = n.size + n2.size
                            if n2 not in nodesToRemove:
                                nodesToRemove.append(n2)
             
        for n in nodesToRemove:
            nodeList.remove(n)
                       

        # now we have fewer nodes left, lets add new ids:
        nid = 0
        for n in nodeList:
            nid = nid + 1
            if(analyteType == 'compound'):
                n.id = "cmpd_src_set_" + str(nid)
            else:
                n.id = "gene_src_set_" + str(nid)
        
        if format == 'json':
            jsonRes = json.dumps(nodeList, default=lambda o: o.__dict__, 
            sort_keys=True, indent=None)
            print(jsonRes)
            return jsonRes
        
        return nodeList
    
       
       
    def collectEntityIntersects(self, analyteType='compound', format='json', filterMets=False):
        sourceInfo = pd.read_table('../misc/sql/analytesource.txt', sep = '\t', header=None, dtype=str)
        sourceInfo = pd.DataFrame(sourceInfo)
        sourceInfo.columns = ['sourceId','rampId', 'idType', 'analyteType', 'commonName', 'status', 'dataSource']
        #sourceInfo.replace('hmdb_kegg', value='kegg', inplace=True)
        #sourceInfo.replace('wikipathways_kegg', value='kegg', inplace=True)
  
        
        # drop mapping to display name for now
#        for source in self.sourceDisplayNames:
#            sourceInfo.replace(source, value=self.sourceDisplayNames[source], inplace=True)
        
        if filterMets:
            sourceInfo = sourceInfo[~sourceInfo['status'].isin(['predicted', 'expected'])]
        
        
        sourceInfo = sourceInfo[sourceInfo['analyteType'] == analyteType]
        sourceSet = set(sourceInfo['dataSource'])
        smSourceData = sourceInfo[['rampId', 'dataSource', 'analyteType']]            
        smSourceData = smSourceData.drop_duplicates()
        combos = []
        nodeList = list()
        for r in range(1,len(sourceSet)+1):
            currCombos = itertools.combinations(sourceSet, r)
            combos += list(currCombos)
        
        intersectIndex = 0
        
        for comb in combos:
            intersectIndex += 1
            combSet = set()
            combIndex = 0

            for s in comb:
                if combIndex == 0:
                    currCombSet = set(smSourceData.loc[smSourceData['dataSource'] == s]['rampId'])
                else:
                    currCombSet = currCombSet.intersection(set(smSourceData.loc[smSourceData['dataSource'] == s]['rampId']))
                    
                combIndex += 1
                                
            restData = smSourceData.loc[~smSourceData['dataSource'].isin(comb)]
            restSet = set(restData['rampId'])
            
            node = intersectNode()
            node.sets = list(comb)
            node.size = len(currCombSet-restSet)
            if(analyteType == 'compound'):
                node.id = "cmpd_src_set_" + str(intersectIndex)
            else:
                node.id = "gene_src_set_" + str(intersectIndex)
            
            if(node.size > 0):
                nodeList.append(node) 
                
        # need to deal with wikipathways_kegg and hmdb_kegg entities
        # if it's one of these, it should overlap hmdb or wiki
        # remove kegg subsource if only kegg subsource, add to the parent source tally         
        for node in nodeList:
            #print(node.id)
            #print(node.sets)
            #print(str(node.size))
            for s in node.sets:
                if len(node.sets) == 1 and s in self.keggSubSources:                                        
                    mainSource = s.replace("_kegg", "")                    
                    #print("Node sets with kegg subsource: " + s)
                    if s == 'wikipathways_kegg':
                        nodeList.remove(node)
                        for n2 in nodeList:
                            if len(n2.sets) == 2 and 'wiki' in n2.sets and 'wikipathways_kegg' in n2.sets:
                                n2.size = n2.size + node.size
                    if s == 'hmdb_kegg':
                        nodeList.remove(node)
                        for n2 in nodeList:
                            if len(n2.sets) == 2 and 'hmdb' in n2.sets and 'hmdb_kegg' in n2.sets:
                                n2.size = n2.size + node.size
                                                            
        # drop in display names        
        for node in nodeList:
            sourceIndex = 0
            for source in node.sets: 
                node.sets[sourceIndex] = self.sourceDisplayNames[source] 
                sourceIndex = sourceIndex + 1              
#                source.replace(source, value=self.sourceDisplayNames[source], inplace=True)
        
        
        # Now we have to correct node sets that have two KEGG nodes. 
        # The node level accounting is the same, it's just that wikipathways_kegg and hmdb_kegg are both 
        # changed just to KEGG and now some nodes have two KEGG's listed.
        # So it's fair to drop a kegg entry within a node that has two keggs in it's set... but the set may or will become redundant, perhaps...
        for node in nodeList:
            keggCount = node.sets.count("KEGG")
            #print("keggCount = "+str(keggCount))
            if keggCount == 2:
                for s in node.sets:
                    if(s == "KEGG"):
                        node.sets.remove(s)
                        break
                #print("keggCountFixed = "+str(node.sets.count("KEGG")))
                #print(node.sets)
                
        # OK, redundant KEGG nodes have removed, BUT now we may have replicates sets
        # Hmmmm doesn't seem like we have redundant source nodes now. 
        nodesToRemove = list()
        touchedNodePairs = list()
        eqCount = 0
        for n in nodeList:
            for n2 in nodeList:
                if n.id != n2.id:
                    if (n.id+n2.id) not in touchedNodePairs and (n2.id+n.id) not in touchedNodePairs:
                        touchedNodePairs.append(n.id+n2.id)
                        touchedNodePairs.append(n2.id+n.id)
                        if n.sets == n2.sets:
                            eqCount = eqCount + 1
                            n.size = n.size + n2.size
                            if n2 not in nodesToRemove:
                                nodesToRemove.append(n2)
             
        for n in nodesToRemove:
            nodeList.remove(n)
                       
        # now we have fewer nodes left, lets add new ids:
        nid = 0
        for n in nodeList:
            nid = nid + 1
            if(analyteType == 'compound'):
                n.id = "cmpd_src_set_" + str(nid)
            else:
                n.id = "gene_src_set_" + str(nid)        

                
        if format == 'json':
            jsonRes = json.dumps(nodeList, default=lambda o: o.__dict__, 
            sort_keys=True, indent=None)
            print(jsonRes)
            return jsonRes
        
        return nodeList


    def create_column_if_not_exists(self, table, column):
        inspector = inspect(self.engine)
        columns = inspector.get_columns(table)
        column_exists = any(col['name'] == column for col in columns)
        if not column_exists:
            with self.engine.connect() as conn:
                conn.execute(f"ALTER TABLE {table} ADD COLUMN {column} INTEGER DEFAULT 0")

    def updateSourcePathwayCount(self):
        print("Started: updating pathway counts in source table")

        self.create_column_if_not_exists('source', 'pathwayCount')
    
        sql = "select count(distinct(ap.pathwayRampId)) as pathwayCount, ap.rampId from analytehaspathway ap "\
        "where ap.pathwaySource != 'hmdb' group by ap.rampId"
        
        sql2 = "update source set pathwayCount = :pathwayCount where rampId = :rampId"
        
        with self.engine.connect() as conn:
            df = conn.execute(sql).all()
            df = pd.DataFrame(df)
            df.columns = ["pathwayCount", "rampId"]

            print("setting pw count... shape=")
            print(df.shape)
            print(df.head(10))

            k = 0
            for i,row in df.iterrows():
                k = k + 1
                conn.execute(sql2, row)

            conn.close()

        print("Finished: updating pathway counts in source table")
    
    
    def updateOntologyMetaboliteCounts(self):
        print("Started: updating metabolite counts in ontology table")

        self.create_column_if_not_exists('ontology', 'metCount')

        sql = "select count(distinct(rampCompoundId)) as metCount, rampOntologyId from analytehasontology group by rampOntologyId"

        sql2 = "update ontology set metCount = :metCount where rampOntologyId = :rampOntologyId"
        
        with self.engine.connect() as conn:
            result = conn.execute(sql).all()
            df = pd.DataFrame(result)
            df.columns = ["metCount", "rampOntologyId"]
            
            for i,row in df.iterrows():
                conn.execute(sql2, row)
            conn.close()
    
        print("Finished: updating metabolite counts in ontology table")
        
    def updateCurrentDBVersionDumpURL(self, dumpUrl):
        self.dbConf.dumpConfig()
   
        print("Updating DB Version")
                
        with self.engine.connect() as conn:
            ts = conn.execute("select max(load_timestamp) from db_version")
            print(ts)
            ts = pd.DataFrame(ts)
            print(ts.shape)
            print(ts)
            print(ts.iloc[0,0])
            dbDumpURLUpdateSQL = "update db_version set db_sql_url = '"+dumpUrl+"' where load_timestamp = '"+str(ts.iloc[0,0])+"'"
            print(dbDumpURLUpdateSQL)
            conn.execute(dbDumpURLUpdateSQL)
            conn.close()
            
    def truncateTables(self, tablesToSkip):
        # self.dbConf.dumpConfig()
      
        print("Updating DB Version")
                
        with self.engine.connect() as conn:
            result = conn.execute("select name from sqlite_master where type = 'table' and name not like 'sqlite_%'")
            
            for row in result:
                print("delete from "+row[0])
                tableName = row[0]
                if tableName in tablesToSkip:
                    continue
                conn.execute("delete from " + tableName)
         
            conn.close()

    def generateAndLoadRampSupplementalData(self):
        start = time.time()

        dataBuilder = RampSupplementalDataBuilder(sqliteCreds = self.sqliteFileName)

        merged_matrix_dictionary = dataBuilder.getMergedSimilarityMatrix()
        duplicates = dataBuilder.getPathwaysWithSameAnalytes()

        sql = "insert into pathway_similarity (pathwayRampId, metabolite_count, gene_count, analyte_blob, metabolite_blob, gene_blob) VALUES (?,?,?,?,?,?)"
        dup_sql = "insert into pathway_duplicates (pathwayRampId1, pathwayRampId2) VALUES (?,?)"

        with self.engine.connect() as conn:
            conn.execute(sql, [
                (
                    key,
                    vals['metabolite_count'],
                    vals['gene_count'],
                    vals['analyte_blob'],
                    vals['metabolite_blob'],
                    vals['gene_blob']
                ) for key, vals in merged_matrix_dictionary.items()])

            conn.execute(dup_sql, duplicates)
            conn.close()

        end = time.time()
        elapsed_time = end - start
        print(f'elapsed_time: {elapsed_time} seconds')


class dbConfig(object):
    
    def __init__(self, configFile):
                
        dbConfig = Properties()
        
        with open(configFile, 'rb') as config_file:
            dbConfig.load(config_file)

        
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

class intersectList(object):

    def __init__(self):
        self.sets = list()
        self.size = 0
        self.id = "" 
      
class intersectNode(object):

    def __init__(self):
        self.sets = list()
        self.size = 0
        self.id = ""              

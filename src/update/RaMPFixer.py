from update.RaMPDatabase import RaMPDatabase
import pymysql.err
from builtins import str
import pandas as pd
import numpy as np
class RaMPFixer(RaMPDatabase):
    '''
    This class simulates the function of the C# code in RaMP
    Get the original RaMP data from MySQL, then do the following things:
    1) For each table, find the missing data cell, and remove the entire row since the map
    is not successful. (Some columns are deliberated to left empty with 'NA'.)
    2) Remove special character such as ' ', and wrong data in the cell. Consider drop entire
    row based on conditions
    3) Remapping corrected RaMP ID relations based on corrected data.
    4) Create and import new RaMP data to the database.
    '''
    
    def __init__(self):
        super().__init__()
        
        self.new_tables = dict()
        
        
    

    def remove_wrong_items_from_tb(self,table = None,Originaldb="mathelabramp"):
        '''
        This function removes error items from exited table.
        The errors include:
            1) Empty cell that is not reprented by None or 'NA'
            2) Wrong items that apparently should not be in here
        param str table string value for the table name
        param str Originaldb string value for the database that currently store the raw data
        ''' 
        assert table is not None and Originaldb is not None, "You should provide table name"
        assert type(table) is str and type(Originaldb) is str,"Table name should be a string"
        
        # Read sql to a dataframe here
        sql = "select * from {};".format(table)
        conn = self.connectToRaMP(dbname = Originaldb)
        df = pd.read_sql(sql,conn)
        print("Now dataframe has size {}".format(df.shape))
        df[["sourceId","rampId","IDtype","geneOrCompound"]] = df[["sourceId","rampId","IDtype","geneOrCompound"]].replace({None:np.NaN})
        df["commonName"] = df["commonName"].replace({None:"NA"})
        df = df.dropna()
        print("After Removing NA, the size is {}".format(df.shape))
        
        return df
from update.RaMPDatabase import RaMPDatabase
import pymysql.err
from builtins import str
import pandas as pd
import numpy as np
from schema import session,Source,Analyte,Pathway,Analytehasontology,\
Analytehaspathway,Analytesynonym,Ontology,Catalyzed

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
            
    


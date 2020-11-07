'''
Created on Nov 6, 2020

@author: braistedjc
'''
from pprint import pprint

class loadResource(object):
    '''
    classdocs
    '''
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


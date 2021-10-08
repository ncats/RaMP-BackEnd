'''
Created on Oct 8, 2021

@author: braistedjc
'''
import pandas as pd
from rampConfig.SourceConfig import SourceConfig

class RampConfig(object):
    '''
    RampConfig manages a collection of sourceConfig objects holding information for local and remote file resources
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.configDict = dict()
        
    def loadConfig(self, configFilePath):
        
        config = pd.read_table(configFilePath)
        
        print("reading resource config file: " + configFilePath)
        for i, confRow in config.iterrows():
            conf = SourceConfig()
            conf.resourceName, conf.sourceFetchMethod, conf.sourceURL, conf.sourceFileName, conf.extractFileName, conf.localDir, conf.compressType, conf.resourceType  = confRow
            self.configDict[conf.resourceName] = conf
            print("loading config:" + conf.resourceName)
        print("finished loading resource config")
        
    def getConfig(self, configKey):
        return self.configDict.get(configKey, None)    



        
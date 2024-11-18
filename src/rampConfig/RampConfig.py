'''
Created on Oct 8, 2021

@author: braistedjc
'''
from dataclasses import field

import pandas as pd
import yaml
from rampConfig.SourceConfig import SourceConfig

class RampConfig(object):
    '''
    RampConfig manages a collection of sourceConfig objects holding information for local and remote file resources
    '''
    configDict: dict = {}
    optionsDict: dict = {}

    def __init__(self, configFilePath, optionsFile = None):
        '''
        Constructor
        '''
        self.loadConfig(configFilePath=configFilePath)
        self.loadOptions(optionsFile=optionsFile)

    def loadOptions(self, optionsFile):
        if optionsFile is None:
            self.optionsDict = {}
        else:
            with open(optionsFile, "r") as file:
                self.optionsDict = yaml.safe_load(file)

    def getOptions(self, *args):
        configTree = self.optionsDict
        for arg in args:
            if arg in configTree:
                configTree = configTree[arg]
            else:
                return None
        return configTree

    def termIsOnOntologyDenyList(self, ontology, term):
        if not self.optionsDict:
            return False
        denyList = self.getOptions('ontology', 'denylist', ontology)
        return denyList is not None and term in denyList

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



        
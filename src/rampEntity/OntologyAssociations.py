'''
Created on Nov 6, 2020

@author: braistedjc
'''

class OnologyAssocations(object):
    '''
    Container class for pathway information.
    '''

    def __init__(self):
        
        self.sourceId
        
        self.origins = []
        
        self.biofluid = []
        
        self.cellularLoc = []
                
        self.tissueLoc = []
        
        self.helthEffect = []
        
        self.metApp = []
        
        self.termDict = dict()
        
        # pre-load theese into termDict to have central access
        self.termDict['Source'] = self.origins
        self.termDict['Biofluid and excreta'] = self.biofluid
        self.termDict['Subcellular'] = self.cellularLoc
        self.termDict['Tissue location'] = self.tissueLoc
        self.ontolDict["Organ and components"] = dict() 
        self.termDict['Industrial application'] = self.metApp
        self.termDict['Health condition'] = self.healthEffect
        
        
    def addOrigin(self, origin):        
        if origin not in self.origins:
            self.origins.append(origin)

    def addBiofluid(self, biofluid):        
        if biofluid not in self.biofluid:
            self.biofluid.append(biofluid)

    def addCellularLocation(self, cellLocation):        
        if cellLocation not in self.cellularLoc:
            self.cellularLoc.append(cellLocation)

    def addTissueLocation(self, tissue):        
        if tissue not in self.tissueLoc:
            self.tissueLoc.append(tissue)

    def addMetaboliteApplication(self, application):        
        if application not in self.metApp:
            self.metApp.append(application)
            
    def addHealthEffect(self, healthEffect):        
        if healthEffect not in self.healthEffect:
            self.helthEffect.append(healthEffect)        
        
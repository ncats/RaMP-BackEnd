import time
import pymysql
from sqlalchemy import *
from sqlalchemy_utils import *
from schema import RaMP_schema
class RampUpdater():
    def __init__(self,dbSource):
        '''
        dbSource is a data object from hmdbData, wikipathwayRDF, reactomeData, keggData 
        '''
        # key: ID for metabolites Value: Common Name (the only name in this database)
        self.metaboliteCommonName = dbSource.metaboliteCommonName
        #key: ID for pathway, Value: pathway name
        self.pathwayDictionary = dbSource.pathwayDictionary
        
        #key: ID for pathway, category: various categories such as cellular process, metabolic process 
        self.pathwayCategory = dbSource.pathwayCategory
        
        #key: gene, value: gene mapping
        self.geneInfoDictionary = dbSource.geneInfoDictionary
        
        #key: metabolite, value: metabolite mapping
        self.metaboliteIDDictionary = dbSource.metaboliteIDDictionary
        
        #key: pathwayID, value: list of genes
        self.pathwaysWithGenesDictionary = dbSource.pathwaysWithGenesDictionary
        
        #key: pathwayId, value: list of metabolites
        self.pathwayWithMetabolitesDictionary = dbSource.pathwayWithMetabolitesDictionary
        
        #empty for reactome
        self.metabolitesWithSynonymsDictionary = dbSource.metabolitesWithSynonymsDictionary
        
        #only not empty when a catalyzed class exists 
        #empty
        self.metabolitesLinkedToGenes = dbSource.metabolitesLinkedToGenes
        
        #key: metabolite, value: list of pathways
        self.metabolitesWithPathwaysDictionary = dbSource.metabolitesWithPathwaysDictionary
        
        #key: current pathway, value: list of pathways linked to this pathway
        self.pathwayOntology = dbSource.pathwayOntology
        
        #stays empty for this class
        self.biofluidLocation = dbSource.biofluidLocation
        
        #stays empty for this class
        self.biofluid = dbSource.biofluid
        
        #stays empty for this class
        self.cellularLocation = dbSource.cellularLocation
        
        #stays empty for this class
        self.cellular = dbSource.cellular
        
        #stays empty
        self.exoEndoDictionary = dbSource.exoEndoDictionary
        self.exoEndo = dbSource.exoEndo
        #tissue location stay empty
        self.tissue = dbSource.tissue
        self.tissueLocation = dbSource.tissueLocation
            
    def checkNewEntry(self):
        '''
        The data object should be from class hmdbData,KeggData,wikipathwayRDF,ReactomeData
        '''
        
        ramp_db = RaMP_schema()
        sess = ramp_db.session
        source_table = ramp_db.Source
        newCompound = set()
        for root,mapping in self.metaboliteIDDictionary.items():
            ids = self.getAllIDsFromMapping(mapping,root)
            (isoverlap, overlap) = self.isoverlap(ids)
            if isoverlap:
                rampid = self.findRampIDFromSource(overlap)
                #print('Root {} overlap with rampid {}'.format(root,rampid))
                #time.sleep(1)
            else:
                self.addNewRampIdToAnalytes('compound')
                time.sleep(3)
                newCompound.add(root)
        print('{} new compound from this update'.format(len(newCompound)))
        
        return 'Hello world'
    def isoverlap(self,listofids):
        sess = RaMP_schema().session
        sourcetb = RaMP_schema().Source
        isoverlap = False
        overlap = set()
        for each in listofids:
            (res,), = sess.query(exists().where(sourcetb.sourceId == each))
            if res:
                overlap.add(each)
        if len(overlap) > 0 :
            isoverlap = True
        
        return (isoverlap,overlap)
        
    def getAllIDsFromMapping(self,mapping,root):
        listofids = set()
        listofids.add(root)
        for source,ids in mapping.items():
            if ids is not 'NA':
                for each in ids:
                    listofids.add(each)
        
        return list(listofids)
    def findRampIDFromSource(self,setofids):
        sess = RaMP_schema().session
        sourcetb = RaMP_schema().Source
        rampids = set()
        for id in setofids:
            rampid = sess.query(sourcetb.rampId).filter(sourcetb.sourceId == id).all()
            rampids.add(rampid[0][0])
            
        return rampids
    def addNewRampIdToAnalytes(self,type):
        sess = RaMP_schema().session
        analytetb = RaMP_schema().Analyte
        rampid = sess.query(analytetb.rampId).filter(analytetb.type == type).order_by(desc(analytetb.rampId)).first()
        print(rampid[0])
        rampid = str(rampid[0])
        rampnumber = int(rampid[7:])
        
        print(rampnumber)
        
        #analyte = RaMP_schema().Analyte(rampId = rampid,type = type)
        #sess.add(analyte)
        #sess.commit()
        
        
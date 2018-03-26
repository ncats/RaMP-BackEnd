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
        # new rampid entry
        self.newRampCompound = dict()
        self.newRampGene = dict()
        self.newRampPathway = dict()
        self.newRampOntology = dict()
            
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
            print('Analytes {} is overlap? {}.'.format(root,overlap))
            if isoverlap:
                rampid = self.findRampIDFromSource(overlap)
                print('Root {} overlap with rampid {}'.format(root,rampid))
                if len(rampid) == 1:
                    self.oneidOverlap(rampid, ids, type = 'C')
                time.sleep(1)
            else:
                self.addNewRampIdToAnalytes('compound')
                time.sleep(3)
                newCompound.add(root)
        print('{} new compound from this update'.format(len(newCompound)))
        
        return 'Hello world'
    def isoverlap(self,listofids):
        '''
        Return true if the listofids has overlap ids with database
        - param list listofids a list of ids from source file or dictionary
        return:
        a 2-tuple (isoverlap, overlap) if isoverlap = True, the overlap will contain the id that has overlaps.
        '''
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
        '''
        Put all ids in mapping to a list.
        '''
        listofids = set()
        listofids.add(root)
        for source,ids in mapping.items():
            if ids not in ['NA','N','A']:
                for each in ids:
                    listofids.add(each)
        
        return list(listofids)
    def findRampIDFromSource(self,setofids):
        '''
        Find overlaped rampid from given source ids.
        '''
        
        sess = RaMP_schema().session
        sourcetb = RaMP_schema().Source
        rampids = set()
        for id in setofids:
            rampid = sess.query(sourcetb.rampId).filter(sourcetb.sourceId == id).all()
            rampids.add(rampid[0][0])
            
        return rampids
    def addNewRampIdToAnalytes(self,type = 'C'):
        assert type in ['C','G'], 'Please input a correct type, "G" for gene, "C" from compound'
        prefix = 'RAMP_' +type +'_'
        numberpart = 9
        sess = RaMP_schema().session
        analytetb = RaMP_schema().Analyte
        rampid = sess.query(analytetb.rampId).filter(analytetb.type == type).order_by(desc(analytetb.rampId)).first()
        print(rampid[0])
        rampid = str(rampid[0])
        rampnumber = int(rampid[7:])
        newid = prefix + '0' * (numberpart - len(str(rampnumber))) + str(rampnumber)
        print(newid)
        
        #analyte = RaMP_schema().Analyte(rampId = rampid,type = type)
        #sess.add(analyte)
        #sess.commit()
        
    def oneidOverlap(self,rampid,sourceids,type =None):
        assert len(rampid) == 1,'Call this function when you only has one rampid overlapped.'
        assert type in ['C','G'],'Analytes type should be C or G'
        sess = RaMP_schema().session
        sourcetb = RaMP_schema().Source
        otherids = sess.query(sourcetb.sourceId).filter(sourcetb.rampId == list(rampid)[0]).all()
        otherids = [i[0] for i in otherids if otherids is not None]
        totaloverlap = set(otherids).union(set(sourceids))
        print('The ids already in db is {}'.format(totaloverlap))
        if type is 'C':
            for each in totaloverlap:
                id_value = list(rampid)[0]
                if each not in self.newRampCompound:
                    self.newRampCompound[each] = id_value
                    
                
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
        self.pathwaysWithMetabolitesDictionary = dbSource.pathwaysWithMetabolitesDictionary
        
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
            
    def checkNewAnalyteEntry(self,analyte_type):
        '''
        The data object should be from class hmdbData,KeggData,wikipathwayRDF,ReactomeData
        - param str analyte_type string that specifies which type of analytes looking for.Should be 'compound' or 'gene' 
        '''
        ramp_db = RaMP_schema()
        # First get the last item from RaMP
        numberpart = 9
        sess = ramp_db.session
        analytetb = RaMP_schema().Analyte
        rampid = sess.query(analytetb.rampId).filter(analytetb.type == analyte_type).order_by(desc(analytetb.rampId)).first()
        rampid = str(rampid[0])
        # rampnumber is the number of analytes in a certain type (gene or compound) 
        rampnumber = int(rampid[7:])
        sess = ramp_db.session
        source_table = ramp_db.Source
        
        if analyte_type is 'compound':
            totalprocess = len(self.metaboliteIDDictionary)
            i = 0
            for root,mapping in self.metaboliteIDDictionary.items():
                i = i + 1
                if i % 100 is 0:
                    print('Processing {}/{}'.format(i,totalprocess))
                ids = self.getAllIDsFromMapping(mapping,root)
                (isoverlap, overlap) = self.isoverlap(ids)
                #print('Analytes {} is overlap? {}.'.format(root,isoverlap))
                if isoverlap:
                    rampid = self.findRampIDFromSource(overlap)
                    #print('Root {} overlap with rampid {}'.format(root,rampid))
                    if len(rampid) == 1:
                        # specific function to add entry that overlaps with only one ramp Id
                        self.oneidOverlap(rampid, ids, analyte_type = 'C')
                    elif len(rampid) > 1:
                        #print('Two rampId overlaped ... ')
                        disjointed_id = self.getDisjointIDsetFromRamp(rampid)
                        disjointed_id = disjointed_id.union(ids)
                        #print(disjointed_id)
                        # these disjointed set are unioned then assigned to a new ramp id. When writing the older one is dropped 
                        (rampnumber,newrampId) = self.findNewRampIdToAnalytes(rampnumber = rampnumber, analyte_type = analyte_type)
                        self.addAllidsToRamp(ids = disjointed_id, analyte_type='C', rampid = newrampId)
                        #time.sleep(2)
                else:
                    (rampnumber,newrampId) = self.findNewRampIdToAnalytes(analyte_type = analyte_type,rampnumber = rampnumber)
                    self.addAllidsToRamp(ids = ids, analyte_type = 'C', rampid = newrampId)
        elif analyte_type is 'gene':
            totalprocess = len(self.geneInfoDictionary)
            i = 0
            for root,mapping in self.geneInfoDictionary.items():
                i = i + 1
                if i % 100 is 0 :
                    print('Processing gene... {}/{}'.format(i,totalprocess))
                ids = self.getAllIDsFromMapping(mapping,root)
                (isoverlap, overlap) = self.isoverlap(ids)
                #print('Analytes {} is overlap? {}.'.format(root,isoverlap))
                if isoverlap:
                    rampid = self.findRampIDFromSource(overlap)
                    #print('Root {} overlap with rampid {}'.format(root,rampid))
                    if len(rampid) == 1:
                        # specific function to add entry that overlaps with only one ramp Id
                        self.oneidOverlap(rampid, ids, analyte_type = 'G')
                    elif len(rampid) > 1:
                        #print('Two rampId overlaped ... ')
                        disjointed_id = self.getDisjointIDsetFromRamp(rampid)
                        disjointed_id = disjointed_id.union(ids)
                        #print(disjointed_id)
                        # these disjointed set are unioned then assigned to a new ramp id. When writing the older one is dropped 
                        (rampnumber,newrampId) = self.findNewRampIdToAnalytes(rampnumber = rampnumber, analyte_type = analyte_type)
                        self.addAllidsToRamp(ids = disjointed_id, analyte_type='G', rampid = newrampId)
                else:
                    (rampnumber,newrampId) = self.findNewRampIdToAnalytes(analyte_type = analyte_type,rampnumber = rampnumber)
                    self.addAllidsToRamp(ids = ids, analyte_type = 'G', rampid = newrampId)
                    #time.sleep(3)
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
        - param set setofids a set containing all source ids that already have a rampId
        Return:
        The set of all rampIds that are from these source ids 
        '''
        sess = RaMP_schema().session
        sourcetb = RaMP_schema().Source
        rampids = set()
        for id in setofids:
            rampid = sess.query(sourcetb.rampId).filter(sourcetb.sourceId == id).all()
            rampids.add(rampid[0][0])
            
        return rampids
    def findNewRampIdToAnalytes(self,rampnumber,analyte_type = 'compound'):
        '''
        If the overlap set is an empty set, this function will return a new rampId based on the current condition of RaMP database
        - param string type a string that represents analytes type. It must be 'compound' or 'gene'.
        return:
            a new ramp Id as a string.
        '''
        assert analyte_type in ['compound','gene'], 'Please input a correct type, "G" for gene, "C" from compound'
        numberpart = 9 # Length of the number part of ramp ID
        analyte_type = 'C' if analyte_type is 'compound' else 'G'
        prefix = 'RAMP_' +analyte_type +'_'
        rampnumber = rampnumber + 1
        newid = prefix + '0' * (numberpart - len(str(rampnumber))) + str(rampnumber)
        #print(newid)
        return (rampnumber,newid)
        
    def oneidOverlap(self,rampid,sourceids,analyte_type =None):
        assert len(rampid) == 1,'Call this function when you only has one rampid overlapped.'
        assert analyte_type in ['C','G'],'Analytes type should be C or G'
        sess = RaMP_schema().session
        sourcetb = RaMP_schema().Source
        otherids = sess.query(sourcetb.sourceId).filter(sourcetb.rampId == list(rampid)[0]).all()
        # get source id from each tuple
        otherids = [i[0] for i in otherids if otherids is not None]
        totaloverlap = set(otherids).union(set(sourceids))
        #print('The ids already in db is {}'.format(totaloverlap))
        self.addAllidsToRamp(totaloverlap, analyte_type,rampid)
    
    def addAllidsToRamp(self,ids,analyte_type,rampid):
        '''
        This functions add all id in ids to the ramp dict of this class based on the given analyte_type
        - param set ids all ids to added to dictionary
        - param str analyte_type the string that represents analytes' type
        - param rampid the ramp id that will be assigned to them.
        '''
        if type(rampid) is set:
            id_value = list(rampid)[0]
        elif type(rampid) is str:
            id_value = rampid 
        oldrampid= id_value
        if analyte_type is 'C':
            for each in ids:
                # check if the entry overlaped with previous updating entries
                if each in self.newRampCompound:
                    id_value = self.newRampCompound[each]
                    break
            for each in ids:
                self.newRampCompound[each] = id_value
            
        elif analyte_type is 'G':
            for each in ids:
                # check if the entry overlaped with previous updating entries
                if each in self.newRampGene:
                    id_value = self.newRampGene[each]
                    break
            for each in ids:
                self.newRampGene[each] = id_value

    def getDisjointIDsetFromRamp(self,rampids):
        '''
        This functions is called when rampids is a set with size greater than 2.
        When the two or more disjointed id mapping is found. These source ids are extracted by this fucntion.
        '''
        sess = RaMP_schema().session
        sourcetb = RaMP_schema().Source 
        sourceids = sess.query(sourcetb.sourceId).filter(sourcetb.rampId.in_(rampids)).all()
        sourceids = {i[0] for i in sourceids if sourceids is not None}
        #print(sourceids)
        return sourceids           
    def checkNewPathwayEntry(self):
        ramp_db = RaMP_schema()
        # First get the last item from RaMP
        numberpart = 9
        sess = ramp_db.session
        pathwaytb = RaMP_schema().Pathway
        rpid = sess.query(pathwaytb.pathwayRampId).order_by(desc(pathwaytb.pathwayRampId)).first()
        #print('last rpid {}'.format(rpid))
        numberpart = int(rpid[0][7:])
        #print(numberpart)
        for pathwayid,name in self.pathwayDictionary.items():
            (res,), = sess.query(exists().where(pathwaytb.sourceId == pathwayid))
            if not res:
                #print('Pathway {}:{} in the database'.format(pathwayid,name))
                #print('Pathway {}:{} not in the database'.format(pathwayid,name))
                numberpart = numberpart + 1
                newrampid = 'RAMP_P_' + '0' * (9 - len(str(numberpart))) + str(numberpart)
                #print(newrampid)
                self.newRampPathway[pathwayid] = newrampid
    def writeToRamp(self,database):
        assert database in ['kegg','wiki','hmdb','reactome'],'Wrong database type specified.'
        # write to pathway first
        ramp_db = RaMP_schema()
        sess = ramp_db.session
        # initialize all useful tables
        pathwaytb = RaMP_schema().Pathway
        sourcetb = RaMP_schema().Source
        analytetb = RaMP_schema().Analyte
        pathwayhasanalytetb = RaMP_schema().Analytehaspathway
        # First, add pathways to the database
        all_pathway = []
        for pathwayid,rpid in self.newRampPathway.items():
            #print('Add pathway {}'.format(pathwayid))
            new_pathway = ramp_db.Pathway(pathwayRampId = rpid,
                                          sourceId = pathwayid,
                                          type = database,
                                          pathwayCategory = self.pathwayCategory[pathwayid],
                                          pathwayName = self.pathwayDictionary[pathwayid])

            all_pathway.append(new_pathway)
        print('Add {} pathways into RaMP'.format(len(all_pathway)))
        sess.add_all(all_pathway)
        sess.commit()
        del all_pathway
        
        for compoundid,rcid in self.newRampCompound.items():
            if ':' not in compoundid:
                continue
            (res_analyte,), = sess.query(exists().where(analytetb.rampId == rcid))
            (res_source,), = sess.query(exists().where(sourcetb.sourceId == compoundid))
            # not in both analyte and source table so it's new compound
            if not res_analyte and not res_source:
                print('Brand new metabolites {}:{}'.format(compoundid,rcid))
                new_analyte = ramp_db.Analyte(rampId = rcid,
                                type = 'compound'
                                )
                sess.add(new_analyte)
                sess.commit()
                if compoundid in self.metaboliteCommonName:
                    commonName = self.metaboliteCommonName[compoundid]
                else:
                    commonName = 'NA'
                new_source = ramp_db.Source(sourceId = compoundid,
                                            rampId = rcid,
                                            IDtype = compoundid[:compoundid.find(':')],
                                            geneOrCompound = 'compound',
                                            commonName = commonName
                                            )
                sess.add(new_source)
                sess.commit()
            
            # if in the source table but not in analyte table, it's disjointed idmapping 
            # Delete previous entry and update the new one
            elif not res_analyte and res_source:
                print('Old entry of row {}'.format(compoundid))
                sess.query(sourcetb).filter(sourcetb.sourceId == compoundid).delete()
                sess.commit()
                new_analyte = ramp_db.Analyte(rampId = rcid,
                                type = 'compound'
                                )
        
                sess.add(new_analyte)
                sess.commit()
                if compoundid in self.metaboliteCommonName:
                    commonName = self.metaboliteCommonName[compoundid]
                else:
                    commonName = 'NA'
                new_source = ramp_db.Source(sourceId = compoundid,
                                            rampId = rcid,
                                            IDtype = compoundid[:compoundid.find(':')],
                                            geneOrCompound = 'compound',
                                            commonName = commonName
                                            )
                sess.add(new_source)
                sess.commit()
            elif res_analyte and not res_source:
                print('New id {} with old entry'.format(compoundid))
                if compoundid in self.metaboliteCommonName:
                    commonName = self.metaboliteCommonName[compoundid]
                else:
                    commonName = 'NA'
                new_source = ramp_db.Source(sourceId = compoundid,
                                            rampId = rcid,
                                            IDtype = compoundid[:compoundid.find(':')],
                                            geneOrCompound = 'compound',
                                            commonName = commonName
                                            )
                sess.add(new_source)
                sess.commit()
        
        # writing gene to the Ramp
        print("##### Writing gene into ramp #####")

        for geneid,rgid in self.newRampGene.items():
            print('Adding {}/{}'.format(geneid,rgid))
            if ':' not in geneid:
                print('Wrong id')
                continue
            (res_analyte,), = sess.query(exists().where(analytetb.rampId == rgid))
            (res_source,), = sess.query(exists().where(sourcetb.sourceId == geneid))
            commonName = 'NA'
            # not in both analyte and source table so it's new compound
            if not res_analyte and not res_source:
                print('Brand new gene {}:{}'.format(geneid,rgid))
                new_analyte = ramp_db.Analyte(rampId = rgid,
                                type = 'gene'
                                )
                sess.add(new_analyte)
                sess.commit()
                
                IDtype = geneid[:geneid.find(':')]
                if IDtype is 'EN':
                    IDtype = 'enzymeNomenclature'
                new_source = ramp_db.Source(sourceId = geneid,
                                            rampId = rgid,
                                            IDtype = IDtype,
                                            geneOrCompound = 'gene',
                                            commonName = commonName
                                            )
                sess.add(new_source)
                sess.commit()
            elif not res_analyte and res_source:
                print('Old entry of row {}'.format(geneid))
                sess.query(sourcetb).filter(sourcetb.sourceId == geneid).delete()
                sess.commit()
                new_analyte = ramp_db.Analyte(rampId = rgid,
                                type = 'gene'
                                )
                sess.add(new_analyte)
                sess.commit()
                
                new_source = ramp_db.Source(sourceId = geneid,
                                            rampId = rgid,
                                            IDtype = geneid[:geneid.find(':')],
                                            geneOrCompound = 'gene',
                                            commonName = commonName
                                            )
                sess.add(new_source)
                sess.commit()
        # update analytehaspathway table
        # update metabolite-pathway relationship
        pathway_metabolite_list =[]
        for pathwayid,metaboliteList in self.pathwaysWithMetabolitesDictionary.items():
            rpid = sess.query(pathwaytb).filter(pathwaytb.sourceId == pathwayid).first()
            #print(rpid.pathwayRampId)
            
            linked_metabolite = sess.query(sourcetb).filter(sourcetb.sourceId.in_(metaboliteList)).all()
            linked_metabolite = {i.rampId for i in linked_metabolite}
            #print(linked_metabolite)
            for each in linked_metabolite:
                new_meta_path = pathwayhasanalytetb(rampId = each,pathwayRampId = rpid.pathwayRampId,
                                    pathwaySource = database)
                pathway_metabolite_list.append(new_meta_path)
        print('Add {} entries to analytehaspathway table'.format(len(pathway_metabolite_list)))
        sess.add_all(pathway_metabolite_list)
        sess.commit()
        pathway_metabolite_list = []
        for metaid, pids in self.metabolitesWithPathwaysDictionary.items():
            rcid = sess.query(sourcetb).filter(sourcetb.sourceId == metaid).first()
            linked_pathways = sess.query(pathwaytb).filter(pathwaytb.sourceId.in_(pids)).all()
            linked_pathways = [i.pathwayRampId for i in linked_pathways]
            for each in linked_pathways:
                new_meta_path = pathwayhasanalytetb(rampId = rcid.rampId,pathwayRampId = each,
                                                    pathwaySource = database)
                pathway_metabolite_list.append(new_meta_path)
        print('Add {} entries to analytehaspathway table'.format(len(pathway_metabolite_list)))
        sess.add_all(pathway_metabolite_list)
        sess.commit()
class QualityControl():
    def __init__(self):
        self.analytesIDtype = [
            'CAS',
            'chebi',
            'chemspider',
            'EnzymeNomenclature',
            'ensembl',
            'entrez',
            'hmdb',
            'kegg',
            'LIPIDMAPS',
            'pubchem',
            'uniprot'
            ]
        self.pathwaySource = [
            'hmdb',
            'kegg',
            'wiki',
            'reactome'
            ]
        self.analytesType = [
            'gene',
            'compound']
    def checkAllColumns(self):
        '''
        This function checks couple columns to see if there are unexpected bugs due to update
        It queries each table to get the unique value of each column. Find distinct value then compare it with our desired value.
        The function will also print the wrong items.
        '''
        db = RaMP_schema()
        sess = db.session
        sourcetb = db.Source
        pathwaytb = db.Pathway
        analytetb = db.Analyte
        sourceIDtype = sess.query(sourcetb.IDtype).distinct().all()
        print(sourceIDtype)
        pathwayIDtype = sess.query(pathwaytb.type).distinct().all()
        print( pathwayIDtype)
        analyteType1 = sess.query(analytetb.type).distinct().all()
        print(analyteType1)
        analyteType2 = sess.query(sourcetb.geneOrCompound).distinct().all()
        print(analyteType2) 
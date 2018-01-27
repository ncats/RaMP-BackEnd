import time
from MetabolomicsData import MetabolomicsData

from builtins import str

class writeToSQL(MetabolomicsData):
    
    '''This class takes the information gathered in database classes (such as hmdbData, keggData) and formats it
    properly for writing to .sql, which are used to create the mySQL database. 
    
    '''
    
    
    def __init__(self):
        super().__init__()
        #key: compound ID, value: rampID
        self.rampCompoundIDdictionary = dict()

        #key: gene ID, value: rampID
        self.rampGeneIDdictionary = dict()
        #key pathway ID value RampId 
        self.rampPathwayIDdictionary = dict()
        #key: RAMP_ID, value: set of databases the RAMPID is present in
        #After this dictionary is made it is not used in this class but will be used in the getStatistics class
        self.rampCompoundIdInWhichDatabases = dict()
        
        #key: RAMP_ID, value: set of databases the RAMPID is present in
        #After this dictionary is made it is not used in this class but will be used in the getStatistics class
        self.rampGeneIdInWhichDatabases = dict()

    
    def checkForWithinDatabaseDuplicatesCompound(self, metaboliteIDDictionary, database):
        
        '''
        The purpose of this function is to remove any duplicate metabolites present in the metaboliteIDDictionary. This may occur
        if it is listed in the database under different but synonymous IDs, for example, once as chebiID and once as HMDBID.
        
        param dict metaboliteIDDictionary: found in the four database classes, links metabolites to other IDs
        param str database: which database you are running this function on for example "hmdb" 
        
        '''
        duplicates = 0
        IDsDict = dict()
        
        listToPop = set()
         
        for key in metaboliteIDDictionary:
            
            
            mapping = metaboliteIDDictionary[key]
            listOfIds = set()
            for source in mapping:
                ids = mapping[source]
                if ids != 'NA':
                    if type(ids) is list:
                        for id in ids:
                            if id !='NA':
                                listOfIds.add(id)
                    else:
                        listOfIds.add(ids)
                        
            isnewCompound = False
            
            for each in listOfIds:
                if each not in IDsDict:
                    isnewCompound = True
                    
            
            if not isnewCompound:
                listToPop.add(key)
                print(key)
            
            if isnewCompound:
                for each in listOfIds:   
                    IDsDict[each] = key  
        duplicates = len(listToPop)
        print("There are "+ str(duplicates) +"items in " + database)            
        for each2 in listToPop:
            metaboliteIDDictionary.pop(each2)
            
    def checkForWithinDatabaseDuplicatesGene(self, geneInfoDictionary, database):
        
        '''
        
        The purpose of this function is to remove any duplicate genes present in the geneInfoDictionary. This may occur
        if it is listed in the database under different but synonymous IDs.
        
        param dict geneInfoDictionary: found in four database classes, links genes to other geneIDS
        param dict database: which database you are running this function on for example "hmdb" 
        '''
        
        IDsDict = dict()
        
        listToPop = set()
         
        for key in geneInfoDictionary:
            mapping = geneInfoDictionary[key]
            #print(geneInfoDictionary[key])
            #time.sleep(3)
            uniprotid = mapping['UniProt']
            hmdbgeneid = mapping['HMDB_protein_accession']
            entrez = mapping["Entrez"]
            enzymeNomenclature = mapping["Enzyme Nomenclature"]
            ensembl = mapping["Ensembl"]
            kegggeneid = mapping["kegg"]
            
           
            
            listOfIDs = set()
            if entrez is not "NA":
                listOfIDs.add(str(entrez))
            if enzymeNomenclature is not "NA":
                listOfIDs.add(enzymeNomenclature)
            if hmdbgeneid is not "NA":
                listOfIDs.add(hmdbgeneid)
            if ensembl is not "NA":
                for eachid in ensembl:
                    if ensembl is not "NA":
                        listOfIDs.add(eachid)
            if uniprotid is not "NA":
                for eachid in uniprotid:
                    if uniprotid is not "NA":
                        listOfIDs.add(eachid)
            if kegggeneid is not "NA":
                listOfIDs.add(kegggeneid)

            
            isnewGene = False
            
            for each in listOfIDs:
                if each not in IDsDict:
                    isnewGene = True
                    
            
            if not isnewGene:
                print("Going to remove: "+key)
                listToPop.add(key)
            
            if isnewGene:
                for each in listOfIDs:   
                    IDsDict[each] = key  

        for each2 in listToPop:
            geneInfoDictionary.pop(each2)
        
        
    def find_source_from_id(self,id,idmapping):
        for key,value in idmapping.items():
            if type(value) is list and id in value:
                return key
            else:
                if id == value:
                    return key
        return None        
             
    def connect_id_mapping(self,metaboliteIDdictionary,
                           entry_metabolite,overlap,
                           root_of_idmapping):
        '''
        if overlap between entry metabolites and ramp c id is found, these ids are considered same
        metabolites. This function looks for the ids that should have same ramp id
        param dict metaboliteIDdictionary metabolites ID mapping that has {source id: {source:other ids}} pair
        param dict entry_metabolite a id mapping that has {source:other ids} id pair
        param set overlap a set that contains overlapping ids
        param dict root_of_idmapping track which original id in original database this id is from
        return:
            set of ramp ids that should be changed their ramp cpd mapping
        '''
        # First compare both id and source
        rampid = set()
        filter_overlap =set()
        for id in overlap:
            origin = root_of_idmapping[id]
            if type(origin) is not list:
                mapping = metaboliteIDdictionary[origin]
                # find the source type from given overlap id
                source_of_entry = self.find_source_from_id(id, entry_metabolite)
                # Remove this id from overlap set becasue source is not same
                if source_of_entry is not None and mapping[source_of_entry] != id:
                    print('source of {} is not overlapped with previous entry'.format(id))
                elif source_of_entry is None:
                    raise ValueError('This id {} is not found from overlap id mapping'.format(id))
                else:
                    rampid.add(self.rampCompoundIDdictionary[id])
        
        
        
        return rampid
        
        
    def assign_rampC_id(self,list_of_mapping,rampid):
        '''
        code for assigning id if mapping is a list or mapping is a dict
        '''
        assert type(list_of_mapping) is list,'The mapping should be a list of metabolitesIDdict[key]'
        for mapping in list_of_mapping:
            if type(mapping) is not dict:
                raise TypeError('mapping should be one metabolitesIDdict') 
            for key,value in mapping.items():
                if value != 'NA':
                    if type(value) is str:
                        self.rampCompoundIDdictionary[value] = rampid
                    elif type(value) is list:
                        self.rampCompoundIDdictionary.update(dict.fromkeys(value, rampid))
        
        return
                
    def createRampCompoundID(self, metaboliteIDDictionary, database, rampCompoundIDnumber = 0):
        
        '''
        This function creates RAMPIDs for the compounds in the the RaMP database. This is a bit complicated of a process,
        since each RAMPID can have multiple "other" ids, spanning across multiple databases. 
        
        The self.rampGeneIDdictionary keeps track of this and prevents duplicates. 
        
        
        param dict metaboliteIDDictionary: found in the four database classes, links metabolites to other IDs
        param dict database: which database you are running this function on for example "hmdb" 
        param string rampCompoundIDnumber: what number to start the compounds on. If another database has already been run through this function, then use
        the function return value to know where to start (this function returns the number of the last compound run through the function)
        return string: returns the number of the last compound run through the function
        
        '''
        assert database in ['hmdb','kegg','reactome','wiki'],'database name should be exactly matched, please see createRampCompoundID function'
        rampCompoundID = "RAMP_C_000000000"
        lengthOfID = len(rampCompoundID)
        lengthOfIndex = len(str(rampCompoundIDnumber))
        prefix = lengthOfID - lengthOfIndex
        rampCompoundIDToFile = str(rampCompoundID[:prefix]) + str(rampCompoundIDnumber) 
    
        # key: source Id from {value} of metaboliteIDDict
        # value: dictionary that has {source database name : root id from the metabolites ID dictionary } 
        root_of_idmapping = dict()
        database_to_root ={'hmdb':'hmdb_id',
                           'kegg':'kegg_id',
                           'reactome':'chebi_id',
                           'wiki':['chebi_id','hmdb_id','kegg_id','CAS','pubchem_compound_id','chemspider_id']}
        for key in metaboliteIDDictionary:  
            # isThisNewCompound = False   
            listOfIDs = dict()
            mapping = metaboliteIDDictionary[key]
            
            # for each id != 'NA' in id mapping
            # they are added to the list
            for source,ids in mapping.items():
                if ids != 'NA':
                    if type(ids) is list:
                        for id in ids:
                            listOfIDs[id] = source
                    elif type(ids) is str:
                        listOfIDs[ids] = source
                    else:
                        raise TypeError('Wrong data type for the source id.')
            print('Creating ramp id for {}'.format(list(listOfIDs.keys())))  
            
            # each id is from id mapping     
            # This links each compoundID to a rampID but ONLY one rampID 
            # no matter if there is all three: chebi/hmdb/kegg
            # every id is the list
            # For relational database: unique item here is each id [source]
            
            # check if the set of ids overlaps with the previous entry
            # if not overlapping with any entry, this metabolite is considered a new metabolite
            # a new ramp compound id is therefore generated
            # set intersection gives shared ids
            overlap = set(listOfIDs.keys()) & set(self.rampCompoundIDdictionary.keys())
            
                
            
            if len(overlap) == 0 :
                # if there is no overlap, then it is a new compound
                rampCompoundIDnumber = rampCompoundIDnumber + 1
                for id in listOfIDs:
                    # keep record for all ids' root id
                    # root id is the key for original metabolite ID dict
                    lengthOfID = len(rampCompoundID)
                    lengthOfIndex = len(str(rampCompoundIDnumber))
                    prefix = lengthOfID - lengthOfIndex
                    rampCompoundIDToFile = str(rampCompoundID[:prefix]) + str(rampCompoundIDnumber)
                
                    self.rampCompoundIDdictionary[id] = rampCompoundIDToFile
                    
                #print('Assign this new compound {} ramp id {}'.format(key,rampCompoundIDnumber))    
                    
            else:
                print()
                notRealOverlap = set()   
                for id in overlap:
                    source = listOfIDs[id]
                    root_id_of_overlap_entry = root_of_idmapping[id]
                    for each in root_id_of_overlap_entry:
                        otherIDmapping = metaboliteIDDictionary[each]
                        if type(otherIDmapping[source]) is list and id not in otherIDmapping[source]:
                            notRealOverlap.add(id)
                        elif type(otherIDmapping[source]) is str and id is not otherIDmapping[source]:
                            notRealOverlap.add(id)
                
                # Remove the set of ids that are not really overlaped with previous entry
                overlap = overlap - notRealOverlap
                         
                # Find all ramp ID that are overlapped
                # inside this function, the overlap set will be checked both by id and source
                #print('Overlap is {} for {}'.format(overlap,key))

                # use the smallest ramp ID for all overlapping items
                # get all old and new ids that should be mapped together
                if len(overlap) == 0:
                    rampCompoundIDnumber = rampCompoundIDnumber + 1
                    #print('After source comparison, {} is a new compound'.format(list(listOfIDs.keys())))
                    for id in listOfIDs:
                        # keep record for all ids' root id
                        # root id is the key for original metabolite ID dict
                        lengthOfID = len(rampCompoundID)
                        lengthOfIndex = len(str(rampCompoundIDnumber))
                        prefix = lengthOfID - lengthOfIndex
                        rampCompoundIDToFile = str(rampCompoundID[:prefix]) + str(rampCompoundIDnumber)
                        if id not in self.rampCompoundIDdictionary:
                            self.rampCompoundIDdictionary[id] = rampCompoundIDToFile
                        else:
                            self.rampCompoundIDdictionary[listOfIDs[id]+':'+id] = rampCompoundIDToFile
                    #print('Assign this new compound {} ramp id {}'.format(key,rampCompoundIDnumber))
                    
                elif len(overlap) > 0:
                    ramp_id = set()
                    for each in overlap:
                        ramp_id.add(self.rampCompoundIDdictionary[each])
                    
                    # note string comparison here
                    #print('After source comparison, assign {} old ramp id {}'.format(key,ramp_id))
                    if len(ramp_id) == 1:
                        #time.sleep(1)
                        ramp_id = list(ramp_id)[0]
                        print('Old ramp id {} is assigned to {}'.format(ramp_id,key))
                        for id in listOfIDs:
                            if id not in self.rampCompoundIDdictionary:
                                self.rampCompoundIDdictionary[id] = ramp_id
                            else:
                                self.rampCompoundIDdictionary[listOfIDs[id]+':'+id] = ramp_id    
                    else:
                        print('Rare cases, this metabolite {} overlaps with more than two previous entry {}'.format(key,ramp_id))
                        time.sleep(10)
            for id in listOfIDs:
                if id not in root_of_idmapping:
                    root_of_idmapping[id] = [key]
                else:
                    root_of_idmapping[id].append(key)
           
           
    
        return rampCompoundIDnumber
    
    def createRampGeneID(self, geneInfoDictionary, database, rampGeneIDnumber = 0):
        
        '''
        This function creates RAMPIDs for the genes in the the RaMP database. This is a bit complicated of a process,
        since each RAMPID can have multiple "other" ids, spanning across multiple databases. 
        
        The rampCompoundIDdictionary keeps track of this and prevents duplicates. 
        
        param dict geneInfoDictionary: found in four database classes, links genes to other geneIDS
        param dict database: which database you are running this function on for example "hmdb" 
        param string rampGeneIDnumber: what number to start the compounds on. If another database has already been run through this function, then use
        the function return value to know where to start (this function returns the number of the last gene run through the function)
        return string: returns the number of the last gene run through the function
        
        '''
        
        rampGeneID = "RAMP_G_000000000"
        lengthOfID = len(rampGeneID)
        lengthOfIndex = len(str(rampGeneIDnumber))
        prefix = lengthOfID - lengthOfIndex
        rampGeneIDToFile = str(rampGeneID[:prefix]) + str(rampGeneIDnumber)

        
        for key in geneInfoDictionary:
            isThisNewGene = False    
            mapping = geneInfoDictionary[key]
            uniprotid = mapping['UniProt']
            hmdbgeneid = mapping['HMDB_protein_accession']
            entrez = mapping["Entrez"]
            enzymeNomenclature = mapping["Enzyme Nomenclature"]
            ensembl = mapping["Ensembl"]
            kegggeneid = mapping["kegg"]
            
           
            
            listOfIDs = []
            if entrez is not "NA":
                listOfIDs.append(str(entrez))
            if enzymeNomenclature is not "NA":
                listOfIDs.append(enzymeNomenclature)
            if hmdbgeneid is not "NA":
                listOfIDs.append(hmdbgeneid)
            if ensembl is not "NA":
                for eachid in ensembl:
                    if ensembl is not "NA":
                        listOfIDs.append(eachid)
            if uniprotid is not "NA":
                for eachid in uniprotid:
                    if uniprotid is not "NA":
                        listOfIDs.append(eachid)
            if kegggeneid is not "NA":
                listOfIDs.append("hsa:"+kegggeneid)
            
            
            for eachid in listOfIDs:
                if eachid not in self.rampGeneIDdictionary:
                    if not isThisNewGene:
                        isThisNewGene = True
                        
                        rampGeneIDnumber = rampGeneIDnumber + 1
                        lengthOfID = len(rampGeneID)
                        lengthOfIndex = len(str(rampGeneIDnumber))
                        prefix = lengthOfID - lengthOfIndex
                        rampGeneIDToFile = str(rampGeneID[:prefix]) + str(rampGeneIDnumber)
                        self.rampGeneIDdictionary[eachid] = rampGeneIDToFile
                        
                        setOfDatabases = set()
                        setOfDatabases.add(database)
                        self.rampGeneIdInWhichDatabases[rampGeneIDToFile] = setOfDatabases 
                        
                    else:
                        isThisNewGene = True
                        self.rampGeneIDdictionary[eachid] = rampGeneIDToFile
                      
                else:
                    OLDrampGeneIDToFile = self.rampGeneIDdictionary[eachid]
                    setOfCurrentDatabases = self.rampGeneIdInWhichDatabases[OLDrampGeneIDToFile]
                    if database not in setOfCurrentDatabases:
                        setOfCurrentDatabases.add(database)
                        self.rampGeneIdInWhichDatabases[OLDrampGeneIDToFile] = setOfCurrentDatabases  
                
                
                   

                    
                
                        
        rampGeneOut = open("../misc/output/" + str(database) + "rampGeneID.txt", 'wb')
        
        for key in self.rampGeneIDdictionary:
            rampGeneOut.write(str(key).encode("utf-8") + b":" + self.rampGeneIDdictionary[key].encode("utf-8")+b"\n")
        
        
        return rampGeneIDnumber
    
    def write_source(self,file,source_id,rampId,database,geneOrCompound,commonName):
        '''
        This functions write the input to the designated source file
        param _io.BufferedWriter file the file-like object that open the designated file
        param list|str source_id the source id of the analyte from original database (could be list or str)
        param str rampId rampId that is mapped with this source id
        param str database database name where the id is from
        param str geneOrCompound if this id is gene or compound
        param str commonName the common name that analyte has  
        
        '''
        assert geneOrCompound == "gene" or geneOrCompound == "compound","Wrong type of analytes"
        
        if type(source_id) is list:
            for id in source_id:
                file.write(id.encode("utf-8") +
                           b'\t'+rampId.encode('utf-8')+
                           b'\t'+database.encode('utf-8') +b'\t' + geneOrCompound.encode('utf-8') +
                           b'\t' + commonName.encode('utf-8') + b'\n')
        else:
            file.write(source_id.encode("utf-8") +
                       b'\t'+rampId.encode('utf-8')+
                       b'\t'+database.encode('utf-8') +b'\t' + 
                       geneOrCompound.encode('utf-8') +
                       b'\t' + commonName.encode('utf-8') + b'\n')        
        
            
        
        
    
    
    def write(self,
              metaboliteCommonName, 
              pathwayDictionary, 
              pathwayCategory,
              metabolitesWithPathwaysDictionary,
              metabolitesWithSynonymsDictionary,
              metaboliteIDDictionary,
              pathwaysWithGenesDictionary,
              metabolitesLinkedToGenes,
              geneInfoDictionary,
              biofluidLocation,
              biofluid,
              cellularLocation,
              cellular,
              pathwayOntology,
              exoEndoDictionary,
              exoEndo,
              tissueLocation,
              tissue,
              database, 
              rampPathwayIDnumber = 0,
              rampOntologyLocationIDnumber = 0):
        '''
        The function writeToFiles takes all the information gathered in the database and writes the required information to files.
        
        The information gathered in previous functions is stored in the dictionary objects that are passed to the function as parameters
        This information is formatted in this function and then written to files that can be used to create the new RaMP database. 
        
        param dict pathwayDictionary: see class for database 
        param dict pathwayCategory: see class for database
        param dict metabolitesWithPathwaysDictionary: see class for database
        param dict metabolitesWithSynonymsDictionary: see class for database
        param dict metaboliteIDDictionary: see class for database
        param dict pathwaysWithGenesDictionary: see class for database
        param dict geneInfoDictionary: see class for database
        param dict biofluidLocation: see class for database (may be empty for some classes)
        param dict biofluid: see class for database (may be empty for some classes)
        param dict cellularLocation: see class for database (may be empty for some classes)
        param dict cellular: see class for database (may be empty for some classes)
        param dict pathwayOntology: see class for database (may be empty for some classes)
        param dict endoExoDictionary: see class for database (may be empty for some classes)
        param dict endoExo: see class for database (may be empty for some classes)
        param dict tissueLocation: see class for database (may be empty for some classes)
        param dict tissue: see class for database (may be empty for some classes)
        param str  database: name of the database (e.g. "kegg") 
        param str rampPathwayIDnumber
        param str rampOntologyLocationIDnumber
        
        
        '''
        
        ###This section of the wrtieToFiles function maps every pathway, compound, and gene to a RAMPID that is pathway, compound, and gene specific 
        # ID dictionary
        # key :common name 
        # value : ramp ID
        rampPathwayIDdictionary = dict()
        rampBiofluidIDdictionary = dict()
        rampCellularIDdictionary = dict()
        rampExoEndoIDdictionary = dict() 
        rampTissueIDdictionary = dict()
        rampPathwayID = "RAMP_P_000000000"
        rampOntologyLocationID = "RAMP_OL_000000000"
                           

        
        for key in pathwayDictionary:
            rampPathwayIDnumber = rampPathwayIDnumber + 1
            lengthOfID = len(rampPathwayID)
            lengthOfIndex = len(str(rampPathwayIDnumber))
            prefix = lengthOfID - lengthOfIndex
            rampPathwayIDToFile = str(rampPathwayID[:prefix]) + str(rampPathwayIDnumber) 
            rampPathwayIDdictionary[key] = rampPathwayIDToFile
            
            #time.sleep(1)
            
        # output pathway id to see..
        '''
        pathwayFile = open('../misc/output/ramp/hmdbPathwayRampId.txt','wb')
        for key in rampPathwayIDdictionary:
            value = rampPathwayIDdictionary[key]
            for listitem in pathwaysWithGenesDictionary[key]:
                pathwayFile.write(key.encode('utf-8') +b'\t' +value.encode('utf-8') +b'\t'+listitem.encode('utf-8') +b'\n')
            
        pathwayFile.close()
        '''
        #return
        
            
        
        for key in biofluid:
            rampOntologyLocationIDnumber = rampOntologyLocationIDnumber + 1
            lengthOfID = len(rampOntologyLocationID)
            lengthOfIndex = len(str(rampOntologyLocationIDnumber))
            prefix = lengthOfID - lengthOfIndex
            rampOntologyIDToFile = str(rampOntologyLocationID[:prefix]) + str(rampOntologyLocationIDnumber) 
            rampBiofluidIDdictionary[key] = rampOntologyIDToFile
          
         
        for key in cellular:
            rampOntologyLocationIDnumber = rampOntologyLocationIDnumber + 1
            lengthOfID = len(rampOntologyLocationID)
            lengthOfIndex = len(str(rampOntologyLocationIDnumber))
            prefix = lengthOfID - lengthOfIndex
            rampOntologyIDToFile = str(rampOntologyLocationID[:prefix]) + str(rampOntologyLocationIDnumber)
            rampCellularIDdictionary[key] = rampOntologyIDToFile 
        
        for key in exoEndo:
            rampOntologyLocationIDnumber = rampOntologyLocationIDnumber + 1
            lengthOfID = len(rampOntologyLocationID)
            lengthOfIndex = len(str(rampOntologyLocationIDnumber))
            prefix = lengthOfID - lengthOfIndex
            rampOntologyIDToFile = str(rampOntologyLocationID[:prefix]) + str(rampOntologyLocationIDnumber)
            rampExoEndoIDdictionary[key] = rampOntologyIDToFile
        
        for key in tissue:
            rampOntologyLocationIDnumber = rampOntologyLocationIDnumber + 1
            lengthOfID = len(rampOntologyLocationID)
            lengthOfIndex = len(str(rampOntologyLocationIDnumber))
            prefix = lengthOfID - lengthOfIndex
            rampOntologyIDToFile = str(rampOntologyLocationID[:prefix]) + str(rampOntologyLocationIDnumber)
            rampTissueIDdictionary[key] = rampOntologyIDToFile
                        
        finalRAMPIDnumbers = [rampPathwayIDnumber, rampOntologyLocationID]
        
        self.check_path("../misc/sql/")
        analyteOutFile = open("../misc/sql/" + str(database) + "analyte.sql", 'wb')
        analyteSynonymOutFile = open("../misc/sql/" + str(database) + "analyteSynonym.sql", 'wb')                
        analyteHasPathwayOutFile = open("../misc/sql/" + str(database) + "analyteHasPathway.sql", 'wb')
        pathwayOutFile = open("../misc/sql/" + str(database) + "pathway.sql", 'wb')
        catalyzedOutFile = open("../misc/sql/" + str(database) + "catalyzed.sql", "wb")
        geneCrossLinksOutFile = open("../misc/sql/" + str(database) + "geneCrossLinks.sql", "wb")
        compoundCrossLinksOutFile = open("../misc/sql/" + str(database) + "compoundCrossLinks.sql", "wb")
        sourceOutFile = open("../misc/sql/" + str(database) + "source.sql", 'wb')
        ontologyLocationOutFile = open("../misc/sql/" + str(database) + "OntologyLocation.sql", "wb")
        analyteHasOntologyLocationOutFile = open("../misc/sql/" + str(database) + "analyteHasOntologyLocation.sql", "wb")
        pathwayOntologyOutFile = open("../misc/sql/" + str(database) + "PathwayOntology.sql", "wb")
        endoExoOutFile = open("../misc/sql/" + str(database) + "EndoExo.sql", "wb")
        
        #METABOLITE
        #analyte
        #analytehaspathway
        print("I'm analyte +analytehaspathway")
        for key in metabolitesWithPathwaysDictionary:
                value = metabolitesWithPathwaysDictionary[key]
                for listItem in value:
                    #This if statement is kinda a "hacky" fix...not sure why there is an empty key in this dictionary in the first place
                    if key is not "":
                        try:
                            analyteHasPathwayOutFile.write(str(self.rampCompoundIDdictionary[key]).encode('utf-8') + b"\t"
                                                                                         +  str(rampPathwayIDdictionary[listItem]).encode('utf-8') + b"\t"
                                                                                         + str(database).encode('utf-8') + b"\n") 
                        except KeyError:
                            print(str(KeyError) + " When writing analytehaspathways ...")
                            print(key)
        #GENE
        #analytehaspathway 
        print("Im analytehaspathway + Gene")                                                                                
        for key in pathwaysWithGenesDictionary:
            value = pathwaysWithGenesDictionary[key]
            for listItem in value:

                try:
                    analyteHasPathwayOutFile.write(str(self.rampGeneIDdictionary[listItem]).encode('utf-8') + b"\t" +  str(rampPathwayIDdictionary[key]).encode('utf-8') + b"\t" + str(database).encode('utf-8') + b"\n")

                except KeyError:
                    print(str(KeyError) + " when writing genes has pathways ..." +listItem)

                    
        #GENE
        #analyte
        for key in geneInfoDictionary:
            
            if key in self.rampGeneIDdictionary:
                analyteOutFile.write(self.rampGeneIDdictionary[key].encode('utf-8') + b"\t" 
                                                                 + b"gene" + b"\n")
        
        for key in metaboliteIDDictionary: 
            #This if statement is kinda a "hacky" fix...not sure why there is an empty key in this dictionary in the first place
            if key is not "":
                analyteOutFile.write(self.rampCompoundIDdictionary[key].encode('utf-8') + b"\t" 
                                                                 + b"compound" + b"\n")
        
        
        
        #METABOLITE
        #source
        for key in metaboliteIDDictionary:    
            mapping = metaboliteIDDictionary[key]
                        #there are multiple chebi for every compound so you cannot 
            #just call mapping["chebi_id"]
            commonName = metaboliteCommonName[key]
            if commonName is None:
                commonName = "NA"
            id_to_write = {'chebi_id':'chebi','hmdb_id':'hmdb',
                           'kegg_id':'kegg','CAS':'CAS',
                           'pubchem_compound_id':'pubchem',
                           'chemspider_id':'chemspider',
                           'LIPIDMAPS':'LIPIDMAPS'}
            for id_key in id_to_write:
                if mapping[id_key] is not 'NA':
                    self.write_source(sourceOutFile, 
                                      mapping[id_key], 
                                      self.rampCompoundIDdictionary[key], 
                                      id_to_write[id_key], 
                                      'compound', 
                                      commonName)
        #GENE
        #Source
        for key in geneInfoDictionary:
            mapping = geneInfoDictionary[key]
            #key = key.replace("HMDBP","HMDB")
            uniprotid = mapping['UniProt']
            hmdbgeneid = mapping['HMDB_protein_accession']
            entrez = mapping["Entrez"]
            enzymeNomenclature = mapping["Enzyme Nomenclature"]
            ensembl = mapping["Ensembl"]
            kegggeneid = mapping["kegg"]
            commonName = mapping["common_name"]
            NameForSource = None
            if type(commonName) is not list:
                commonName = commonName.replace("\n", "")
                commonName = commonName.replace("\"", "")
                commonName = commonName.replace(" ", "")
                if commonName is not "NA":
                    NameForSource = commonName
                else:
                    NameForSource = "NA"
                        
            else:
                item = commonName[0]
                item = item.replace("\n", "")
                item = item.replace("\"", "")
                item = item.replace(" ", "")
                if item is not "NA":
                    NameForSource = item
                else:
                    NameForSource = "NA"
                            
            if key in self.rampGeneIDdictionary:
                if uniprotid is not "NA":
                    for eachid in uniprotid:
                        if eachid is not "NA":
                            sourceOutFile.write(eachid.encode('utf-8') + b"\t" 
                                                + self.rampGeneIDdictionary[key].encode('utf-8') + 
                                                b"\t" + b"uniprot" + b"\t" + b"gene"
                                                + b"\t" + NameForSource.encode("utf-8") + b"\n")
                
                if hmdbgeneid is not "NA":
                    sourceOutFile.write(hmdbgeneid.encode('utf-8') + b"\t" +
                                         self.rampGeneIDdictionary[key].encode('utf-8') +
                                          b"\t" + b"hmdb" + b"\t" + b"gene"+
                                          b"\t" + NameForSource.encode("utf-8") + b"\n")
                    
                if entrez is not "NA":
                    sourceOutFile.write(str(entrez).encode('utf-8') + b"\t" + 
                                        self.rampGeneIDdictionary[key].encode('utf-8')
                                         + b"\t" + b"entrez" + b"\t" + b"gene"
                                         + b"\t" + NameForSource.encode("utf-8") + b"\n")
                    
                    
                    
                if enzymeNomenclature is not "NA":
                    sourceOutFile.write(enzymeNomenclature.encode('utf-8') + 
                                        b"\t" + self.rampGeneIDdictionary[key].encode('utf-8') 
                                        + b"\t" + b"enzymeNomenclature" + b"\t" + b"gene"
                                         + b"\t" + NameForSource.encode("utf-8") + b"\n")
                    
                if ensembl is not "NA":
                    for eachid in ensembl:
                        if eachid is not "NA":
                            sourceOutFile.write(eachid.encode('utf-8') + b"\t"
                                                 + self.rampGeneIDdictionary[key].encode('utf-8')
                                                  + b"\t" + b"ensembl" + b"\t" + b"gene"+
                                                    b"\t" + NameForSource.encode("utf-8") + b"\n")
                    
                if kegggeneid is not "NA":
                    if type(kegggeneid) is not list:
                        #kegggeneid = "hsa:" + kegggeneid
                        sourceOutFile.write(kegggeneid.encode('utf-8')
                                             + b"\t" + self.rampGeneIDdictionary[key].encode('utf-8') 
                                             + b"\t" + b"kegg" + b"\t" + b"gene"
                                              + b"\t" + NameForSource.encode("utf-8") + b"\n")
                       
            else:
                print("The kegg gene {} does not have Ramp Gene Id".format(key))            
                #time.sleep(0.1)
        #METABOLITE
        #metabolite with synonym (synonyms for the common name)
        for key in metabolitesWithSynonymsDictionary:
            value = metabolitesWithSynonymsDictionary[key]
            for listItem in value:
                listItem = listItem.replace("\n", "")
                listItem = listItem.replace("\"", "")
                listItem = listItem.replace(" ", "")
                listItem = listItem.replace(";", "")
                listItem.lower()
                if listItem is not "NA":
                    analyteSynonymOutFile.write(listItem.encode('utf-8') + b"\t" + 
                                                self.rampCompoundIDdictionary[key].encode('utf-8') + 
                                                b"\t" + b"compound" +b"\t"+
                                                database.encode("utf-8") +  b"\n")
        
        
        #Gene
        #gene with synonym
        for key in geneInfoDictionary:
            if key in self.rampGeneIDdictionary:
                
                mapping = geneInfoDictionary[key]
                commonName = mapping['common_name']
                if type(commonName) is not list:
                    commonName = commonName.replace("\n", "")
                    commonName = commonName.replace("\"", "")
                    commonName = commonName.replace(" ", "")
                    if commonName is not "NA":
                        analyteSynonymOutFile.write(commonName.encode('utf-8') + b"\t" + self.rampGeneIDdictionary[key].encode('utf-8') +
                                                     b"\t" + b"gene" +b"\t"+database.encode("utf-8")+ b"\n")
                else:
                    for item in commonName:
                        item = item.replace("\n", "")
                        item = item.replace("\"", "")
                        item = item.replace(" ", "")
                        if item is not "NA":
                            analyteSynonymOutFile.write(item.encode('utf-8') + b"\t" + 
                                                        self.rampGeneIDdictionary[key].encode('utf-8') 
                                                        + b"\t" + b"gene"+b"\t" +database.encode("utf-8") + b"\n")
                        
        #PATHWAY
        #pathway name   
        #rampid, resourcename, originaltype, originalid
        # Error
        for key in pathwayDictionary:
            try:
                pathwayOutFile.write(str(rampPathwayIDdictionary[key]).encode('utf-8') + b"\t"+ str(key).encode('utf-8') + b"\t"
                                                                            + str(database).encode('utf-8') + b"\t"  
                                                                            + str(pathwayCategory[key]).encode('utf-8') + b"\t"  
                                                                            + str(pathwayDictionary[key]).encode('utf-8') + b"\n")
                #print("Key Right !!")
            except KeyError:
                print('Key Error')
                print(key)
                print(pathwayDictionary[key])
                time.sleep(10)
        
        print("Metabolites linked to genes......................")
        for key in metabolitesLinkedToGenes:
            value = metabolitesLinkedToGenes[key]
            for listItem in value:
                try:
                    catalyzedOutFile.write(str(self.rampCompoundIDdictionary[key]).encode('utf-8') + b"\t" + str(self.rampGeneIDdictionary[listItem]).encode('utf-8') + b"\n")
                except:
                    pass 
       # Crosslink file         
        '''
        for key in metaboliteIDDictionary:
            mapping = metaboliteIDDictionary[key]
            
            #there are multiple chebi for every compound so you cannot 
            #just call mapping["chebi_id"]
            chebiJoined = "NA"
            if mapping["chebi_id"] is not "NA":
                chebiJoined = ":".join(mapping["chebi_id"])
                if chebiJoined == "":
                    chebiJoined = "NA"
            hmdbidJoined = "NA"
            if mapping["chebi_id"] is not "NA":
                hmdbidJoined = ":".join(mapping["hmdb_id"])
   
                if hmdbidJoined == "":
                    hmdbidJoined = "NA"
            
            compoundCrossLinksOutFile.write(self.rampCompoundIDdictionary[key].encode('utf-8') + b"\t" 
                                    + chebiJoined.encode('utf-8')  + b"\t"
                                    + mapping["knapsack_id"].encode('utf-8') + b"\t"
                                    + mapping["pubchem_compound_id"].encode('utf-8') + b"\t"
                                    + mapping["drugbank_id"].encode('utf-8')  + b"\t"
                                    + mapping["drugbank_metabolite_id"].encode('utf-8')  + b"\t"
                                    + mapping["phenol_explorer_compound_id"].encode('utf-8')  + b"\t"
                                    + mapping["phenol_explorer_metabolite_id"].encode('utf-8') + b"\t"
                                    + mapping["foodb_id"].encode('utf-8') + b"\t"
                                    + mapping["chemspider_id"].encode('utf-8') + b"\t"
                                    + mapping["kegg_id"].encode('utf-8') + b"\t"
                                    + mapping["biocyc_id"].encode('utf-8') + b"\t"
                                    + mapping["bigg_id"].encode('utf-8') + b"\t"
                                    + mapping["wikipidia"].encode('utf-8') + b"\t"
                                    + mapping["nugowiki"].encode('utf-8') + b"\t"
                                    + mapping["metagene"].encode('utf-8') + b"\t"
                                    + mapping["metlin_id"].encode('utf-8') + b"\t" 
                                    + mapping["het_id"].encode('utf-8') + b"\t"
                                    + hmdbidJoined.encode('utf-8') + b"\t"
                                    + mapping["CAS"].encode('utf-8') + b"\n")
        
        
        

        
        for key in geneInfoDictionary:
           
            mapping = geneInfoDictionary[key]
            
            ensemblJoined = "NA"
            if mapping["Ensembl"] is not "NA":
                ensemblJoined = ":".join(mapping["Ensembl"])
                if ensemblJoined == "":
                    ensemblJoined = "NA"
                    
            uniprotJoined = "NA"
            if mapping["UniProt"] is not "NA":
                uniprotJoined = ":".join(mapping["UniProt"])
                if uniprotJoined == "":
                    uniprotJoined = "NA"
            
            if key in self.rampGeneIDdictionary:
                if type(mapping["common_name"]) is not list:
                    geneCrossLinksOutFile.write(self.rampGeneIDdictionary[key].encode('utf-8') + b"\t"
                                                + mapping['common_name'].encode('utf-8') + b"\t" 
                                                + mapping['kegg'].encode('utf-8') + b"\t" 
                                                + ensemblJoined.encode('utf-8') + b"\t" 
                                                + mapping['HGNC'].encode('utf-8') + b"\t" 
                                                + mapping['HPRD'].encode('utf-8') + b"\t" 
                                                + mapping['NCBI-GeneID'].encode('utf-8') + b"\t"  
                                                + mapping['NCBI-ProteinID'].encode('utf-8') + b"\t"  
                                                + mapping['OMIM'].encode('utf-8') + b"\t" 
                                                + uniprotJoined.encode('utf-8') + b"\t" 
                                                + mapping['Vega'].encode('utf-8') + b"\t" 
                                                + mapping['miRBase'].encode('utf-8') + b"\t" 
                                                + mapping['HMDB_protein_accession'].encode('utf-8') + b"\t" 
                                                + str(mapping['Entrez']).encode('utf-8') + b"\t" 
                                                + mapping['Enzyme Nomenclature'].encode('utf-8')
                                                + b"\n")
                else:
                    for item in mapping["common_name"]:
                        geneCrossLinksOutFile.write(self.rampGeneIDdictionary[key].encode('utf-8') + b"\t"
                                                + item.encode('utf-8') + b"\t" 
                                                + mapping['kegg'].encode('utf-8') + b"\t" 
                                                + ensemblJoined.encode('utf-8') + b"\t" 
                                                + mapping['HGNC'].encode('utf-8') + b"\t" 
                                                + mapping['HPRD'].encode('utf-8') + b"\t" 
                                                + mapping['NCBI-GeneID'].encode('utf-8') + b"\t"  
                                                + mapping['NCBI-ProteinID'].encode('utf-8') + b"\t"  
                                                + mapping['OMIM'].encode('utf-8') + b"\t" 
                                                + uniprotJoined.encode('utf-8') + b"\t" 
                                                + mapping['Vega'].encode('utf-8') + b"\t" 
                                                + mapping['miRBase'].encode('utf-8') + b"\t" 
                                                + mapping['HMDB_protein_accession'].encode('utf-8') + b"\t" 
                                                + str(mapping['Entrez']).encode('utf-8') + b"\t" 
                                                + mapping['Enzyme Nomenclature'].encode('utf-8')
                                                + b"\n")
        '''
        # construct sql file that has rampOLId/CommonName/BiofluidORCellular
        # key is biofluid location (string) 
        # Total : biofluid, cellular, exo/endo , tissue
    
        for key in biofluidLocation:
            value = biofluidLocation[key]
            for listItem in value:
                analyteHasOntologyLocationOutFile.write(self.rampCompoundIDdictionary[key].encode('utf-8') + b"\t" 
                                                                                                    + rampBiofluidIDdictionary[listItem].encode('utf-8') + b"\n")
        
        for key in biofluid:
            ontologyLocationOutFile.write(rampBiofluidIDdictionary[key].encode('utf-8') + b"\t" 
                                                                           + key.encode('utf-8') + b"\t" 
                                                                    + b"biofluid" + b"\n")
             
        for key in cellularLocation:
            value = cellularLocation[key]
            for listItem in value:
                analyteHasOntologyLocationOutFile.write(self.rampCompoundIDdictionary[key].encode('utf-8') + b"\t" 
                                                                                                            + rampCellularIDdictionary[listItem].encode('utf-8') + b"\n")
      
        for key in cellular:
            ontologyLocationOutFile.write(rampCellularIDdictionary[key].encode('utf-8') + b"\t" + key.encode('utf-8') + b"\t" + b"cellular location" + b"\n")
        
        for key in pathwayOntology:
            listOfPathways =  pathwayOntology[key] 
            for listItem in listOfPathways:
                pathwayOntologyOutFile.write(key.encode('utf-8')  + b"\t" + listItem.encode('utf-8') + b"\n")
        for key in exoEndo:
            ontologyLocationOutFile.write(rampExoEndoIDdictionary[key].encode('utf-8') + b"\t" + key.encode('utf-8') + b"\t" + b"origins" + b"\n")  
        
        for key in exoEndoDictionary:
            listOfExoEndo = exoEndoDictionary[key]
            for item in listOfExoEndo:
                analyteHasOntologyLocationOutFile.write(self.rampCompoundIDdictionary[key].encode('utf-8') + b"\t" + rampExoEndoIDdictionary[item].encode('utf-8') + b"\n" )
        
        for key in tissue:
            ontologyLocationOutFile.write(rampTissueIDdictionary[key].encode('utf-8') + b"\t" + key.encode('utf-8') + b"\t" +b"tissue location" + b"\n")    
        for key in tissueLocation:
            listOfTissue = tissueLocation[key]
            for item in listOfTissue:
                analyteHasOntologyLocationOutFile.write(self.rampCompoundIDdictionary[key].encode('utf-8') + b"\t" + rampTissueIDdictionary[item].encode('utf-8') + b"\n")
        # Close all files ...
        analyteOutFile.close()
        analyteSynonymOutFile.close()
        analyteHasPathwayOutFile.close()
        pathwayOutFile .close()
        catalyzedOutFile.close()
        geneCrossLinksOutFile.close()
        compoundCrossLinksOutFile.close()
        sourceOutFile.close()
        ontologyLocationOutFile.close()
        analyteHasOntologyLocationOutFile.close()
        pathwayOntologyOutFile.close()
        endoExoOutFile.close()
        return finalRAMPIDnumbers
        
    def writeIdInWhichdatabase(self):
        rampCompIdToDbFile = open("../misc/sql/" + "compoundIdInWhichDB.sql", 'wb')
        for key in self.rampCompoundIdInWhichDatabases:
            dbset = self.rampCompoundIdInWhichDatabases[key]
            str = ""
            for db in dbset:
                if str =="":
                    str = str + db
                else:
                    str = str + ',' +db
            rampCompIdToDbFile.write(key.encode('utf-8') + b'\t' + str.encode('utf-8') +b'\n')

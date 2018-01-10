import mygene

class IDconversion():

    
    
    def MetaboliteChebiToHMDB(self, otherMetaboliteIDDictionary, hmdbMetaboliteIDDictionary, database):
        
        '''
        The purpose of this function is find HMDBids for the metabolites that do not have them, but could. This is done by looking through 
        all the compounds in a dictionary and looking to see if there is a chebi that are the same in both hmdb and the "other" database 
        (kegg, reactome, or wikipathways). Then, HMDBid can be given to the mapping for the "other" database. 
        
        Or, in peusdocode:
        
        if hmdbCompoundCHEBI== otherCompoundCHEBI:
            THEY ARE THE SAME COMPOUND
            give hmdbCompoundChebi to otherCompoundChebi
        
        The converted IDs will be places in a mapping that will be later written to sql files for the RaMP database, as well as written
        to an output file (a text file) listing the converted IDs. 
        
        :param dict otherMetaboliteIDDictionary: a dictionary of metabolites, obtained by using the keggData. reactomeData, or wikipathwaysData class 
        :param dict hmdbMetaboliteIDDictionary: a dictionary of hmdb metabolites, obtained by using the hmdbData class
        :param str database: a string to identify the output files containing the converted IDs
          
        '''
        
        #open outfile
        metaboliteIDConversionOutFile = open("../misc/output/" + str(database) + "MetaboliteIDConversionChebiToHMDB.txt", 'wb')
        
        for key in otherMetaboliteIDDictionary:
            othermapping = otherMetaboliteIDDictionary[key]
            otherchebiid = othermapping["chebi_id"]
            for otherchebi in otherchebiid:
                for hmdbkey in hmdbMetaboliteIDDictionary:
                    hmdbmapping = hmdbMetaboliteIDDictionary[hmdbkey]
                    hmdbchebiid = hmdbmapping["chebi_id"]
                    for hmdbchebi in hmdbchebiid:
                        if otherchebi == hmdbchebi:
                            othermapping["hmdb_id"] = hmdbmapping["hmdb_id"] 
                            
                            
                            ###formatting for writing to files##################
                            if othermapping["hmdb_id"] is not "NA":
                                hmdbidJoined = ":".join(othermapping["hmdb_id"])
   
                            if hmdbidJoined == "":
                                hmdbidJoined = "NA"
                            #####################################################
                            
                            
                            #write to file listing converted IDs
                            metaboliteIDConversionOutFile.write(hmdbidJoined.encode('utf-8') + b":" + hmdbchebi.encode('utf-8') + b"\n")
                        
                        
            
            otherMetaboliteIDDictionary[key] = othermapping    
            
    
    def MetaboliteKeggIDToChebi(self, keggMetaboliteIDDictionary, hmdbMetaboliteIDDictionary, database):
        
        '''
        The purpose of this function is find additional chebiIDs for the hmdb metabolites via KeggIDs (the keggID acts as a bridge). 
        This is done by looking through all the kegg metabolites and looking for keggIDs that are the same in both the kegg compound and in the hmdb compound. 
        Then, if the kegg compound also has a chebi, but the hmdb compound is missing the chebi, the chebi can be "imputed" this way.
        
        Or, in peusdocode:
        
        if hmdbCompoundKeggID == keggCompoundKeggID:
            THEY ARE THE SAME COMPOUND
            give keggCompoundChebi to hmdbCompoundChebi
        
        The converted IDs will be places in a mapping that will be later written to sql files for the RaMP database, as well as written
        to an output file (a text file) listing the converted IDs. 
        
        :param dict keggMetaboliteIDDictionary: a dictionary of kegg metabolites, obtained by using the keggData class 
        :param dict hmdbMetaboliteIDDictionary: a dictionary of hmdb metabolites, obtained by using the hmdbData class
        :param str database: a string to identify the output files containing the converted IDs
          
        '''
        #open outfile 
        metaboliteIDConversionOutFileKeggToChebi = open("../misc/output/" + str(database) + "MetaboliteIDConversionKeggtoChebi.txt", 'wb')
        
        for keggkey in keggMetaboliteIDDictionary:
            keggmapping = keggMetaboliteIDDictionary[keggkey]
            keggkeggid = keggmapping["kegg_id"]
            
            for hmdbkey in hmdbMetaboliteIDDictionary:
               hmdbmapping = hmdbMetaboliteIDDictionary[hmdbkey]
               hmdbkeggid = hmdbmapping["kegg_id"]
               if keggkeggid == hmdbkeggid and keggmapping["chebi_id"] is not "NA":
                   hmdbmapping["chebi_id"] = keggmapping["chebi_id"]
                   
                   #the found chebi is placed in mapping here
                   hmdbMetaboliteIDDictionary[hmdbkey] = hmdbmapping
                   
                   #write to file listing converted IDs
                   metaboliteIDConversionOutFileKeggToChebi.write(keggkeggid.encode('utf-8') + b":" + ":".join(keggmapping["chebi_id"]).encode('utf-8') + b"\n")
                   
            
            
                
            
    def GeneUniprotToHMDBP(self, otherGeneInfoDictionary, hmdbGeneInfoDictionary, database):
        
        '''
        The purpose of this function is to find HMDB identifiers for uniprot ids by looking to find a uniprot id match between
        the hmdb database and another database. If there is a match then the HMDB id (HMDB_protein_accession) is added to the 
        mapping for the other database. 
        
        This function uses a python module that converts gene IDs, found here: https://pypi.python.org/pypi/mygene
        
        :param dict otherGeneInfoDictionary: a dictionary of genes, obtained by using any database class except hmdb 
        :param dict hmdbGeneInfoDictionary: a dictionary of genes, obtained by using the hmdbData class
        :param str database: a string to identify the output files containing the converted IDs
        
        
        '''
        
        geneIDConversionOutFile = open("../misc/output/" + str(database) + "GeneIDConversion.txt", 'wb')
        
        
        for key in otherGeneInfoDictionary:
            othermapping = otherGeneInfoDictionary[key]
            otherUniprotidList = othermapping["UniProt"]
            for eachOtherUniprot in otherUniprotidList:
                if eachOtherUniprot is not "NA":
                    for hmdbkey in hmdbGeneInfoDictionary:
                        hmdbmapping = hmdbGeneInfoDictionary[hmdbkey]
                        hmdbUniprotList = hmdbmapping["UniProt"]
                        for eachUniprot in hmdbUniprotList:
                            if eachOtherUniprot == eachUniprot:
                                othermapping["HMDB_protein_accession"] = hmdbmapping["HMDB_protein_accession"] 
                                geneIDConversionOutFile.write(key.encode('utf-8') + b":" + hmdbmapping["HMDB_protein_accession"].encode('utf-8') + b"\n")  
            otherGeneInfoDictionary[key] = othermapping
    
    def GeneConvert(self, GeneInfoDictionaryToConvert, database):
        
        '''
        This function converts among entrez, uniprot, and ensembl geneids -- if there is one of the above, it will find the other two. 
        This will be helpful later on, when comparing genes in databases since they all use different IDs. 
        
        :param dict GeneInfoDictionaryToConvert: a dictionary of genes to convert
        :param str database: "hmdb", "reactome", etc.
        
        '''
        
        listOfGenesToQuery = []
        mapIDtoKey = {}
        
        for key in GeneInfoDictionaryToConvert:
            print('convert HMDBP ID ' + key )
            value = GeneInfoDictionaryToConvert[key]
            entrez = value["Entrez"]
            uniprot = value["UniProt"]
            ensembl = value["Ensembl"]
            
            if entrez is not "NA":
                listOfGenesToQuery.append(entrez)
                mapIDtoKey[entrez] = key
            
            if uniprot is not "NA":
                for eachuni in uniprot:
                    if eachuni is not "NA":
                        listOfGenesToQuery.append(eachuni)
                        mapIDtoKey[eachuni] = key
            if ensembl is not "NA":
                for eachensembl in ensembl:
                    if eachensembl is not "NA":
                        listOfGenesToQuery.append(eachensembl)
                        mapIDtoKey[eachensembl] = key
            
            
            
        mg = mygene.MyGeneInfo()
        queryResult = mg.querymany(listOfGenesToQuery, scopes='entrezgene,uniprot,ensembl.gene', fields='entrezgene,uniprot,ensembl.gene', species='human')
        for result in queryResult:
            original = result['query']
            notFoundTrue = False
            if "notfound" in result:
                notFoundTrue = result["notfound"]
            
            if not notFoundTrue:
                listOfUniprotID = []
                listOfEnsemblID = []
                entrezResult = entrez
                if "uniprot" in result:
                    if "Swiss-Prot" in result['uniprot']: 
                        swissProt = result['uniprot']['Swiss-Prot']
                        
                        if not isinstance(swissProt, list):
                            listOfUniprotID.append(swissProt)
                        else: 
                            
                            for each in swissProt:
                                listOfEnsemblID.append(each)
                if "ensembl" in result:
                    if len(result['ensembl']) == 1:
                        listOfEnsemblID.append(result['ensembl']['gene'])
                    else: 
                        for ensemblid in result['ensembl']:
                            listOfEnsemblID.append(ensemblid['gene'])
                if "entrezgene" in result:
                    entrezResult = result['entrezgene']
                    
                
                if len(listOfUniprotID) == 0:
                    listOfUniprotID = "NA"
                if len(listOfEnsemblID) == 0:
                    listOfEnsemblID = "NA"
                
                
                if original in listOfGenesToQuery:
                    key = mapIDtoKey[original]
                    mapping = GeneInfoDictionaryToConvert[key]
                    if mapping['Ensembl'] is "NA":
                        mapping['Ensembl'] = listOfEnsemblID
                    
                    if mapping['UniProt'] is "NA":
                        mapping['UniProt'] = listOfUniprotID
                    
                    if mapping['Entrez'] is "NA":
                        mapping['Entrez'] = entrezResult
                    GeneInfoDictionaryToConvert[key] = mapping
            
            
            
            
                
        
            
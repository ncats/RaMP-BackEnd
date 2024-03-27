'''
Created on Nov 24, 2020

@author: braistedjc
'''

import pandas as pd

class Gene(object):
    """
    The Gene class holds key fields used for defining a gene within ramp.
    This entity class was specifically created to hold a collection of ids, common names and synonyms.    
    """
    
    def __init__(self):
        '''
        Constructor
        '''
        # data source id for the gene
        self.sourceId = ""
        
        # gene ramp id
        self.rampId = ""
        
        # uniuqe list of ids
        self.idList = list()

        # source: id dictionary
        self.idDict = dict()
        
        # keys are source, values are dictionaries of source id to commonName
        self.commonNameDict = dict()
        
        # synonym dictionary, organized by source, id, synonym
        self.synonymDict = dict()
               
        # primary data source for the gene entity.        
        self.primarySource = ""

        # pathway objects organized by source, [source]:list(pathway)
        self.pathways = dict()
        
        # list of sources associated with the gene.          
        self.sources = list()
        
        self.proteinType = "Unknown"
        
        # status of review by uniprot.
        # TrEMBL uniprot entries are not fully curated (0)
        # SwissProt uniprot entries are completely curated (1)
        self.isReviewed = 0
    
    def __eq__(self, other):
        """
        Returns true if the genes either share ramp id, or if the associated ids have at least one overlap
        """
        if self.rampId and other.rampId:
            return self.rampId == other.rampId
        return len(set(self.idList).intersection(set(other.idList))) > 0
    
    
    def shareAltIds(self, other):
        """
        Returns true if two genes share an associaetd id
        """
        return len(set(self.idList).intersection(set(other.idList))) > 0
       
       
    def __hash__(self):
        """
        Returns a hashcode value used for evaluating equality in python collection methods.
        """
        if self.rampId:
            return hash(str(self.rampId))
        else:
            return hash("&".join(self.idList))
        
        
    def addId(self, id, source):
        """
        Adds an id to the general id list and the source-specific id lists
        """
        # keep this a unique id list
        if id not in self.idList:
            self.idList.append(id)
        if source not in self.idDict:
            self.idDict[source] = list()
        if id not in self.idDict[source]:
            self.idDict[source].append(id)
        
            
    def addSource(self, sourceName):
        """
        Adds a new source on the gene.
        """
        # maintain as a unique list
        if sourceName not in self.sources:
            self.sources.append(sourceName)
        
            
    def addCommonNameAndSynonym(self, id, commonName, source, annoType):
        """
        Adds common names as a <source | id | common_name> triple.
        """
        if annoType == 'common_name' or annoType == 'gene_name' or annoType == 'protein_name':
            if source not in self.commonNameDict:            
                self.commonNameDict[source] = dict()
            self.commonNameDict[source][id] = commonName
        self.addSynonym(commonName, source)
        
            
    def addPathway(self, pathway, source):
        """
        Adds a pathway record object and associated source to the gene.
        """
        if source not in self.pathways:
            self.pathways[source] = list()
        if pathway not in self.pathways[source]:
            self.pathways[source].append(pathway)
            
            
    def addSynonym(self, synonym, source):
        """
        Adds a gene synonym and related source to the gene
        """
        if source not in self.synonymDict:
            self.synonymDict[source] = list()
        if synonym not in self.synonymDict[source]:
            self.synonymDict[source].append(synonym)  
    
    
    def toSourceString(self):
        """
        Utility method to create a string suitable for source data output to file.
        The format is suitable for populating the source table.
        """
        lines = 0
        s = ""
        for source in self.commonNameDict:
            
            haveKeggPathwayMap = self.checkPathwaySourceLink(source, "kegg")
            
            for id in self.commonNameDict[source]:
                
                currId = id
                
                if isinstance(id, float) or isinstance(id, int):
                    print(self.printGene())
                
                idSplit = id.split(":")
                if len(idSplit) > 1:
                    idType = idSplit[0]
                else:
                    idType = "gene_symbol"                    
                    currId = "gene_symbol:" + id

                lines = lines + 1
                
                currSource = source
                                
                commonName = str(self.commonNameDict[source][id])
                commonName = commonName.replace("gene_symbol:", "", 1)
                
                if haveKeggPathwayMap:
                    if source == 'hmdb':
                        currSource = 'hmdb_kegg'             
                    if source == 'wiki':
                        currSource = 'wikipathways_kegg'
                    # add a row for current source, embedded kegg
                    s = s + str(currId) + "\t" + str(self.rampId) + "\t" + str(idType) + "\tgene\t" + commonName + "\t" + "N/A" + "\t" + str(currSource) + "\n"
                        
                             
                s = s + str(currId) + "\t" + str(self.rampId) + "\t" + str(idType) + "\tgene\t" + commonName + "\t" + "N/A" + "\t" + str(source) + "\n"
                
        return s
       
       
    def toPathwayMapString(self):
        """
        Utility method that returns a string for outputing pathway map information.
        """
        s = ""
        for source in self.pathways:
            for pathway in self.pathways[source]:
                s = s + self.rampId + "\t" + pathway.pathRampId + "\t" + str(pathway.pathSource) + "\n"
        return s
       
    def toSynonymsString(self):
        """
        Utility method to return a string for outputting synonyms information prior to db upload.
        """
        s = ""
        for source in self.synonymDict:
            for syn in self.synonymDict[source]:
                if(syn.startswith("gene_symbol:")):
                    syn = syn.split(":")[1]
                s = s + syn + "\t" + self.rampId + "\tgene\t" + source + "\n" 
        
        return s
    
    
    def printGene(self):
        """
        Utility to print a gene summary to standard output
        """
        s = "rampId: " + self.rampId + "\n"        
        for source in self.idDict:
            for id in self.idDict[source]:
                s = s + "id: " + str(id) + ", source: " + source + "\n"
        for source in self.commonNameDict:
            for id in self.commonNameDict[source].keys():
                s = s + "id: " + str(id) + ", source: " + str(source) + ", name:" + str(self.commonNameDict[source][id]) + "\n"
        
        for source in self.synonymDict:
            for syn in self.synonymDict[source]:
                s = s + "synonym: " + str(syn) + ", source: " + str(source) + "\n"
                
        print(s)
        
        print("pathway count: " + str(len(self.pathways)))
        
        for source in self.pathways.keys():
            for pathway in self.pathways[source]:
                pathway.printPathway()
            
             
    
    
    
#     The original code in writeToSQL would steal the commonName for a gene
#     from another common name when one instance of the gene (from the same source) lacked a common name.
#     It's a bit of a hack but for genes it will complete this field for the source table.
#     Also, all common names will have distinct rows with the same rampId.
#     This means that common name will resolve to a ramp id from any common name available.
#     
    def resolveCommonNames(self):
        for source in self.idDict:
            for id in self.idDict[source]:
                if source in self.commonNameDict and id not in self.commonNameDict[source]:
                    # now we know we have a common name dictionary for the source
                    # and our id doesn't have a commmon name entry.
                    
                    #grab a key
                    keyId = list(self.commonNameDict[source].keys())[0]
                    self.commonNameDict[source][id] = self.commonNameDict[source][keyId]
                
                else:
                    # what if we don't have a common name dictionary for the source, no common name
                    # just make an entry for the id but with no common name. This handles genes without common names
                    if source not in self.commonNameDict:
                        self.commonNameDict[source] = dict()
                        self.commonNameDict[source][id] = ""
        
    def subsumeGene(self, gene):
        """
        Copies data/objects from the passed gene into the invoking gene object.
        This method is used when merging gene entities.
        """
        # copy ids
        for source in gene.idDict:
            for id in gene.idDict[source]:
                self.addId(id, source)
        for source in gene.commonNameDict:
            for id in gene.commonNameDict[source]:
                self.addCommonName(id, gene.commonNameDict[source][id] ,source)
        for source in gene.sources:
            self.addSource(source)
        for source in gene.synonymDict:
            for syn in gene.synonymDict[source]:
                self.addSynonym(syn, source)        
            
            
        # check if a pathway of a given data source and category exists.        
    def checkPathwaySourceLink(self, source, category):
        jointMembership = False
        # check for pathway membership
        if source in self.pathways:
            for pathway in self.pathways[source]:                
                jointMembership = pathway.checkPathwaySourceAndCategory(source, category)
                if jointMembership:
                    return jointMembership
        return jointMembership
            
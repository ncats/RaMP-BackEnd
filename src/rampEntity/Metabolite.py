'''
Created on Nov 6, 2020

@author: braistedjc
'''
from statistics import median
import math

class Metabolite(object):
    '''
    Ramp metabolite object data container. A ramp metabolite represents one or more chemical entities representing
    a metabolite or metabolite class.
    '''
    
    # Class variable to indicate the equality policy
    # 0 = ID based, 1 = full lychi-based, 2 = lychi H3 based, 3 = full InchiKey based, 4 = InchiKey prefix match
    __metaboliteEqualityMetric = 0
    
    def __init__(self):
        
        self.sourceId = ""
        
        self.rampId = ""
        
        # uniuqe list of ids
        self.idList = list()

        # source: id dictionary
        self.idDict = dict()
        
        # keys are source, values are dictionaries of source id to commonName
        self.commonNameDict = dict()
        
        self.synonymDict = dict()
                
        self.primarySource = ""

        # source to pathway map
        self.pathways = dict()
     
        self.sources = list()
        
        self.associatedGenes = list()
        
        # source: [source_id: Molecule]  Molecule = molecule props
        self.chemPropsMolecules = dict()
        
        # layered dictionary source:dict(classLevel:list(className))
        self.metClasses = dict()
      
        self.ontologyTerms = list()
        
        self.hmdbStatus = None
        
        self.isCofactor = 0
                 
        self.inchiPrefixNeigbors = list()         
                 
    def __eq__(self, other):
        """
        Equality check. Checks for shared ramp_id or a non-empty intersection of ids lists.
        """ 
        if self.rampId and other.rampId:
            return self.rampId == other.rampId
        return len(set(self.idList).intersection(set(other.idList))) > 0
    
    
    def shareAltIds(self, other):
        """
        Check for shared ids
        """
        return len(set(self.idList).intersection(set(other.idList))) > 0
       
    def __hash__(self):
        """
        Hash code generator for ramp Metabolite objects
        """
        if self.rampId:
            return hash(str(self.rampId))
        else:
            return hash("&".join(self.idList))
                
    def addId(self, id, source):
        """
        Adds a new id if it doesn't exist
        The id is added to a the idList, a simple list, as well as an id dictionary that 
        groups ids by source.
        """
        # keep this a unique id list
        if id not in self.idList:
            self.idList.append(id)
        if source not in self.idDict:
            self.idDict[source] = list()
        if id not in self.idDict[source]:
            self.idDict[source].append(id)
    
                
    def printMet(self):
        """
        Utility method to print metabolite record to standard out
        """
        print("METABOLITE RECORD\n")
        s= "RampID: " + self.rampId + "\n"
        s= s + "sourceId " + self.sourceId + "\n"
        
        print(s)
        
        s = ""
        for src in self.idDict:
            for id in self.idDict[src]:
                s = s + "source: " + src + " altId: " + id + "\n"
        print(s)
        
        s=""
        for src in self.commonNameDict:
            for id in self.commonNameDict[src]:
                s = s + "source: " + src + " id: " + id + " commonName: " + self.commonNameDict[src][id] + "\n"
            
        print(s)
        
        s=""
        for src in self.synonymDict:
            for syn in self.synonymDict[src]:
                s = s + "source: " + src + " syn: " + syn + "\n"
        
        print("pathways")
        for source in self.pathways.keys():
            print(source + " pathways")
            for pathway in self.pathways[source]:
                pathway.printPathway()
            print("")
        
        
    def addSource(self, sourceName):
        """
        Adds source name to list of sources
        """
        # maintain as a unique list
        if sourceName not in self.sources:
            self.sources.append(sourceName)
    
    
    def getSortedSources(self):
        """
        Returns the list of sources associated with the metabolite
        """
        srcs = set(self.sources)
        list(srcs).sort()
        return srcs
    
    def mergeMets(self, otherMet):
        """
        Deprecated: Merges ID lists and sources for two metabolites.
        """
        for id in otherMet.idList:
            self.idList.append(id)
        for src in otherMet.sources:
            self.sources.append(src)
        
            
    def addPathway(self, pathway, source):
        """
        Adds a pathway
        Pathways are separated by source.
        """
        if source not in self.pathways:
            self.pathways[source] = list()
        if pathway not in self.pathways[source]:
            self.pathways[source].append(pathway)
            
         
    def getPathwayCount(self):
        """
        Returns the number of pathways to which the metabolite is associated 
        """
        count = 0
        for source in self.pathways:
            count = len(self.pathways[source])
        return count
            
    def addCommonName(self, id, commonName, source):
        """
        Adds a metabolite common names as a triple source, id, common name 
        """
        if source not in self.commonNameDict:            
            self.commonNameDict[source] = dict()
        self.commonNameDict[source][id] = commonName
        
    def addSynonym(self, synonym, source):
        """
        Adds a metoblite synonym, separated by source
        """
        if source not in self.synonymDict:
            self.synonymDict[source] = list()
        if synonym not in self.synonymDict[source]:
            self.synonymDict[source].append(synonym)
    
    def addAssociatedGene(self, gene):
        if gene not in self.associatedGenes:
            self.associatedGenes.append(gene)
            
            
    def subsumeMetabolite(self, metabolite):
        """
        Utility method to transfer contents from the passed Metabolite to the invoking Metabolite.
        This method is used to merge metabolite entities.
        """
        # copy ids
        for source in metabolite.idDict:
            for id in metabolite.idDict[source]:
                self.addId(id, source)
        for source in metabolite.commonNameDict:
            for id in metabolite.commonNameDict[source]:
                self.addCommonName(id, metabolite.commonNameDict[source][id] ,source)
        for source in metabolite.sources:
            self.addSource(source)
        for source in metabolite.synonymDict:
            for syn in metabolite.synonymDict[source]:
                self.addSynonym(syn, source)
    
        if metabolite.isCofactor == 1:
            self.isCofactor = 1
    
    def resolveCommonNames(self):
        """
        Verifies that all ids are associated with common names for export to the source table.
        1/12/2021 - This can be refined to be more robust, id -> source_id -> source_common name.
        """
        for source in self.idDict:
            for id in self.idDict[source]:
                if source not in self.commonNameDict:
                    self.commonNameDict[source] = dict()
                              
                if id not in self.commonNameDict[source]:
                    # now we know we have a common name dictionary for the source
                    # and our id doesn't have a commmon name entry.                    
                    # grab a key
                    if len(list(self.commonNameDict[source].keys())) < 1:
                        keyId = None
                    else:                            
                        keyId = list(self.commonNameDict[source].keys())[0]
                    
                    if keyId is None:
                        self.commonNameDict[source][id] = "NA"
                    else:
                        self.commonNameDict[source][id] = self.commonNameDict[source][keyId]
    
    def addChemProps(self, molecule):
        """
        Adds a molecule object to the metabolite.
        The molecule is added based on molecular information source and the molecule's id.
        Chemical properties are stored as source-id-molecule.` 
        """
        if molecule.source not in self.chemPropsMolecules:
           self.chemPropsMolecules[molecule.source] = dict()
        self.chemPropsMolecules[molecule.source][molecule.id] = molecule
        
        # we should pass on molecule names as synonyms
        for name in molecule.names:
            self.addSynonym(name, molecule.source)
    
    def addMetClass(self, source, sourceId, classLevel, className):
        if className == "Triradylcglycerols":
            className = "Triacylglycerols"
        
        sourceClasses = self.metClasses.get(source, None)
        if sourceClasses is None:
            self.metClasses[source] = dict()
            self.metClasses[source][sourceId] = dict()
            self.metClasses[source][sourceId][classLevel] = list()
            self.metClasses[source][sourceId][classLevel].append(className)
        else:
            sourceIdDict = self.metClasses[source].get(sourceId, None)
            
            if(sourceIdDict is None):
                self.metClasses[source][sourceId] = dict()
               
            classLevelList = self.metClasses[source][sourceId].get(classLevel, None)
            if classLevelList is None:
                self.metClasses[source][sourceId][classLevel] = list()
                self.metClasses[source][sourceId][classLevel].append(className)
            else:
                if className not in self.metClasses[source][sourceId][classLevel]:
                    self.metClasses[source][sourceId][classLevel].append(className)

    def addOntologyTerm(self, ontology):
        if ontology not in self.ontologyTerms:
            self.ontologyTerms.append(ontology)
    
    
    def toSourceString(self):
        """
        Utility method to return a tab delimited string suitable for source export.
        """
        lines = 0
        s = ""
        if self.hmdbStatus is not None:
            status = self.hmdbStatus
        else:
            status = "no_HMDB_status"
        for source in self.commonNameDict:
            for id in self.commonNameDict[source]:
                idSplit = id.split(":")
                if len(idSplit) > 1:
                    idType = idSplit[0]
                else:
                    idType = "None"
                lines = lines + 1
                currSource = source
                if idType == 'kegg':
                    if source == 'hmdb':
                        currSource = 'hmdb_kegg'
                    if source == 'wiki':
                        currSource = 'wikipathways_kegg'    
                s = s + str(id) + "\t" + str(self.rampId) + "\t" + str(idType) + "\tcompound\t" + str(self.commonNameDict[source][id]) + "\t" + status + "\t" + str(currSource) + "\n"
            #s = s.strip()

#         for source in self.idDict:
#             for id in self.idDict[source]:
#                 idSplit = id.split(":")
#                 if len(idSplit) > 1:
#                   idType = idSplit[0]
#                 else:
#                     idType = "NA"
#                 s = s + str(id) + "\t" + str(self.rampId) + "\tcompound\t" + str(idType) + "\t" + "CommonNamePlaceHolder" + "\t" + str(source) +"\n"'''
        return s;    

    
    def toPathwayMapString(self):
        """
        Utility method to return a tab delimited string for analyte to pathway export.
        """
        s = ""
        for source in self.pathways:
            for pathway in self.pathways[source]:
                s = s + self.rampId + "\t" + pathway.pathRampId + "\t" + str(pathway.pathSource) + "\n"
        return s
    

    def toSynonymsString(self):
        """
        Utility method to return a tab delimited string for analyte to pathway export.
        """
        s = ""
        for source in self.synonymDict:
            for syn in self.synonymDict[source]:                    
                s = s + str(syn) + "\t" + self.rampId + "\tcompound\t" + source + "\n" 
        
        return s
    
    def toMetaboliteOntologyString(self):
        s = ""
        for ontology in self.ontologyTerms:
            s = s + self.rampId + "\t" + ontology.ontolRampId + "\n"
        return s

    def toChemPropsString(self):
        """
        Utility method to return a chemical property string suitable for exporting chemical properties.
        """
        s = ""
        for source in self.chemPropsMolecules:
            for id in self.chemPropsMolecules[source]:
                mol = self.chemPropsMolecules[source][id]                
                chemProps = mol.toChemPropsString()
                if chemProps is not None:
                    s = s + self.rampId + "\t" + chemProps
        return s
    
    def toMetToGeneAssociationString(self):
        s = ""
        
        for gene in self.associatedGenes:
            s = s + self.rampId + "\t" + gene.rampId + "\t" + gene.proteinType + "\n"
        
        return s
    
    def toMetaboliteClassString(self):
        s = ""
        for src in self.metClasses:
            for metSourceId in self.metClasses[src]:
                for classLevel in self.metClasses[src][metSourceId]:          
                    for className in self.metClasses[src][metSourceId][classLevel]:
                        s = ( s + str(self.rampId) + 
                        "\t" + str(metSourceId) + 
                        "\t" + str(classLevel) + 
                        "\t" + str(className) + 
                        "\t" + str(src) + "\n")
        return s                
    
    
    def checkMWParity(self, mwTolerance, pctOrAbs):
        """
        Returns 0.0 if the molecular weights between contained molecules is within the mass tolerance
        based on pct (fraction) or absolute mw difference.
        """
        mwList = list() 
        
        for source in self.chemPropsMolecules:
            for id in self.chemPropsMolecules[source]:
                mol = self.chemPropsMolecules[source][id]
                if mol.mw is not None and len(mol.mw) > 0:
                    mwList.append(float(mol.mw))
        
        if len(mwList) < 2:
            return 0.0
        
        if pctOrAbs == "pct":
            pct = (max(mwList) - min(mwList))/min(mwList)
            if pct > mwTolerance:
                return pct
            else:
                return 0.0
        
        if (max(mwList) - min(mwList)) > mwTolerance:
            return (max(mwList) - min(mwList))
        else:
            return 0
    
    
    def checkInchiBaseParity(self):
        """
        Returns the number of distinct inchikey base values.
        If the return is 1, then all molecules have the same constituents and connectivity. 
        """
        inchiDict = dict() 
        
        for source in self.chemPropsMolecules:
            for id in self.chemPropsMolecules[source]:
                mol = self.chemPropsMolecules[source][id]
                if mol.inchiKey is not None and len(mol.inchiKey) > 10:
                    inchiBase = mol.inchiKey.split("-")[0]
                    inchiDict[inchiBase] = inchiBase
        
        return len(inchiDict)
    
    
    def toCommonNameJoinString(self):
        """
        Utility method to supply a concatenated common name with semicolon separator. 
        """
        cname = list()
        for source in self.commonNameDict:
            for id in self.commonNameDict[source]:
                cname.append(self.commonNameDict[source][id])
                
        return "; ".join(cname)        
            
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
    
            
    #def setStatus(self):            
    def getInchiPrefixes(self):
        inchiPrefixes = []
        for source in self.chemPropsMolecules:
            molDict = self.chemPropsMolecules[source]
            for sourceId in molDict:
                mol = molDict[sourceId]
                if mol.inchiKeyPrefix is not "" and mol.inchiKeyPrefix not in inchiPrefixes:
                    inchiPrefixes.append(mol.inchiKeyPrefix)
        return inchiPrefixes
    
    def getInchiKeys(self):
        inchiKeys = []
        for source in self.chemPropsMolecules:
            molDict = self.chemPropsMolecules[source]
            for sourceId in molDict:
                mol = molDict[sourceId]
                if mol.inchiKey is not "" and mol.inchiKey not in inchiKeys:
                    inchiKeys.append(mol.inchiKey)
        return inchiKeys
    
    def getInchiKeyDuplexes(self):
        inchiKeyDuplexes = []
        for source in self.chemPropsMolecules:
            molDict = self.chemPropsMolecules[source]
            for sourceId in molDict:
                mol = molDict[sourceId]
                if mol.inchiKeyDuplex is not "" and mol.inchiKeyDuplex not in inchiKeyDuplexes:
                    inchiKeyDuplexes.append(mol.inchiKeyDuplex)
        return inchiKeyDuplexes
        
    def addInchiNeighbor(self, otherMet):
        if self is not otherMet:
            if otherMet not in self.inchiPrefixNeigbors:
                # if the other met hasn't already been added as a neighbor
                for neighbor in self.inchiPrefixNeigbors:
                    # add neighbors to the other met - introductions...
                    neighbor.addInchiNeighbor(otherMet)
                    # introduce other neighbor to existing neighbors
                    otherMet.addInchiNeighbor(neighbor)
                
                # finally add the new neighbor to the neighbor list    
                self.inchiPrefixNeigbors.append(otherMet)
            
    
    def getInchiNeighborhood(self):
        
        neighbors = self.inchiPrefixNeigbors
        
        for neighbor in self.inchiPrefixNeigbors:
            neighbor.getNeighbors(neighbors)

        return neighbors

    # recursive get neighbors
    def getNeighbors(self, neighbors):
               
        for neighbor in self.inchiPrefixNeigbors:
            # just work on new neighbors, add the neighbor and get their neighbors
            if(neighbor not in neighbors):
                neighbors.append(neighbor)
                neighbor.getNeighbors(neighbors)

     
    def getAveMW(self):
        mws = []
        medMw = 0.0
        for source in self.chemPropsMolecules:
            molDict = self.chemPropsMolecules[source]
            for sourceId in molDict:
                mol = molDict[sourceId]
                
               if(mol.mw is not None and mol.mw != ""):
                    mw = float(mol.mw)
                    if not math.isnan(mw):
                        mws.append(mw)
                                        
        if(len(mws) > 0):
            medMw = median(mws)
               
        return medMw
    
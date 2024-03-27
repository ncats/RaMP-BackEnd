'''
Created on Jun 30, 2022

@author: braistedjc
'''

class RheaReaction(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        
        self.rxnRampId = ""
        
        # Reaction metadata
        self.rhea_id = ""
        
        self.parent_rhea_id = ""
        
        self.rhea_label = ""
        
        self.rhea_equation = ""
        
        self.rhea_html_eq = ""
        
        #UN, LR, RL, BD
        self.direction = ""

        self.status = 1
        
        # 1:n associations
        
        # mapping rhea_id to directional id records
        self.rhea_directional_ids = []
        
        # should the related directional records just link to other reaction entities
        self.rhea_directional_rxns = []
        
        # possible mapping... are ids to ecs one to many?
        self.ec = None
        
        # mapping rhea id to uniprot ids, can be 1:n
        self.proteins = []
        
        self.proteinDict = dict()
        
        self.isTransport = 0        
                
        # compounds ids for left side, populate using rhea2ec tsv file.
        self.left_comp_ids = []
        
        # compound ids for right side
        self.right_comp_ids = []
        
        # compounds ids for left side, populate using rhea2ec tsv file.
        self.left_comps = []
        
        # compound ids for right side
        self.right_comps = []
        
        self.hasHumanEnzyme = False
        
        self.hasOnlyHumanMetabolites = True
        
        self.hasAHumanMetabolite = False
                
        self.ecAssociationBlock = []     
                
                
    def getBasicRecordString(self):
        ec = self.ec
        
        if ec is None:
            ecVal = ""
        else:
            eCount = 0
            ec = sorted(ec)
            for e in ec:
                if eCount == 0:
                    ecVal = e
                else:
                    ecVal = ecVal + "; " + e
                        
        dir = self.direction
        if dir is None:
            dir = ""
        
        humanEnzyme = self.hasHumanEnzyme * 1
        onlyHumanMets = self.hasOnlyHumanMetabolites * 1
        
        s = (self.rhea_id + "\t" + str(self.status) + "\t" + str(self.isTransport) + "\t" +self.direction + "\t" + self.rhea_label + "\t" + 
             self.rhea_equation + "\t" + self.rhea_html_eq + "\t" + ecVal + "\t" + str(humanEnzyme) + "\t" + str(onlyHumanMets) +"\n")
       
        return s
    
    
    def getMainRecordString(self):
        ec = self.ec
        
        if ec is None or type(ec) == float:
            ecVal = ""        
        else:
            if len(ec) == 1:
                ecVal = ec[0]
            else:
                eCount = 0
                ec = sorted(ec)
                for e in ec:
                    if eCount == 0:
                        ecVal = e
                    else:
                        ecVal = ecVal + "; " + e
            
        direction = self.direction
        
        if direction is None:
            direction = ""
        
        humanEnzyme = self.hasHumanEnzyme * 1
        onlyHumanMets = self.hasOnlyHumanMetabolites * 1
        
        s = str(self.rxnRampId) + "\t" + str(self.rhea_id) + "\t" + str(self.status) + "\t" + str(self.isTransport) + "\t"
        s = s + str(direction) + "\t" + str(self.rhea_label) + "\t" 
        s = s + str(self.rhea_equation) + "\t" + str(self.rhea_html_eq) + "\t" + str(ecVal) + "\t" + str(humanEnzyme) + "\t" + str(onlyHumanMets) + "\n"
       
        return s
    
    
    def getMainReactionToMetString(self, source):
        
        s = ""
        
        for c in self.left_comps:
            ids = c.idDict.get(source, None)
            if ids is not None and len(ids) > 0:
                metId = ids[0]
            else:
                metId = ""
                
            namesDict = c.commonNameDict.get(source, None)
            name = ""
            if namesDict is not None:
                cid = list(namesDict.keys())[0]
                name = namesDict[cid]
                         
            s = s + self.rxnRampId + "\t" + self.rhea_id + "\t" + c.rampId + "\t0\t" + metId + "\t" + name + "\t" + str(c.isCofactor) + "\n"
    
        for c in self.right_comps:
            ids = c.idDict.get(source, None)
            if ids is not None and len(ids) > 0:
                metId = ids[0]
            else:
                metId = ""
                
            namesDict = c.commonNameDict.get(source, None)
            name = ""
            if namesDict is not None:
                cid = list(namesDict.keys())[0]
                name = namesDict[cid]
            
            s = s + self.rxnRampId + "\t" + self.rhea_id + "\t" + c.rampId + "\t1\t" + metId + "\t" + name + "\t" + str(c.isCofactor) + "\n"
            
        return s    

    
    def getMainReactionToProteinString(self, source):
        
        s = ""
                
        for p in self.proteins:

            uniprot = p.sourceId
#            ids = p.idDict.get(source, None)
#            if ids is not None and len(ids) > 0:                
#                uniprot = ids[0]
#            else:
#                uniprot = 'NA'            

            names = p.commonNameDict.get(source, None)
            name = " "
            
            if names is not None:
                name = names.get(uniprot, None)
                if name is None:
                    name = "UNK"
                    #print("export rxn to prot, HAVE NAME DICT, BUT NO NAME for uniprot: "+uniprot + " DICT LEN: " + str(len(list(names.keys()))))
                else:
                    # print("Have a name in rxt to prot... but it's an empty string **|"+name+"|**")
                    if name == "":
                        name = "UNK5"
                        
                        #print("Dumping Names...")
                        #for i in names:
                        #    print(str(i) + "---" + str(names[i]))
 
            else:
                name = "UNK2"
                #print("export rxn to prot, NO NAME DICT")
            
            if name == "":
                name = "UNK3"
            
            if name == " ":
                name = "UNK4"
                
            s = s + self.rxnRampId + "\t" + self.rhea_id + "\t" + p.rampId + "\t" + uniprot + "\t" + name + "\t" + str(p.isReviewed) + "\n"
        
        return s    
           
           
           
    def getMainReactionToProteinStringAllIds(self, source, secondaryIdSet, swissProtIds, tremblProtIds):
        
        s = ""
        
        if self.direction != 'UN':
            return s
        
        uniprotIdSet = set()
        
        outputCount = 0
        
        if swissProtIds is None:
            swissCount = 0
        else:
            swissCount = len(swissProtIds)

        if tremblProtIds is None:
            tremblCount = 0
        else:
            tremblCount = len(tremblProtIds)

        
        totalIdCount = swissCount + tremblCount
        
        for p in self.proteins:

            ids = p.idDict.get(source, None)
            
            if ids is not None and len(ids) > 0:
                
                for uniprot in ids:
                    
                    if uniprot in secondaryIdSet:                        
                        continue
                    
                    if swissProtIds is not None and uniprot in swissProtIds:
                        isReviewed = 1
                        outputCount = outputCount + 1
                    elif tremblProtIds is not None and uniprot in tremblProtIds:
                        isReviewed = 0
                        outputCount = outputCount + 1
                    else:
                        # limits reporting of non-rhea proteins
                        continue
                    
                    #uniprot = ids[0]                    
                          
                    names = p.commonNameDict.get(source, None)
                    name = " "
                    
                    if names is not None:
                        name = names.get(uniprot, None)
                        if name is None:
                            name = "UNK"
                            #print("export rxn to prot, HAVE NAME DICT, BUT NO NAME for uniprot: "+uniprot + " DICT LEN: " + str(len(list(names.keys()))))
                        else:
                            # print("Have a name in rxt to prot... but it's an empty string **|"+name+"|**")
                            if name == "":
                                name = ""
                                
                                #print("Dumping Names...")
                                #for i in names:
                                #    print(str(i) + "---" + str(names[i]))
         
                    else:
                        name = ""
                        #print("export rxn to prot, NO NAME DICT")
                    
#                     if name == "":
#                         name = "UNK3"
#                     
#                     if name == " ":
#                         name = "UNK4"

                    if uniprot not in uniprotIdSet:                        
                        s = s + self.rxnRampId + "\t" + self.rhea_id + "\t" + p.rampId + "\t" + uniprot + "\t" + name + "\t" + str(isReviewed) + "\n"
                        uniprotIdSet.add(uniprot)
        
        if outputCount < totalIdCount:
            print("rhea protein output mismatch. totalIdCount = " + str(totalIdCount), + " outputCount = " + str(outputCount))
        
        
        return s    
    
           
           
    def getReactionProteinToMetString(self, source):
        
        s = ""
        
        hitProteins = []
        for p in self.proteins:
            if p.rampId not in hitProteins:
                hitProteins.append(p.rampId)
                
                ids = p.idDict.get(source, None)
                if ids is not None and len(ids) > 0:
                    uniprot = ids[0]
                else:
                    uniprot = 'NA'
                    
                hitMets = []
                for met in self.left_comps:
                    if met.rampId not in hitMets:
                        hitMets.append(met.rampId)
                        
                        namesDict = met.commonNameDict.get(source, None)
                        name = ""
                        if namesDict is not None:
                            cid = list(namesDict.keys())[0]
                            name = namesDict[cid]
                
                        s = s + self.rxnRampId + "\t" + self.rhea_id + "\t" + p.rampId + "\t" + uniprot + "\t0\t" + met.rampId + "\t" + cid + "\t" + name + "\t" + str(met.isCofactor) + str(p.isReviewed) + "\n"
                        
                for met in self.right_comps:
                    if met.rampId not in hitMets:
                        hitMets.append(met.rampId)
                        
                        namesDict = met.commonNameDict.get(source, None)
                        name = ""
                        if namesDict is not None:
                            cid = list(namesDict.keys())[0]
                            name = namesDict[cid]
                
                        s = s + self.rxnRampId + "\t" + self.rhea_id + "\t" + p.rampId + "\t" + uniprot + "\t1\t" + met.rampId + "\t" + cid + "\t" + name + "\t" + str(met.isCofactor) + "\n"
                        
        return s                
                
    
    def assignPrimaryFields(self, dataVals):
        self.rhea_id = dataVals[0]
        self.status = dataVals[1]
        self.isTransport = dataVals[2]
        self.direction = dataVals[3]
        self.rhea_label = dataVals[4]
        self.rhea_equation = dataVals[5]
        self.rhea_html_eq = dataVals[6]
        self.ec = [dataVals[7]]
        self.hasHumanEnzyme = dataVals[8]
        self.hasOnlyHumanMetabolites = dataVals[9]
        
    
    def getCompById(self, cid, compList):
        #print("search id" + str(cid))
        for cmp in compList:
            #print("search id: " + str(cid))
            #print("checkCompList have a comp: " + str(cmp.chebiId))
            if cmp.chebiId == cid:
                #print("match")
                return cmp
        return None
    
    # exports reaction to chebi id reaction side and is cofactor
    def getRheaIdToCompMappingString(self):
        s = ""
        for cid in self.left_comp_ids:
            cmpd = self.getCompById(cid, self.left_comps)
            isCofactor = 0
            if cmpd is not None:
                isCofactor = cmpd.isCofactor
                
                s = s + self.rhea_id + "\t" + cid + "\t0\t" + str(isCofactor) + "\n" 
       
        for cid in self.right_comp_ids:
            cmpd = self.getCompById(cid, self.right_comps)
            isCofactor = 0
            if cmpd is not None:
                isCofactor = cmpd.isCofactor
                
                s = s + self.rhea_id + "\t" + cid + "\t1\t" + str(isCofactor) + "\n"            

        return s
    
    def getRheaIdToUniprotMappingString(self):    
        s = ""
        for pid in self.proteins:
            pName = ""
            p = self.proteinDict.get(pid, None)
            protIsReviewed = 0
            if p is not None:
                protIsReviewed = p.isReviewed
            else:
                protIsReviewed = 3
                
            s = s + self.rhea_id + "\t" + pid + "\t0\t" + str(protIsReviewed) + "\n"             
        return s

    def getCompoundToProteinString(self):
        s = ""
        for pid in self.proteins:
            for cid in self.left_comp_ids:
                s = s + self.rhea_id + "\t" + pid + "\t" + cid + "\t0\n"
            for cid in self.right_comp_ids:
                s = s + self.rhea_id + "\t" + pid + "\t" + cid + "\t1\n"
        return s
                
                
    def getRheaReactionToEcString(self):
        ecBlock = ""
        if len(self.ecAssociationBlock) > 0:
            self.ecAssociationBlock = list(set(self.ecAssociationBlock))
            for ecData in self.ecAssociationBlock:
                ecBlock = ecBlock + self.rxnRampId + "\t" + self.rhea_id + "\t" + ecData.strip() + "\n"
        
        return ecBlock
    
    def addEcAssociationBlock(self, ecData):
        self.ecAssociationBlock.append(ecData) 
        
   
    def doIHaveACofactorCompoundCheck(self):
        for c in self.left_comps:
            if c.isCofactor == 1:
                return 1

        for c in self.right_comps:
            if c.isCofactor == 1:
                return 1
        
        return 0    
            
#     def doIHaveACofactorCompoundIdCheck(self):
#         for cid in self.left_comp_ids:
#             c = self.getCompById(cid, self.left_comps)
#             if c is not None:
#                 if c.isCofactor == 1:
#                     return 1
#             else:
#                 print("have an left ID with no compound???? "+ str(cid))
# 
#         for cid in self.right_comp_ids:
#             c = self.getCompById(cid, self.right_comps)
#             if c is not None:
#                 if c.isCofactor == 1:
#                     return 1
#             else:
#                 print("have an right ID with no compound???? "+ str(cid))
# 
#         return 0   
            
    def addProteinSet(self, accList, accToProteinMap):
        for acc in accList:
            p = accToProteinMap.get(acc, None)
            if p is not None:
                self.proteinDict[acc] = p
                
        

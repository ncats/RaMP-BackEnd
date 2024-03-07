'''
Created on Jun 30, 2022

@author: braistedjc
'''
import sys, os
from os.path import exists
import gzip
import shutil
from parse.MetabolomicsData import MetabolomicsData
from rampEntity.Protein import Protein
from rampConfig.RampConfig import RampConfig
import numpy as np

from rdflib import URIRef, Graph
import rdflib.namespace
from builtins import str
from rampEntity.RheaReaction import RheaReaction
from rampEntity.RheaCompound import RheaCompound
from parse.UniprotParser import UniprotParser
from parse.ChebiOwlParser import ChebiOwlParser

import pandas as pd



class RheaParser(MetabolomicsData):
    '''
    classdocs
    '''


    def __init__(self, resConfig):
        '''
        Constructor
        '''
        # relative dir, use '../' for testing, use "" for production calls
        self.relDir = ""     
        
        self.config = resConfig

        self.graph = None
        
        self.rheaReactionDict = dict()
        
        self.rheaCompoundDict = dict()
        
        self.rheaProteinDict = dict()
        
        self.rheaEcToClassDict = dict()
        
        self.rheaLocalRdfFile = ""
        
        self.rheaLocalRheaToUniprotFile = ""
        
        self.rheaLocalRheaToSwissprotFile = ""
        
        self.rheaLocalRheaToEcFile = ""

        self.rheaLocalRxnDirectionFile = ""
        
        self.expasyLocalEc2ClassFile = ""
        
        self.humanUniprotRecordDict = dict()
        
        self.humanUniprotAccSet = set()
        
        self.chebiHumanIdSet = set()
        
        self.chebiCofactorId = "23357"
        
        self.chebiCofactorSet = set()
                     
    def processRhea(self):
        self.buildSupportingUniprotData()
        self.buildSupportingChebiData()
        
        print("Length of cofactor set = "+str(len(self.chebiCofactorSet)))
                
        self.getRheaFiles()
        self.constructRDF()
                 
        # builds reactions objects
        self.processAllReactions()
        
        # this gets expasy ec to enzyme class
        self.ecToEnzymeClassFromExpasy()        
        
        self.appendUniprotToReaction()
        self.appendEcToReaction()
         
         
        self.setReactionHumanUniprotState()
        self.setReactionHumanChebiState()
         
        self.exportIntermediateFiles()
    


    def buildSupportingUniprotData(self):
        print("building uniprot data store")
        uniParser = UniprotParser(self.config)
        
        # this will parse human uniprot, including SwissProt and TrEMBL accessions.
        uniParser.parseUniprot()
        self.humanUniprotRecordDict = uniParser.uniprotRecords
        
        threadedUniprotDict = dict()
        
        for acc in self.humanUniprotRecordDict:
            self.humanUniprotAccSet.add(acc)
            p = self.humanUniprotRecordDict[acc]
            threadedUniprotDict[acc] = p
            for acc2 in p.secondaryAccs:
                self.humanUniprotAccSet.add(acc2)
                threadedUniprotDict[acc2] = p
                
        self.humanUniprotRecordDict = threadedUniprotDict        
    
        print("in rhea uniprot build")
        print("length of the dict"+str(len(self.humanUniprotRecordDict.keys())))
        print("length of the set"+str(len(self.humanUniprotAccSet)))


    def buildSupportingChebiData(self):
        print("building chebi data store")
        
        cop = ChebiOwlParser(self.config)   
        cop.getChebiFiles()
        
        # use this when testing if a stored graph file exists... for faster graph construction
        # cop.deserializeGraph(None)
        # else.... use buildGraph()
        
        cop.buildGraph()                
        #cop.deserializeGraph(None)        
        
        cop.extractHumanMetaboliteStatus()
        
        self.chebiHumanIdSet = cop.chebiHumanIdSet
        
        self.chebiCofactorSet = cop.extractCofactorStatus(self.chebiCofactorId)
    
    def getRheaFiles(self):

        rdfConf = self.config.getConfig('rhea_rdf')
        uniprotToRheaConf = self.config.getConfig('uniprot_to_rhea')
        swissprotToRheaConf = self.config.getConfig('swissprot_to_rhea')
        rheaToEcConf = self.config.getConfig('rhea_to_ec')
        rheaDirectionConf = self.config.getConfig('rhea_rxn_direction')
        expasyEc2EnzymeClassConf = self.config.getConfig('expasy_ec2class')

        localDir = rdfConf.localDir
        
        if not exists(self.relDir + localDir):
            os.mkdir(self.relDir + localDir)

        # first RDF file
        rdfLocalFile = rdfConf.extractFileName
                
        if not exists(self.relDir + localDir + rdfLocalFile):
            # get the file
            rdfUrl = rdfConf.sourceURL
            rdfRemoteFile = rdfConf.sourceFileName
            
            self.download_files(rdfUrl, self.relDir + localDir + rdfRemoteFile)            
            
            # now gunzip
            with gzip.open(self.relDir + localDir + rdfRemoteFile, 'rb') as f_in:
                with open(self.relDir + localDir + rdfLocalFile, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)       
        else:
            print("Using cached Rhea RDF file.")
        
        
        # rhea to uniprot
        rhea2UniprotFile = uniprotToRheaConf.extractFileName
        
        self.rheaLocalRheaToUniprotFile = self.relDir + localDir + rhea2UniprotFile
        
        if not exists(self.relDir + localDir + rhea2UniprotFile):
            rhea2UniprotUrl = uniprotToRheaConf.sourceURL
            rhea2UniprotRemoteFile = uniprotToRheaConf.sourceFileName
            
            self.download_files(rhea2UniprotUrl, self.relDir + localDir + rhea2UniprotRemoteFile)            
        else:
            print("Using cached Rhea Uniprot-to-Rhea file.")
        
        
        
        # rhea to swissprot
        rhea2SwissprotFile = swissprotToRheaConf.extractFileName
        
        self.rheaLocalRheaToSwissprotFile = self.relDir + localDir + rhea2SwissprotFile
        
        if not exists(self.relDir + localDir + rhea2SwissprotFile):
            rhea2SwissprotUrl = swissprotToRheaConf.sourceURL
            rhea2SwissprotRemoteFile = swissprotToRheaConf.sourceFileName
            
            self.download_files(rhea2SwissprotUrl, self.relDir + localDir + rhea2SwissprotRemoteFile)            
        else:
            print("Using cached Rhea swissprot-to-Rhea file.")
            
        
        
        
        # rhea to ec
        rhea2EcFile = rheaToEcConf.extractFileName
        
        self.rheaLocalRheaToEcFile = self.relDir + localDir + rhea2EcFile
        
        if not exists(self.relDir + localDir + rhea2EcFile):
            rhea2EcUrl = rheaToEcConf.sourceURL
            rhea2EcRemoteFile = rheaToEcConf.sourceFileName
            
            self.download_files(rhea2EcUrl, self.relDir + localDir + rhea2EcRemoteFile)            
        else:
            print("Using cached Rhea Uniprot-to-EC file.")

        #rhea_rxn_direction
        rheaReactionDirectionFile = rheaDirectionConf.extractFileName
        
        self.rheaLocalRxnDirectionFile = self.relDir + localDir + rheaReactionDirectionFile
        
        if not exists(self.relDir + localDir + rheaReactionDirectionFile):
            rheaDirUrl = rheaDirectionConf.sourceURL
            rheaDirRemoteFile = rheaDirectionConf.sourceFileName
            
            self.download_files(rheaDirUrl, self.relDir + localDir + rheaDirRemoteFile)            
        else:
            print("Using cached Rhea reaction direction file.")


        # supporting expasy EC to Enzyme Class file
        ec2classFile = expasyEc2EnzymeClassConf.extractFileName

        self.expasyLocalEc2ClassFile = self.relDir + localDir + ec2classFile

        if not exists(self.relDir + localDir + ec2classFile):
            rheaDirUrl = expasyEc2EnzymeClassConf.sourceURL
            rheaDirRemoteFile = expasyEc2EnzymeClassConf.sourceFileName
            
            self.download_files(rheaDirUrl, self.relDir + localDir + rheaDirRemoteFile)            
        else:
            print("Using cached Expasy ec2enzymeClass file.")
            

     
    def constructRDF(self):
                
        g = Graph()
        rdfConf = self.config.getConfig('rhea_rdf')        
        
        path = self.relDir + rdfConf.localDir + rdfConf.extractFileName
        
        g.parse(path, format="application/rdf+xml")        
        
        self.graph = g        
       
        
        
    def processAllReactions(self):
        
        type_predicate = URIRef('http://www.w3.org/2000/01/rdf-schema#subClassOf')    

        reaction_object = URIRef('http://rdf.rhea-db.org/Reaction')
        dir_reaction_object = URIRef("http://rdf.rhea-db.org/DirectionalReaction")
        bidir_reaction_object = URIRef("http://rdf.rhea-db.org/BidirectionalReaction")

        res = self.graph.subjects(predicate=type_predicate, object=reaction_object)
        dirRes = self.graph.subjects(predicate=type_predicate, object=dir_reaction_object)
        biDirRes = self.graph.subjects(predicate=type_predicate, object=bidir_reaction_object)
        
        # process reactions that are non-directional 
        self.processReactions(self.graph, res)
        # process reactions that are directional
        self.processReactions(self.graph, dirRes)
        # process reactions that are bi-directional
        self.processReactions(self.graph, biDirRes)
        # now add directional info to reactions
        self.processReactionDirectionInfo()
        
    
    def processReactions(self, g, res):    
                
        acc_predicate = URIRef("http://rdf.rhea-db.org/accession")
        label_predicate = URIRef("http://www.w3.org/2000/01/rdf-schema#label")        
        eq_predicate = URIRef("http://rdf.rhea-db.org/equation")
        html_eq_predicate = URIRef("http://rdf.rhea-db.org/htmlEquation")
        part_pred = URIRef("http://rdf.rhea-db.org/contains")
        compound_pred = URIRef("http://rdf.rhea-db.org/compound")
        status_predicate = URIRef("http://rdf.rhea-db.org/status")
        
        name_pred = URIRef("http://rdf.rhea-db.org/name")
        html_name_pred = URIRef("http://rdf.rhea-db.org/htmlName")
        formula_pred = URIRef("http://rdf.rhea-db.org/formula")
        is_transport_pred = URIRef("http://rdf.rhea-db.org/isTransport")
        
        for rid in res:
            s = URIRef(rid)
            
            rxnAcc = g.objects(subject=s, predicate=acc_predicate)
            rxnLabel = g.objects(subject=s, predicate=label_predicate)            
            rxnEq = g.objects(subject=s, predicate=eq_predicate)
            rxnHtmlEq = g.objects(subject=s, predicate=html_eq_predicate)
            status = g.objects(subject=s, predicate=status_predicate)
            isTransport = g.objects(subject=s, predicate=is_transport_pred)
                        
            reaction = RheaReaction()
            
            # accession
            for acc in rxnAcc:
                reaction.rhea_id = acc.strip().lower()
                # append the reaction
                self.rheaReactionDict[reaction.rhea_id] = reaction
            
            for label in rxnLabel:                   
                reaction.rhea_label = label.strip()
            
            # equation
            for eq in rxnEq:                   
                reaction.rhea_equation = eq.strip()
             
            for htmlEq in rxnHtmlEq:    
                reaction.rhea_html_eq = htmlEq.strip()
        
            for s in status:
                statusVal = s.split("/")[-1]
                stat = statusVal.strip()
                if stat == 'Obsolete':
                    reaction.status = -1
                elif stat == 'Preliminary':
                    reaction.status = 0    
                
            for s in isTransport:
                # this is a rdflib.term.Literal, not a string
                # conversion is needed to use == below for equality on 'true'
                s = str(s)

                if s == "true":
                    reaction.isTransport = 1
                    
            # now we need to query for participants
            left_subj = URIRef(rid + "_L")
            right_subj = URIRef(rid + "_R")
        
#             print(left_subj)
#             print(right_subj)
#             print(part_pred)
        
            leftRefs = g.objects(predicate=part_pred, subject=left_subj)
            rightRefs = g.objects(predicate=part_pred, subject=right_subj)            
        
            for leftRef in leftRefs:
                
                leftSubj = URIRef(leftRef)
 
                compRefs = g.objects(predicate=compound_pred, subject = leftSubj)
                
                for compRef in compRefs:
                    
                    compound = RheaCompound()

                    comps = g.objects(predicate = acc_predicate, subject = URIRef(compRef))
                    compNames = g.objects(predicate = name_pred, subject = URIRef(compRef))
                    htmlNames = g.objects(predicate = html_name_pred, subject = URIRef(compRef))
                    formulas = g.objects(predicate = formula_pred, subject = URIRef(compRef))
                                        
                    for comp in comps:
                        cmp = comp.strip().lower()
                        cmpPrefix = cmp.split(":")[0]
                        if cmpPrefix == 'generic':
                            cmp = cmp.replace('generic', 'rhea-comp')
                        reaction.left_comp_ids.append(cmp)
                        compound.chebiId = cmp
                        if compound.chebiId in self.chebiCofactorSet:
                            compound.isCofactor = 1
                    
                    for name in compNames:
                        compound.name = name
                        
                    for htmlName in htmlNames:
                        compound.htmlName = htmlName
                    
                    for formula in formulas:
                        compound.formula = formula

                    self.rheaCompoundDict[compound.chebiId] = compound
                    reaction.left_comps.append(compound)

            for rightRef in rightRefs:
                
                rightSubj = URIRef(rightRef)
                
                compRefs = g.objects(predicate=compound_pred, subject = rightSubj)
                    
                for compRef in compRefs:
                                        
                    compound = RheaCompound()
                    
                    comps = g.objects(predicate = acc_predicate, subject = URIRef(compRef))
                    compNames = g.objects(predicate = name_pred, subject = URIRef(compRef))
                    htmlNames = g.objects(predicate = html_name_pred, subject = URIRef(compRef))
                    formulas = g.objects(predicate = formula_pred, subject = URIRef(compRef))
                    
                    for comp in comps:
                        cmp = comp.strip().lower()
                        cmpPrefix = cmp.split(":")[0]
                        if cmpPrefix == 'generic':
                            cmp = cmp.replace('generic', 'rhea-comp')
                        reaction.right_comp_ids.append(cmp)
                        compound.chebiId = cmp
                        if compound.chebiId in self.chebiCofactorSet:
                            compound.isCofactor = 1
                    
                    for name in compNames:
                        compound.name = name
                        
                    for htmlName in htmlNames:
                        compound.htmlName = htmlName
                    
                    for formula in formulas:
                        compound.formula = formula    
                              
                    self.rheaCompoundDict[compound.chebiId] = compound
                    reaction.right_comps.append(compound)
        
            
            
    def processReactionDirectionInfo(self):
        
        dirTable = pd.read_csv(self.rheaLocalRxnDirectionFile, sep="\t", header=0)
        
        dirMapping = dict()
        
        for idx, row in dirTable.iterrows():
            dirMapping['rhea:'+str(row[0]).strip()] = "UN"
            dirMapping['rhea:'+str(row[1]).strip()] = "LR"
            dirMapping['rhea:'+str(row[2]).strip()] = "RL"
            dirMapping['rhea:'+str(row[3]).strip()] = "BD"
            
            # only the UN reactions have compound id lists
            # others must inherit from the UN.

            rxn = self.rheaReactionDict.get('rhea:'+str(row[0]).strip(), None)
            lrRxn = self.rheaReactionDict.get('rhea:'+str(row[1]).strip(), None)
            rlRxn = self.rheaReactionDict.get('rhea:'+str(row[2]).strip(), None)
            bdRxn = self.rheaReactionDict.get('rhea:'+str(row[3]).strip(), None)
            
            if rxn is not None:
                if lrRxn is not None:
                    lrRxn.left_comp_ids = rxn.left_comp_ids
                    lrRxn.right_comp_ids = rxn.right_comp_ids
                if rlRxn is not None:
                    # note that when rxn is R to L, then the 
                    # formula still goes left to right when written, so relative side of compounds changes.
                    rlRxn.left_comp_ids = rxn.right_comp_ids
                    rlRxn.right_comp_ids = rxn.left_comp_ids
                if bdRxn is not None:
                    bdRxn.left_comp_ids = rxn.left_comp_ids
                    bdRxn.right_comp_ids = rxn.right_comp_ids
            
        
        
        print("In direction processing dirMapping size and reaction dict size")
        print(str(len(dirMapping)))
        print(str(len(self.rheaReactionDict)))
        noDir = 0
        obsol = 0
        other = 0
        noLeft = 0
        noRight = 0
        haveBoth = 0
        hasHumanUniprot = 0
        for rxnId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rxnId]
            dir = dirMapping.get(rxnId, None)
            
            if rxn.status == -1:
                obsol = obsol + 1
        
            if rxn.status == 0:
                other = other + 1
                
            if dir is not None:
                rxn.direction = dir
            else:
                noDir = noDir + 1
        
            if len(rxn.left_comp_ids) == 0:
                noLeft = noLeft + 1

            if len(rxn.right_comp_ids) == 0:
                noRight = noRight + 1
        
            if len(rxn.left_comp_ids) > 0 and len(rxn.right_comp_ids) > 0:
                haveBoth = haveBoth + 1
        
        print("no dir count and obsolete, other count")
        print(str(noDir))
        print(str(obsol))
        print(str(other))
        print("no Left, no right or have both compound status")
        print(str(noLeft))
        print(str(noRight))
        print(str(haveBoth))
        print("")
    
    def setReactionHumanUniprotState(self):
        print("setting uniprot human status")
        
        numHumanUniprot = 0
        numUNHumanUniprot = 0
        for rheaId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rheaId]
            for p in rxn.proteins:
                if p in self.humanUniprotAccSet:
                    if rxn.status == 1:
                        numHumanUniprot = numHumanUniprot + 1
                        if rxn.direction == 'UN':
                            numUNHumanUniprot = numUNHumanUniprot + 1
                            
                    rxn.hasHumanEnzyme = True
                    break
        
        print("number of rxn (status == 1) with human uniprot, and UN uniprot rxn count")
        print(str(numHumanUniprot))
        print(str(numUNHumanUniprot))
        
    def setReactionHumanChebiState(self):
        print("setting chebi human status")
        print("human chebi size = " + str(len(self.chebiHumanIdSet)))
        
        for rheaId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rheaId]
            for c in rxn.left_comp_ids:
                if c not in self.chebiHumanIdSet:
                    rxn.hasOnlyHumanMetabolites = False
                    break
            for c in rxn. right_comp_ids:
                if c not in self.chebiHumanIdSet:
                    rxn.hasOnlyHumanMetabolites = False
                    break    

        for rheaId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rheaId]
            for c in rxn.left_comp_ids:
                if c in self.chebiHumanIdSet:
                    rxn.hasAHumanMetabolite = True
                    break
            for c in rxn. right_comp_ids:
                if c in self.chebiHumanIdSet:
                    rxn.hasAHumanMetabolite = True
                    break

        # human reaction check
        humanRxnCnt = 0
        hasAHumanMet = 0
        for rheaId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rheaId]
            if rxn.hasOnlyHumanMetabolites and rxn.status == 1:
                humanRxnCnt = humanRxnCnt + 1
            if rxn.hasAHumanMetabolite and rxn.status == 1:
                hasAHumanMet = hasAHumanMet + 1
                
        humanUNRxnCnt = 0
        for rheaId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rheaId]
            if rxn.hasOnlyHumanMetabolites and rxn.direction == 'UN':
                humanUNRxnCnt = humanUNRxnCnt + 1
          
        humanLRRxnCnt = 0
        for rheaId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rheaId]
            if rxn.hasOnlyHumanMetabolites and rxn.direction == 'LR':
                humanLRRxnCnt = humanLRRxnCnt + 1

        humanRLRxnCnt = 0
        for rheaId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rheaId]
            if rxn.hasOnlyHumanMetabolites and rxn.direction == 'RL':
                humanRLRxnCnt = humanRLRxnCnt + 1

        humanBDRxnCnt = 0
        for rheaId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rheaId]
            if rxn.hasOnlyHumanMetabolites and rxn.direction == 'BD':
                humanBDRxnCnt = humanBDRxnCnt + 1
                                
        print("chebi human reaction count = " + str(humanRxnCnt))
        print("chebi UN human reaction count = " + str(humanUNRxnCnt))
        print("chebi LR human reaction count = " + str(humanLRRxnCnt))
        print("chebi RL human reaction count = " + str(humanRLRxnCnt))
        print("chebi BD human reaction count = " + str(humanBDRxnCnt))
        print("rxn HAS A (at least one) human chebi: " + str(hasAHumanMet))
          
    def exportIntermediateFiles(self):

        dir = self.relDir + "../misc/output/rhea_reactions/"
        
        
        if not exists(dir):
            os.mkdir(dir)
        
        recordsFile = "rhea_primary_records.txt"
        
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for acc in self.rheaReactionDict:
            reaction = self.rheaReactionDict[acc]
            recordOut.write(reaction.getBasicRecordString())

        recordOut.close()
     
     
        recordsFile = "rhea_rxn_to_chebi_and_dir.txt"
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for acc in self.rheaReactionDict:
            reaction = self.rheaReactionDict[acc]
            recordOut.write(reaction.getRheaIdToCompMappingString())

        recordOut.close()
        
        
        recordsFile = "rhea_compound_info.txt"
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for acc in self.rheaCompoundDict:
            compound = self.rheaCompoundDict[acc]
            recordOut.write(compound.rheaCompoundToRecordString())
            
        recordOut.close()
        
   
        recordsFile = "rhea_uniprot_mapping.txt"
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for acc in self.rheaReactionDict:
            rxn = self.rheaReactionDict[acc]
            recordOut.write(rxn.getRheaIdToUniprotMappingString())    
        
        recordOut.close()


        recordsFile = "rhea_reaction_to_ec.txt"
                
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for acc in self.rheaReactionDict:
            rxn = self.rheaReactionDict[acc]
            ecList = rxn.ec
            if ecList is not None and len(ecList) > 0:
                ecBlock = self.buildRxnEcExportBlock(acc, ecList)
                if len(ecBlock) > 0:
                    recordOut.write(ecBlock)
        
        recordOut.close()


        recordsFile = "rhea_compound_to_protein_mapping.txt"
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for acc in self.rheaReactionDict:
            rxn = self.rheaReactionDict[acc]
            recordOut.write(rxn.getCompoundToProteinString())    
                    
        recordOut.close()


        recordsFile = "rheametaboliteIDDictionary.txt"
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for cid in self.rheaCompoundDict:
            type = cid.split(":")[0]     
            recordOut.write(cid + "\t" + type + "\t" + cid + "\n")  
                    
        recordOut.close()


        recordsFile = 'rheametaboliteCommonName.txt'
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for cid in self.rheaCompoundDict:
            comp = self.rheaCompoundDict[cid]    
            recordOut.write(cid + "\t" + comp.name + "\n") 
                    
        recordOut.close()


        recordsFile = 'rheageneInfoDictionary.txt'
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for rxnId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rxnId]

            for pid in rxn.proteins:
                p = self.humanUniprotRecordDict.get(pid, None)
                if p is not None:
                    if p.recName is not None and p.recName != "":
                        recordOut.write(pid + "\tprotein_name\t" + p.recName + "\n")
                    if p.hgncSymbol is not None and p.hgncSymbol != "":
                        recordOut.write(pid + "\tgene_name\t" + p.hgncSymbol + "\n")
                    for secondaryAcc in p.secondaryAccs:         
                        recordOut.write(pid + "\tuniprot\t" + secondaryAcc + "\n")
                else:
                    print("No protein for id: "+ pid)

                
        recordOut.close()



        
    def appendUniprotToReaction(self):
        #self.rheaLocalRheaToEcFile
        #self.rheaLocalRheaToUniprotFile
        r2u = pd.read_csv(self.rheaLocalRheaToUniprotFile, sep="\t", header=0)
        
        r2uMap = dict()
        
        print(str(r2u.shape))
        
        for idx, row in r2u.iterrows():
            #print(row)
            #print("appending protein accessions to reactions..." + str(row.RHEA_ID)+ "  " +str(row.ID))

            # !!! just adding human uniprot            
            if ("uniprot:" + row.ID) in self.humanUniprotAccSet:
                #print("Have the human id!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                unis = r2uMap.get("rhea:" + str(row.RHEA_ID))
                if unis is None:     
                    unis = ['uniprot:'+row.ID]
                    r2uMap['rhea:'+str(row.RHEA_ID)] = unis
                else:
                    unis.append('uniprot:'+row.ID)



        for rxn in r2uMap:
            #print('adding uniprot')
            #print('reaction '+rxn)
            
            uniSet = r2uMap[rxn]
            currRxn = self.rheaReactionDict.get(rxn, None)
            if currRxn is None:
                currRxn = self.rheaReactionDict.get("rhea:"+rxn, None)
            if currRxn is not None:   
                currRxn.proteins = uniSet
                #print("setting proteins, len:"+str(len(currRxn.proteins)))
        
        
            
        # swiss prot    
        r2u = pd.read_csv(self.rheaLocalRheaToSwissprotFile, sep="\t", header=0)
        
        print(str(r2u.shape))
        
        for idx, row in r2u.iterrows():
            #print(row)
            #print("appending protein accessions to reactions..." + str(row.RHEA_ID)+ "  " +str(row.ID))

            # !!! just adding human uniprot            
            if ("uniprot:" + row.ID) in self.humanUniprotAccSet:
                #print("Have the human id!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                unis = r2uMap.get("rhea:" + str(row.RHEA_ID))
                if unis is None:     
                    unis = ['uniprot:'+row.ID]
                    r2uMap['rhea:'+str(row.RHEA_ID)] = unis
                else:
                    unis.append('uniprot:'+row.ID)

        for rxn in r2uMap:
            #print('adding uniprot')
            #print('reaction '+rxn)
            
            uniSet = r2uMap[rxn]
            currRxn = self.rheaReactionDict.get(rxn, None)
            if currRxn is None:
                currRxn = self.rheaReactionDict.get("rhea:"+rxn, None)
            if currRxn is not None:   
                currRxn.proteins = uniSet
                #print("setting proteins, len:"+str(len(currRxn.proteins)))
                
        
        
    def appendEcToReaction(self):
        #self.rheaLocalRheaToEcFile
        #self.rheaLocalRheaToUniprotFile
        r2u = pd.read_csv(self.rheaLocalRheaToEcFile, sep="\t", header=0)
        
        r2EcMap = dict()
        
        print(str(r2u.shape))
        
        for idx, row in r2u.iterrows():
            rheaRxnId = 'rhea:'+str(row.RHEA_ID)
            ecList = r2EcMap.get(rheaRxnId,None)
            if ecList is None:
                r2EcMap[rheaRxnId] = [row.ID]
            else:
                ecList.append(row.ID)

        for rxn in r2EcMap:
            ecList = r2EcMap[rxn]
            currRxn = self.rheaReactionDict.get(rxn, None)
            if currRxn is not None:
                currRxn.ec = list(set(ecList))
        
    def ecToEnzymeClassFromExpasy(self):
        
        # ec2class = pd.read_csv(self.expasyLocalEc2ClassFile, sep="\t", skiprows=11, skipfooter=5)
        with open(self.expasyLocalEc2ClassFile, 'r') as ec2c:
            ec2classStrings = ec2c.readlines()
        
            # The file has a header and a footer that have to be skipped.
            start = 11
            end = len(ec2classStrings) - 5
        
        for i in range(start, end):
            line = ec2classStrings[i].strip()
            ec_data = line.split("  ")
            ec = ec_data[0]
            enzClass = ec_data[1]
            if len(ec_data) == 3:
                enzClass = ec_data[2]
            ec = ec.replace(" ", "")
            enzClass = enzClass.strip() 
            self.rheaEcToClassDict[ec] = enzClass



    def buildRxnEcExportBlock(self, rxnId, ecList):
        ecBlock = ""
        enzClassJoin = ""
        for ec in ecList:
            ecChildren = self.getEcChildren(ec)
            enzClassJoin = ""
            i = 0
            for ecc in ecChildren:
                enzClass = self.rheaEcToClassDict.get(ecc, None)
                if enzClass is not None:
                    if i == 0:
                        enzClassJoin = enzClass
                        # just mark that we are past the first entry
                        i = 1
                    else:
                        # concatentate the enzyme class info :), I think this is finally correct :)
                        enzClassJoin = enzClassJoin + " | " + enzClass
                    ecLevel = 4 - ecc.count("-")
                    ecBlock = ecBlock + rxnId + "\t" + ecc + "\t" + str(ecLevel) + "\t" + enzClass + "\t" + enzClassJoin + "\n"
                    
                    
        return ecBlock             


    def getEcChildren(self, ec):
        data = ec.split('.')
        ecVariants = [ec]
        ecVariants.append(data[0] + "." + data[1] + "." + data[2] + ".-")
        ecVariants.append(data[0] + "." + data[1] + ".-.-")
        ecVariants.append(data[0] + ".-.-.-")
        ecVariants = sorted(ecVariants)
        return ecVariants

    
   
        
                
#rConf = RampConfig()
#rConf.loadConfig("../../config/external_resource_config.txt")
# # #                         
#rp = RheaParser(rConf)            
#rp.processRhea()
# rp.appendUniprotToReaction()
# rp.appendEcToReaction()
# rp.exportIntermediateFiles()

# rxn = rp.rheaReactionDict["RHEA:10000"]
# print(rxn.rhea_id)
# print("left")
# for leftPart in rxn.left_comps:
#     print(leftPart)
# print("right")
# for rightPart in rxn.right_comps:
#     print(rightPart)

# acc2 = ""
# d2 = dict()
# for acc in rp.rheaReactionDict:
#     #print("**"+acc+"**")
#     d2[acc] = rp.rheaReactionDict[acc]
#     acc2 = acc
#     
# print(str(len(rp.rheaReactionDict)))
#     
# 
# 
# rxn = d2[acc2]
# print(rxn.rhea_id)
# print("left")
# for leftPart in rxn.left_comp_ids:
#     print(leftPart)
# print("right")
# for rightPart in rxn.right_comp_ids:
#     print(rightPart)
# 
# 
# rxn = d2['rhea:31411']
# print(rxn.rhea_id)
# print("left")
# for leftPart in rxn.left_comp_ids:
#     print(leftPart)
# print("right")
# for rightPart in rxn.right_comp_ids:
#     print(rightPart)


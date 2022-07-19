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

from rdflib import URIRef,Graph
import rdflib.namespace
from rdflib.namespace import RDF,FOAF,RDFS,DC,DCTERMS
from builtins import str
from rampEntity.RheaReaction import RheaReaction
from rampEntity.RheaCompound import RheaCompound
from docutils.nodes import compound
from libchebipy._parsers import __FORMULAE

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
        self.relDir = "../"        
        
        self.config = resConfig
        
        self.rheaReactionDict = dict()
        
        self.rheaCompoundDict = dict()
        
        self.rheaLocalRdfFile = ""
        
        self.rheaLocalRheaToUniprotFile = ""
        
        self.rheaLocalRheaToEcFile = ""

        self.rheaLocalRxnDirectionFile = ""
        
        self.distinctStatus = set()
             
    def parseRheaFiles(self):
        self.getRheaFiles()
        self.constructRDF()
        self.exportIntermediatFiles()
    
    def getRheaFiles(self):

        rdfConf = self.config.getConfig('rhea_rdf')
        uniprotToRheaConf = self.config.getConfig('uniprot_to_rhea')
        rheaToEcConf = self.config.getConfig('rhea_to_ec')
        rheaDirectionConf = self.config.getConfig('rhea_rxn_direction')

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
            print("Using cacheed Rhea Uniprot-to-Rhea file.")

            
        # rhea to ec
        rhea2EcFile = rheaToEcConf.extractFileName
        
        self.rheaLocalRheaToEcFile = self.relDir + localDir + rhea2EcFile
        
        if not exists(self.relDir + localDir + rhea2EcFile):
            rhea2EcUrl = rheaToEcConf.sourceURL
            rhea2EcRemoteFile = rheaToEcConf.sourceFileName
            
            self.download_files(rhea2EcUrl, self.relDir + localDir + rhea2EcRemoteFile)            
        else:
            print("Using cacheed Rhea Uniprot-to-EC file.")

        #rhea_rxn_direction
        rheaReactionDirectionFile = rheaDirectionConf.extractFileName
        
        self.rheaLocalRxnDirectionFile = self.relDir + localDir + rheaReactionDirectionFile
        
        if not exists(self.relDir + localDir + rheaReactionDirectionFile):
            rheaDirUrl = rheaDirectionConf.sourceURL
            rheaDirRemoteFile = rheaDirectionConf.sourceFileName
            
            self.download_files(rheaDirUrl, self.relDir + localDir + rheaDirRemoteFile)            
        else:
            print("Using cacheed Rhea Uniprot-to-EC file.")


     
    def constructRDF(self):
        
        
        g = Graph()
        rdfConf = self.config.getConfig('rhea_rdf')
        
        relDir = "../"
        
        path = relDir + rdfConf.localDir + rdfConf.extractFileName
        
        g.parse(path, format="application/rdf+xml")        
        
        #g.parse(path, format="nt")        
       
        
        
        reactions = g.objects("Reaction")

        i = 0
        
        subs = []
        preds = []
        objs = []
#         for subj, pred, obj in g:            
#             i = i + 1
#             subs.append(subj)
#             print(obj)
#             preds.append(pred)
#             objs.append(obj)
#         
#         
#         print(i)
#         print(len(set(subs)))
#         print(len(set(preds)))
#         print(len(set(objs)))
        
        #res = g.subject_objects(predicate="rdfs:subClassOf")
    
        
#         q = """
#     PREFIX foaf: <http://xmlns.com/foaf/0.1/>
# 
#     SELECT ?name
#     WHERE {
#         ?p rdfs:subClassOf rdf:resource "http://rdf.rhea-db.org/Reaction"/ .
# 
#         ?p foaf:name ?name .
#     }
# """
#         q = """
#         PREFIX rh: <http://rdf.rhea-db.org/>
#         PREFIX rdfs:  <http://www.w3.org/2000/01/rdf-schema#>
#         PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
#         
#         SelectQuery ?rh:accession where { ?rdfs:subClassOf  rdf:resource: "http://rdf.rhea-db.org/Reaction .}
#         """
#         res = g.query(q)
#         
#         i = 0
#         for s in res:
#             i = i + 1
#             
#         print("number of reaction subject_object tuples")
#         print(str(i))
        
        type_predicate = URIRef('http://www.w3.org/2000/01/rdf-schema#subClassOf')    

        reaction_object = URIRef('http://rdf.rhea-db.org/Reaction')
        dir_reaction_object = URIRef("http://rdf.rhea-db.org/DirectionalReaction")
        bidir_reaction_object = URIRef("http://rdf.rhea-db.org/BidirectionalReaction")

        res = g.subjects(predicate=type_predicate, object=reaction_object)
        dirRes = g.subjects(predicate=type_predicate, object=dir_reaction_object)
        biDirRes = g.subjects(predicate=type_predicate, object=bidir_reaction_object)
        
        # process reactions that are non-directional 
        self.processReactions(g, res)
        # process reactions that are directional
        self.processReactions(g, dirRes)
        # process reactions that are bi-directional
        self.processReactions(g, biDirRes)
        # now add directional info to reactions
        self.processReactionDirectionInfo()
    
        for status in self.distinctStatus:
            print(status)
    
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
       
        for rid in res:
            s = URIRef(rid)
            
            rxnAcc = g.objects(subject=s, predicate=acc_predicate)
            rxnLabel = g.objects(subject=s, predicate=label_predicate)            
            rxnEq = g.objects(subject=s, predicate=eq_predicate)
            rxnHtmlEq = g.objects(subject=s, predicate=html_eq_predicate)
            status = g.objects(subject=s, predicate=status_predicate)

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
                reaction.status = statusVal.strip()
                self.distinctStatus.add(statusVal)
            
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
                    
                    for name in compNames:
                        compound.name = name
                        
                    for htmlName in htmlNames:
                        compound.htmlName = htmlName
                    
                    for formula in formulas:
                        compound.formula = formula

                    self.rheaCompoundDict[compound.chebiId] = compound
                    reaction.left_comps[compound.chebiId] = compound

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
                        reaction.left_comp_ids.append(cmp)
                        compound.chebiId = cmp
                    
                    for name in compNames:
                        compound.name = name
                        
                    for htmlName in htmlNames:
                        compound.htmlName = htmlName
                    
                    for formula in formulas:
                        compound.formula = formula    
                              
                    self.rheaCompoundDict[compound.chebiId] = compound
                    reaction.right_comps[compound.chebiId] = compound
        
            
            
    def processReactionDirectionInfo(self):
        
        dirTable = pd.read_csv(self.rheaLocalRxnDirectionFile, sep="\t", header=0)
        
        dirMapping = dict()
        
        for idx, row in dirTable.iterrows():
            dirMapping['rhea:'+str(row[0]).strip()] = "UN"
            dirMapping['rhea:'+str(row[1]).strip()] = "LR"
            dirMapping['rhea:'+str(row[2]).strip()] = "RL"
            dirMapping['rhea:'+str(row[3]).strip()] = "BD"
        
        for rxnId in self.rheaReactionDict:
            rxn = self.rheaReactionDict[rxnId]
            dir = dirMapping.get(rxnId, None)
            if dir is not None:
                rxn.direction = dir

    def exportIntermediatFiles(self):

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

        recordsFile = "rhea_compound_to_protein_mapping.txt"
        
        recordOut = open(dir + recordsFile, 'w', encoding="utf-8")
        for acc in self.rheaReactionDict:
            rxn = self.rheaReactionDict[acc]
            recordOut.write(rxn.getCompoundToProteinString())    
                    
        recordOut.close()

        
        
    def appendUniprotToReaction(self):
        #self.rheaLocalRheaToEcFile
        #self.rheaLocalRheaToUniprotFile
        r2u = pd.read_csv(self.rheaLocalRheaToUniprotFile, sep="\t", header=0)
        
        r2uMap = dict()
        
        print(str(r2u.shape))
        
        for idx, row in r2u.iterrows():
            #print(row)
            #print(str(row.RHEA_ID)+ "  " +str(row.ID))
            unis = r2uMap.get("rhea:" + str(row.RHEA_ID))
            if unis is None:
                unis = [row.ID]
                r2uMap['rhea:'+str(row.RHEA_ID)] = unis
            else:
                unis.append(row.ID)

        for rxn in r2uMap:
            uniSet = r2uMap[rxn]
            currRxn = self.rheaReactionDict.get(rxn, None)
            if currRxn is not None:   
                currRxn.proteins = uniSet
        
        
    def appendEcToReaction(self):
        #self.rheaLocalRheaToEcFile
        #self.rheaLocalRheaToUniprotFile
        r2u = pd.read_csv(self.rheaLocalRheaToEcFile, sep="\t", header=0)
        
        r2EcMap = dict()
        
        print(str(r2u.shape))
        
        for idx, row in r2u.iterrows():
            r2EcMap['rhea:'+str(row.RHEA_ID)] = row.ID

        for rxn in r2EcMap:
            ec = r2EcMap[rxn]
            currRxn = self.rheaReactionDict.get(rxn, None)
            if currRxn is not None:
                currRxn.ec = ec

        
rConf = RampConfig()
rConf.loadConfig("../../config/external_resource_config.txt")
                        
rp = RheaParser(rConf)            
rp.parseRheaFiles()
rp.appendUniprotToReaction()
rp.appendEcToReaction()
rp.exportIntermediatFiles()

# rxn = rp.rheaReactionDict["RHEA:10000"]
# print(rxn.rhea_id)
# print("left")
# for leftPart in rxn.left_comps:
#     print(leftPart)
# print("right")
# for rightPart in rxn.right_comps:
#     print(rightPart)

acc2 = ""
d2 = dict()
for acc in rp.rheaReactionDict:
    #print("**"+acc+"**")
    d2[acc] = rp.rheaReactionDict[acc]
    acc2 = acc
    
print(str(len(rp.rheaReactionDict)))
    


rxn = d2[acc2]
print(rxn.rhea_id)
print("left")
for leftPart in rxn.left_comp_ids:
    print(leftPart)
print("right")
for rightPart in rxn.right_comp_ids:
    print(rightPart)


rxn = d2['rhea:31411']
print(rxn.rhea_id)
print("left")
for leftPart in rxn.left_comp_ids:
    print(leftPart)
print("right")
for rightPart in rxn.right_comp_ids:
    print(rightPart)


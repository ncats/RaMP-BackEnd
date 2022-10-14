'''
Created on Jul 6, 2022

@author: braistedjc
'''
from rdflib import *
import rdflib
import sys, os
from os.path import exists
import pickle
import gzip
import shutil
from parse.MetabolomicsData import MetabolomicsData
from rampConfig.RampConfig import RampConfig

class ChebiOwlParser(MetabolomicsData):
    '''
    classdocs
    '''

    def __init__(self, resConfig):
        '''
        Constructor
        '''
        # used '../' for unit testing, else set to ""
        self.relDir = ""
        
        self.config = resConfig
        
        self.onto = None
        
        self.g = None
        
        self.relations = dict()
        
        self.chemEntityRelations = dict()
        
        self.localRelationsFile = ""
        
        self.localOntoFile = ""
        
        self.localOntoDir = ""
        
        self.outputDir = ""
        
        self.chebiStatus = dict()
        
        self.chebiOntoToLabel = dict()
        
        self.chebiOntoParentToChild = dict()
        
        self.chebiHumanIdSet = set()
        
        
    def createOwlGraph(self):
        print("creating owl graph")
            
    def getChebiFiles(self):
        print("fetching files")
        
        owlConfig = self.config.getConfig('chebi_ontology_owl')
        
        localDir = owlConfig.localDir
        
        if not exists(self.relDir + localDir):
            os.mkdir(self.relDir + localDir)

        #make an output dir for parsing
        self.outputDir = self.relDir + "../misc/output/chebi"
        if not exists(self.outputDir):
            os.mkdir(self.outputDir)
        

        # first RDF file
        owlLocalFile = owlConfig.extractFileName
                
        self.localOntoFile = self.relDir + localDir + owlLocalFile
        self.localOntoDir = self.relDir + localDir
                
        if not exists(self.relDir + localDir + owlLocalFile):
            # get the file
            owlUrl = owlConfig.sourceURL
            owlRemoteFile = owlConfig.sourceFileName
            
            self.download_files(owlUrl, self.relDir + localDir + owlRemoteFile)     

            # gunzip the file...
            with gzip.open(self.relDir + localDir + owlRemoteFile, 'rb') as f_in:
                with open(self.relDir + localDir + owlLocalFile, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)  
        
        
        # now get chebi relations file        
        relationsConfig = self.config.getConfig('chebi_to_chebi_relations')
            
        relLocalFile = relationsConfig.extractFileName
        
        self.localRelationsFile = self.relDir + localDir + relLocalFile
        
        print(relLocalFile)
        if not exists(self.relDir + localDir + relLocalFile):
            # get the file
            relUrl = relationsConfig.sourceURL
            relRemoteFile = relationsConfig.sourceFileName
            
            print(relUrl)
            print(relRemoteFile)
            print(localDir)
            self.download_files(relUrl, self.relDir + localDir + relRemoteFile) 
            
            
    def buildGraph(self):
        print("building graph")

        self.g = Graph()
        self.g.parse(self.localOntoFile, format="application/rdf+xml")
        
        print('finished parsing graph....')

        print('serializing graph to file')
        gFile = self.localOntoDir + "/chebi_graph.dat"
        with open(gFile, 'wb') as f:
            pickle.dump(self.g, f)
        
        print('finished graph serialized store')

    def deserializeGraph(self, graphFilePath):
        print("deserializing graph from file")   
        
        gFile = graphFilePath
        
        if gFile is None:
            gFile = self.localOntoDir + "/chebi_graph.dat"
        
        with open(gFile, 'rb') as f:
            self.g = pickle.load(f)
            
        print("finished graph deserialization")


    def buildRoleRelations(self):
        
        with open(self.localRelationsFile, 'r') as relFile:
            for line in relFile:
#            line = relFile.readline()
                #print(line)
                if line is not None and line != "":
                    terms = line.split("\t")
                    
                    if len(terms) == 5 and (terms[1] == 'is_a' or terms[1] == 'has_functional_parent'):
                        relSet = self.relations.get(terms[2], None)
                        if relSet is None:
                            relSet = [terms[3]]                        
                            self.relations[terms[2]] = relSet
                        else:
                            relSet.append(terms[3])
             
        # now lets build a dictionary starting with role 50906 on down.
        newRelations = dict()
        newRelations['50906'] = []
        self.getChildRoles('50906', self.relations, newRelations)
        
        termsSet = set()
        
        # now set the role relations as global relations object
        self.relations = newRelations
        parentTot = 0
        relTot = 0
        for r in self.relations:
            termsSet.add(r)
            parentTot = parentTot + 1
            for c in self.relations[r]:
                termsSet.add(c)
                relTot = relTot + 1
        
#         print(parentTot)
#         print(relTot)
        
#        roleKids = self.relations['75768']
#         for r in roleKids:
#             print(r)
            
        # export chebi role relations ontology terms
        with open(self.outputDir + '\chebi_role_and_function_ontology_relations.txt', 'w') as f:
            for parent in self.relations:
                children = self.relations[parent]
                for child in children:
                    f.write("chebi:"+parent+"\t"+"chebi:"+child+"\n")


    def buildMolecularEntityRelations(self):
        
        with open(self.localRelationsFile, 'r') as relFile:
            for line in relFile:
#            line = relFile.readline()
                #print(line)
                if line is not None and line != "":
                    terms = line.split("\t")
                    
                    if len(terms) == 5 and (terms[1] == 'is_a' or terms[1] == 'has_functional_parent'):
                        relSet = self.chemEntityRelations.get(terms[2], None)
                        if relSet is None:
                            relSet = [terms[3]]                        
                            self.chemEntityRelations[terms[2]] = relSet
                        else:
                            relSet.append(terms[3])
             
        # now lets build a dictionary starting with role 50906 on down.
        newRelations = dict()
        newRelations['24431'] = []
        self.getChildRoles('24431', self.chemEntityRelations, newRelations)
        
        termsSet = set()
        
        # now set the role relations as global relations object
        self.chemEntityRelations = newRelations
        parentTot = 0
        relTot = 0
        for r in self.chemEntityRelations:
            termsSet.add(r)
            parentTot = parentTot + 1
            for c in self.chemEntityRelations[r]:
                termsSet.add(c)
                relTot = relTot + 1
#             
        # export chebi role relations ontology terms
        with open(self.outputDir + '\chebi_chem_entity_ontology_relations.txt', 'w') as f:
            for parent in self.chemEntityRelations:
                children = self.chemEntityRelations[parent]
                for child in children:
                    f.write("chebi:"+parent+"\t"+"chebi:"+child+"\n")
    
    
    def extractCofactorStatus(self, chebiRoleId):
        cofactorData = self.collectMolecularEntitiesUnderRole(chebiRoleId, exportFileName = None, exportFile = False, relationTypes=['is_a', 'has_role'])
        return cofactorData[0]
    
    def collectMolecularEntitiesUnderRole(self, chebiRoleId, exportFileName, exportFile = False, relationTypes = ['is_a', 'has_role', 'has_functional_parent']):
        chemEntityRelations = dict()
        members = set()
        
        with open(self.localRelationsFile, 'r') as relFile:
            for line in relFile:
#            line = relFile.readline()
                #print(line)
                if line is not None and line != "":
                    terms = line.split("\t")
                    
#                    if len(terms) == 5 and (terms[1] == 'is_a' or terms[1] == 'has_functional_parent' or terms[1] == 'has_role'):
                    if len(terms) == 5 and (terms[1] in relationTypes):
                        
                        relSet = chemEntityRelations.get(terms[2], None)                                                        
                        if relSet is None:
                            relSet = [terms[3]]                        
                            chemEntityRelations[terms[2]] = relSet
                        else:
                            relSet.append(terms[3])
             
        # now lets build a dictionary starting with role 50906 on down.
        newRelations = dict()
        newRelations[chebiRoleId] = []
        self.getChildRoles(chebiRoleId, chemEntityRelations, newRelations)
        
        termsSet = set()
        
        chemEntityRelations = newRelations
        parentTot = 0
        relTot = 0
        for r in chemEntityRelations:
            termsSet.add(r)
            parentTot = parentTot + 1
            for c in chemEntityRelations[r]:
                termsSet.add(c)
                relTot = relTot + 1
#             
        # export chebi role relations ontology terms
        if exportFile:
            with open(self.outputDir + exportFileName, 'w') as f:
                for parent in chemEntityRelations:
                    children = chemEntityRelations[parent]
                    for child in children:
                        f.write("chebi:"+parent+"\t"+"chebi:"+child+"\n")
        
        else:
            for parent in chemEntityRelations:
                members.add("chebi:" + str(parent))
                # print("add member: "+str(parent))
                children = chemEntityRelations[parent]
                for child in children:
                    # print("add member: "+str(child))
                    members.add('chebi:' + str(child))
                        
            return [members, chemEntityRelations]            
                    
                                   
    def getChildRoles(self, parent, baseRelations, newRelations):
        children = baseRelations.get(parent, None)
        if children is not None:
            for child in children:
                newChildren = newRelations.get(parent, None)
                if newChildren is None:
                    newChildren = [child]
                    newRelations[parent] = newChildren
                else:
                    if child not in newChildren:
                        newChildren.append(child)    
             
                self.getChildRoles(child, baseRelations, newRelations)
     
    def extractRoleOntologyLabels(self):
        classPrefix = "http://purl.obolibrary.org/obo/"
        idSet = set()
        
        for parent in self.relations:
            idSet.add(parent)
            children = self.relations[parent]
            for child in children:
                idSet.add(child)
        
        query="""
        select $id $label 
        where {
            $id <http://www.w3.org/2000/01/rdf-schema#label> $label
        }"""    
        
        currQ = query
        
        print(currQ)
        
        idToRoleMapping = dict()
        
        res = self.g.query(currQ)
        for r in res:
            print(r[0])
            cid = self.uriToChebiNumber(r[0])
            print(cid)
            if cid is not None and cid in idSet:
                idToRoleMapping["chebi:"+cid] = r[1]

        with open(self.outputDir + '\chebi_role_id_to_labels.txt', 'w') as f:
            for cid in idToRoleMapping:                
                f.write(cid+"\t"+idToRoleMapping[cid]+"\n")
        
    
    def extractChemicalEntityOntologyLabels(self):
        classPrefix = "http://purl.obolibrary.org/obo/"
        idSet = set()
        
        for parent in self.chemEntityRelations:
            idSet.add(parent)
            children = self.chemEntityRelations[parent]
            for child in children:
                idSet.add(child)
        
        query="""
        select $id $label 
        where {
            $id <http://www.w3.org/2000/01/rdf-schema#label> $label
        }"""    
        
        currQ = query
        
        print(currQ)
        
        idToRoleMapping = dict()
        
        res = self.g.query(currQ)
        for r in res:
            print(r[0])
            cid = self.uriToChebiNumber(r[0])
            print(cid)
            if cid is not None and cid in idSet:
                idToRoleMapping["chebi:"+cid] = r[1]

        with open(self.outputDir + '\chebi_chem_entity_to_labels.txt', 'w') as f:
            for cid in idToRoleMapping:                
                f.write(cid+"\t"+idToRoleMapping[cid]+"\n")
                
                   
    def extractMetabolitesInSubclass(self, subclass):
        print("collecting metabolites in subclass")
        
    def extractHumanMetaboliteStatus(self):
        print("extracting human metabolites")
        
        class_pred = URIRef("http://www.w3.org/2002/07/owl#Class")
        subClass_predicate = URIRef('http://www.w3.org/2000/01/rdf-schema#subClassOf')        
        id_pred = URIRef("http://www.geneontology.org/formats/oboInOwl#id")
        restriction_pred = URIRef("http://www.w3.org/2002/07/owl#Restriction")
        
        someValues_pred = URIRef("http://www.w3.org/2002/07/owl#someValuesFrom")
        
        #someValues_pred = "owl:someValuesFrom"

        humanMetsOntol = "http://purl.obolibrary.org/obo/CHEBI_77746"
        humanMetsSerumOntol = "http://purl.obolibrary.org/obo/CHEBI_85234"
        humanMetsUrineOntol = "http://purl.obolibrary.org/obo/CHEBI_84087"
        humanMetsXenoOntol = "http://purl.obolibrary.org/obo/CHEBI_76967"
        
        humanMetsOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_77746")
        humanMetsSerumOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_85234")
        humanMetsUrineOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_84087")
        humanMetsXenoOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_76967")
 
        xenoOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_76206")
        drugMetOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_49103")
           
        humanMetsOntolList = [humanMetsOntol, humanMetsSerumOntol, humanMetsUrineOntol, humanMetsXenoOntol]
        
        humanMetsOntolPredList = [humanMetsOntol_pred, humanMetsSerumOntol_pred, humanMetsUrineOntol_pred, humanMetsXenoOntol_pred, drugMetOntol_pred]
        #xenoOntol_pred, drugMetOntol_pred]
        
        humanMetsOntoDict = dict()
        humanMetsOntoDict[humanMetsOntol_pred] = "human metabolite"
        humanMetsOntoDict[humanMetsSerumOntol_pred] = "human blood serum metabolite"
        humanMetsOntoDict[humanMetsUrineOntol_pred] = "human urinary metabolite"
        humanMetsOntoDict[humanMetsXenoOntol_pred] = "human xenobiotic metabolite"
        humanMetsOntoDict[drugMetOntol_pred] = "drug metabolite"        
        
        t2 = 0
        xlist = []
        xdict = dict()
        for s, p, o in self.g:
            if o in humanMetsOntolPredList:
                par = self.g.subjects(object=s)
                for x in par:
#                     print("have parent...")
#                     print(x)
                    if not isinstance(x, rdflib.BNode):
                        if x not in xlist:
                            xlist.append(x)
                            t2 = t2 + 1
                            xdict[x] = o
        
        
                    
        # get child nodes
        for x in xlist:
            nodes = self.getChildNodes(x, xlist, xdict, subClass_predicate, xdict[x])
        
        xlist = nodes
        
        print("finished extracting human chebi metabolites")
        print(str(len(set(xlist))))
        print(self.localOntoDir)
                
        with open(self.outputDir + '\human_chebi_ids.txt', 'w') as f:
            for line in xlist:
                f.write(self.uriToChebiId(line)+"\n")
                
                # store a list of human chebi's
                self.chebiHumanIdSet.add(self.uriToChebiId(line))

#         with open(self.outputDir + '\human_chebi_ids_with_class.txt', 'w') as f:
#             for chebiId in xdict:
#                 cclass = xdict[chebiId]
#                 cid = "chebi:" + chebiId.split("_")[1].strip()
#                 clabel = humanMetsOntoDict[cclass]
#                 f.write(cid+"\thuman metabolite\n")

        
    def getChildNodes(self, parentNode, xlist, xdict, subClass_predicate, parentClass):    
        
        z = self.g.subjects(predicate=subClass_predicate, object=parentNode) 
                
        for u in z:
            if not isinstance(u, rdflib.BNode) and u not in xlist:
                xlist.append(u)
                xdict[u] = parentClass
                self.getChildNodes(u, xlist, xdict, subClass_predicate, parentClass)
        
        return xlist

    def extractChebiOntologyStructure(self):
        print("extracting chebi ontology mappings")
        
        subset_pred = URIRef("http://www.geneontology.org/formats/oboInOwl#inSubset")
        someValues_pred = URIRef("http://www.w3.org/2002/07/owl#someValuesFrom")
        subClass_predicate = URIRef('http://www.w3.org/2000/01/rdf-schema#subClassOf')
        label_pred = URIRef("http://www.w3.org/2000/01/rdf-schema#label")
        threeStars_obj = URIRef("http://purl.obolibrary.org/obo/chebi#3_STAR")
        twoStars_obj = URIRef("http://purl.obolibrary.org/obo/chebi#2_STAR")
        
        threeStarComps = self.g.subjects(predicate=subset_pred, object=threeStars_obj)
        twoStarComps = self.g.subjects(predicate=subset_pred, object=twoStars_obj)

        onto_pred = URIRef("http://www.geneontology.org/formats/oboInOwl#hasOBONamespace")

        ontologyPropertyObj = URIRef('http://purl.obolibrary.org/obo/RO_0000087')
        


        ########
        #
        # Note: This section gets all ontology terms that are actually linked to compound records
        # This is useful for associations, but will miss high level terms that are not attached to compound
        # Following from child to parent terms will extend this in sections below, but it has to be done 
        # more than one level deep.
        #
        
        q = """
        SELECT DISTINCT ?onto ?label
        WHERE {
            ?onto <http://www.w3.org/2000/01/rdf-schema#label> $label .
            ?x <http://www.w3.org/2002/07/owl#someValuesFrom> $onto .
            ?x <http://www.w3.org/2002/07/owl#onProperty> <http://purl.obolibrary.org/obo/RO_0000087>
        }"""
                
        res = self.g.query(q)
        rc = 0
        idToName = dict()
        ontoIds = []
        parentToChild = dict()
        for r in res:
            rc = rc + 1
            ontoIds.append(r[0])
            idToName[r[0]] =r[1]
            
            # get the subjects (class entry) for each term node
            r2 = self.g.subjects(object=r[0])
            
            children = parentToChild.get(r[0], None)
            if children is None:
                parentToChild[r[0]] = []
                
            # go over all term nodes
            for r3 in r2:
                
                # collect parent nodes
                r4 = self.g.objects(subject=r3, predicate=subClass_predicate)
                
                # for each parent node, add the child node
                #
                # HERE we add the next generation parents... but hey... we are not recursively going up the chain to get
                # the reset of the hierarchy, beyond immediate parents.
                # 
                #
                                
                for r5 in r4:
                    if not isinstance(r5, rdflib.BNode):
                        #print(str(r3) + " " + str(r5))
                        currentChildren = parentToChild.get(r5, None)
                        if currentChildren is None:
                            currentChildren = []
                            parentToChild[r5] = currentChildren
                        if r3 not in currentChildren:
                            currentChildren.append(r3)
                            
                        # need to work on parent r5
                        # this will break into a recursive call to...
                        self.getParentOntologyTerms(r5, parentToChild, subClass_predicate)
                    
        #print("class type = " + str(type(res)))

        termCount = 0
        vals = []
        for val in parentToChild:
            termCount = termCount + 1
            #print(val)
            vals.append(val)
            
        print(termCount)
        print(rc)
        extra = set(vals).difference(set(ontoIds))
        print(str(len(extra)))
        
        # if we have references to extra ontology terms that are not associated to compounds
        # traverse and add to id list and get the label
        for e in extra:
            # print(e)
            ontoIds.append(e)
            labels = self.g.objects(subject=e, predicate=label_pred)
            for label in labels:
                idToName[e] = label

            # hey lets check that we have parent terms for these guys...
            self.getParentOntologyTerms(e, parentToChild, subClass_predicate)

            
        rc = 0
        for t in idToName:
            rc = rc + 1
        
        pcc = 0
        for p in parentToChild:
            pcc = pcc + 1
        
        print(termCount)
        print(rc)
        print(pcc)
        
        # export id to name mapping
        with open(self.outputDir + '\chebi_ontology_terms.txt', 'w') as f:
            for chebiId in idToName:
                term = idToName[chebiId]                                                 
                f.write(self.uriToChebiId(chebiId)+"\t"+self.uriToValue(term)+"\n") 

        # export id to id term to terms
        with open(self.outputDir + '\chebi_ontology_term2terms_parent_child_mapping.txt', 'w') as f:
            for parentId in parentToChild:
                childIds = parentToChild[parentId]
                for child in childIds:                                           
                    f.write(self.uriToChebiId(parentId)+"\t"+self.uriToChebiId(child)+"\n")
                    
        # store to class fields
        self.chebiOntoToLabel = idToName
        
        self.chebiOntoParentToChild = parentToChild        
    
    def uriToChebiId(self, uri):
        return ('chebi:' + uri.split("_")[1]).strip()
    
    def uriToValue(self, uri):
        vals = uri.split("/")
        return vals[len(vals)-1]
    
    def uriToChebiNumber(self, uri):
        val = None
        if '_' in uri: 
            val = (uri.split("_")[1]).strip()    
        return val
    
    def getParentOntologyTerms(self, childNode, parentToChild, subClass_predicate):
        parents = self.g.objects(subject=childNode, predicate=subClass_predicate)
        for parent in parents:
            if not isinstance(parent, rdflib.BNode):
                children = parentToChild.get(parent, None)
                if children is None:
                    children = []
                    parentToChild[parent] = children
                if childNode not in children: 
                    children.append(childNode)
                    
                # now go get the parent's parents
                self.getParentOntologyTerms(parent, parentToChild, subClass_predicate)    
    
    
    def extractChebiOntologyAssociations(self):
        print("extracting chebi ontology mappings")
        
        subset_pred = URIRef("http://www.geneontology.org/formats/oboInOwl#inSubset")
        someValues_pred = URIRef("http://www.w3.org/2002/07/owl#someValuesFrom")
        subClass_predicate = URIRef('http://www.w3.org/2000/01/rdf-schema#subClassOf')
        label_pred = URIRef("http://www.w3.org/2000/01/rdf-schema#label")
        threeStars_obj = URIRef("http://purl.obolibrary.org/obo/chebi#3_STAR")
        twoStars_obj = URIRef("http://purl.obolibrary.org/obo/chebi#2_STAR")
        
        threeStarComps = self.g.subjects(predicate=subset_pred, object=threeStars_obj)
        twoStarComps = self.g.subjects(predicate=subset_pred, object=twoStars_obj)

        ontologyPropertyObj = URIRef('http://purl.obolibrary.org/obo/RO_0000087')
        
        q = """
        SELECT DISTINCT ?onto ?label
        WHERE {
            ?onto <http://www.w3.org/2000/01/rdf-schema#label> $label .
            ?x <http://www.w3.org/2002/07/owl#someValuesFrom> $onto .
            ?x <http://www.w3.org/2002/07/owl#onProperty> <http://purl.obolibrary.org/obo/RO_0000087>
        }"""
        
        res = self.g.query(q)
        rc = 0
        for r in res:
            rc = rc + 1
            print(r)
         
        print("query res size:"+str(rc))   
        
        # allTopStarComps = threeStarComps + twoStarComps

        ontoDict = dict()
        
        # these get primary mappings, secondary mappings are inherited by descendants 
        for c in threeStarComps:
            ontoVals = self.g.objects(subject=c, predicate=subClass_predicate)
            
            onto = ontoDict.get(c, None)
            ontLen = 0
            if onto is None:
                onto = []
                ontoDict[c] = onto
                
            for node in ontoVals:
                if isinstance(node, rdflib.BNode):
                    node
                    ontoIds = self.g.objects(subject=node, predicate=someValues_pred)
                    for ontoVal in ontoIds:
                        onto.append(ontoVal)
                        ontLen = ontLen + 1
        
#             if ontLen > 0:
#                 if c == URIRef('http://purl.obolibrary.org/obo/CHEBI_33448'):
#                     print(c)
#                     print(ontLen)
#                     ontos = ontoDict[c]
#                     for ont in ontos:
#                         print(ont)
        print(str(sum(1 for _ in ontoDict)))
        
        for c in twoStarComps:
            ontoVals = self.g.objects(subject=c, predicate=subClass_predicate)
            ontLen = 0
            onto = ontoDict.get(c, None)
            if onto is None:
                onto = []
                ontoDict[c] = onto
                        
            for node in ontoVals:
                if isinstance(node, rdflib.BNode):
                    ontoIds = self.g.objects(subject=node, predicate=someValues_pred)
                    for ontoVal in ontoIds:
                        onto.append(ontoVal)
                        ontLen = ontLen + 1
        
            #print(ontLen)
        
        print(str(sum(1 for _ in ontoDict)))
        
        # now traverse for descendants of these primary nodes
        # get a copy of the keys since ontoDict will be updated dynamically
        ontoDictCopy = ontoDict.copy()
        keysCopy = ontoDictCopy.keys() 
        for parent in ontoDictCopy:
            self.getChildOntologies(parent, ontoDict, subClass_predicate, keysCopy)
        
                
#         res = ontoDict.get(URIRef("http://purl.obolibrary.org/obo/CHEBI_84095"), None)
#         if res is not None:
#             for r in res:
#                 print(r)
            
        print(str(sum(1 for _ in ontoDict)))
        
        numWithOntology = 0
        ontologySet = set()
        for c in ontoDict:
            if len(ontoDict[c]) > 0:
                numWithOntology = numWithOntology + 1;
                for ontology in ontoDict[c]:
                    ontologySet.add(ontology)
    
        print(str(numWithOntology))
        print(str(len(ontologySet)))
        
        with open(self.outputDir + '\chebi_ontologies.txt', 'w') as f:
            for chebiId in ontoDict:
                onto = ontoDict[chebiId]
                cid = "chebi:" + chebiId.split("_")[1].strip() 
                for term in onto:                 
                    f.write(cid+"\t"+term+"\n")
        
    def getChildOntologies(self, parentNode, ontoDict, subClass_predicate, compInclusionList):    
        
        z = self.g.subjects(predicate=subClass_predicate, object=parentNode)
        
        parentOntol = ontoDict[parentNode]       
        #print("hey we have parent ontolgoy of length...:" + str(len(parentOntol)))
        for u in z:
            if not isinstance(u, rdflib.BNode) and u in compInclusionList:
                currList = ontoDict.get(u, None)
                if currList is None:
                    currList = parentOntol
                    ontoDict[u] = currList
                    # recursion if we haven't visited here before...
                    self.getChildOntologies(u, ontoDict, subClass_predicate)
                    print("Hey we're recursing to get child ontologies? " + u)                
                else:
                    for ontol in parentOntol:
                        if ontol not in currList:
                            currList.append(ontol)
            
#            if u in allComps:               
                

        return


    def extractHumanMetaboliteStatusTest(self):
        
        print("extracting human metabolites")
        
        class_pred = URIRef("http://www.w3.org/2002/07/owl#Class")
        subClass_predicate = URIRef('http://www.w3.org/2000/01/rdf-schema#subClassOf')        
        id_pred = URIRef("http://www.geneontology.org/formats/oboInOwl#id")
        restriction_pred = URIRef("http://www.w3.org/2002/07/owl#Restriction")
        
        someValues_pred = URIRef("http://www.w3.org/2002/07/owl#someValuesFrom")
        
        #someValues_pred = "owl:someValuesFrom"

        humanMetsOntol = "http://purl.obolibrary.org/obo/CHEBI_77746"
        humanMetsSerumOntol = "http://purl.obolibrary.org/obo/CHEBI_85234"
        humanMetsUrineOntol = "http://purl.obolibrary.org/obo/CHEBI_84087"
        humanMetsXenoOntol = "http://purl.obolibrary.org/obo/CHEBI_76967"
        
        humanMetsOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_77746")
        humanMetsSerumOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_85234")
        humanMetsUrineOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_84087")
        humanMetsXenoOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_76967")
 
        xenoOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_76206")
        drugMetOntol_pred = URIRef("http://purl.obolibrary.org/obo/CHEBI_49103")
           
        humanMetsOntolList = [humanMetsOntol, humanMetsSerumOntol, humanMetsUrineOntol, humanMetsXenoOntol]
        
        humanMetsOntolPredList = [humanMetsOntol_pred, humanMetsSerumOntol_pred, humanMetsUrineOntol_pred, humanMetsXenoOntol_pred]
        #, xenoOntol_pred, drugMetOntol_pred]

        
        humanMetsOntoDict = dict()
        humanMetsOntoDict[humanMetsOntol] = "human metabolite"
        humanMetsOntoDict[humanMetsSerumOntol] = "human blood serum metabolite"
        humanMetsOntoDict[humanMetsUrineOntol] = "human urinary metabolite"
        humanMetsOntoDict[humanMetsXenoOntol] = "human xenobiotic metabolite"
        
        #res = self.g.subjects(predicate=type_predicate, object=reaction_object)
                
        s = self.g.objects(subject=URIRef("http://purl.obolibrary.org/obo/CHEBI_100"), predicate=URIRef("http://purl.obolibrary.org/obo/chebi/inchikey"))
        
        subset_pred = URIRef("http://www.geneontology.org/formats/oboInOwl#inSubset")
        threeStars_obj = URIRef("http://purl.obolibrary.org/obo/chebi#3_STAR")
        
        allComps = self.g.subjects(predicate=subset_pred, object=threeStars_obj)
       
       
        # try a sparql query
        q = """
        SELECT DISTINCT ?Class 
        WHERE {
            ?Class <http://www.w3.org/2002/07/owl#someValuesFrom> <http://purl.obolibrary.org/obo/CHEBI_77746>
        }"""
       
        rc = 0
        res = self.g.query(q)
        for row in res:
            # print(row)
            rc = rc + 1
       
        # print(str(rc))
        
        humanMetsOntol_obj = URIRef("http://purl.obolibrary.org/obo/CHEBI_77746")
        humanMetsSerumOntol_obj = URIRef("http://purl.obolibrary.org/obo/CHEBI_85234")
        humanMetsUrineOntol_obj = URIRef("http://purl.obolibrary.org/obo/CHEBI_84087")
        humanMetsXenoOntol_obj = URIRef("http://purl.obolibrary.org/obo/CHEBI_76967")

        humanMets = self.g.subjects(predicate = someValues_pred, object = humanMetsOntol_obj)
        # print("HMs")
        #hmc = sum(1 for _ in humanMets)
        #print(hmc)
        #print(type(humanMets))
        hmc = 0
        for met in humanMets:            
            #print(met)
            #print(type(met))
            hmc = hmc + 1
            
            secondMets = self.g.subjects(object=met, predicate=class_pred)
            
            for m in secondMets:
                #print(m)
                hmc = hmc + 1
#             sc = self.g.subjects(object=met, predicate=class_pred)
#             for s in sc:
#                 hmc = hmc + 1

        print(hmc)


        humanMets = self.g.subjects(predicate = someValues_pred, object =humanMetsSerumOntol_obj)
        print("HMs Serum")
        print(sum(1 for _ in humanMets))
        
        humanMets = self.g.subjects(predicate = someValues_pred, object =humanMetsUrineOntol_obj)
        print("HMs Urine")
        print(sum(1 for _ in humanMets))
        
        humanMets = self.g.subjects(predicate = someValues_pred, object = humanMetsXenoOntol_obj)
        print("HMs Xeno")
        print(sum(1 for _ in humanMets))

        
        i = 0
        numHuman = 0
        numNonHuman = 0
        hc = 0
        restrictionCount = 0
        
        for ssub in allComps:
            i = i + 1
            #print(ssub)
            
            subClasses = self.g.objects(subject=ssub, predicate=subClass_predicate)
            chebiIds = self.g.objects(subject=ssub, predicate=id_pred)
            restrictions = self.g.objects(subject=ssub, predicate=restriction_pred)
            someValuesFrom = self.g.objects(subject=ssub, predicate = someValues_pred)

            #for sv in someValuesFrom:
            #    print(sv)
            
#             if restrictions is not None:
#                 print("restriction")
#                 print(type(restrictions))
#                 print(str(sum(1 for _ in restrictions)))
#                 print(restrictions)
                # try to use the restriction as a subject.... :)
#                 
#                 someValuesFrom = self.g.objects(subject=restrictions, predicate = someValues_pred)
#                 
#                 for r in someValuesFrom:
#                     print("values????")
#                     print(type(r))
#                     print(r)
#                     restrictionCount = restrictionCount + 1
            

            chebiId = ""
            for cid in chebiIds:
                chebiId = cid
                #print(chebiId)
            
            for subc in subClasses:
                #print(subc)
                #print(type(subc))                 
#                 if isinstance(subc, rdflib.BNode):
#                     for node in subc:
#                         print(type(node))
#                         print(node)

                
#                 print(ssub)
#                 print(subc)
#                 print("svf size")
#                 print(someValuesFrom)
#                 print(str(sum(1 for _ in someValuesFrom)))
                
#                 for s in someValuesFrom:
#                     # print("hey some values..")
#                     # print(s)                    
#                     # c = humanMetsOntoDict.get(s, None)
#                     #if c is not None:
#                     if s in humanMetsOntolList:
#                         print("hey human met???????")
#                         print(humanMetsOntoDict[s])
#                         hc = hc + 1
                
                metType = humanMetsOntoDict.get(subc, None)
                if chebiId != None:
                    if metType == None:
                        self.chebiStatus[chebiId] = "None"
                        numNonHuman = numNonHuman + 1                    
                    else:
                        self.chebiStatus[chebiId] = metType
                        #print(chebiId+" "+metType)
                        numHuman = numHuman + 1
                        break
        

                
# rConf = RampConfig()
# rConf.loadConfig("../../config/external_resource_config.txt")
#                         
# cop = ChebiOwlParser(rConf)   
# cop.getChebiFiles()
# #cop.buildGraph()
# cop.deserializeGraph(None)
# #cop.extractChebiOntologyAssociations()
# #cop.extractChebiOntologyStructure()
# cop.extractHumanMetaboliteStatus()

#cop.buildRoleRelations()
#cop.buildMolecularEntityRelations()
#cop.extractRoleOntologyLabels()
#cop.extractChemicalEntityOntologyLabels()



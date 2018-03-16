import urllib.request as RE
import time
import os
import sys
from MetabolomicsData import MetabolomicsData
import zipfile
import path
from rdflib import URIRef,Graph
import rdflib.namespace
from rdflib.namespace import RDF, FOAF,RDFS,DC,DCTERMS

class WikipathwaysRDF(MetabolomicsData):
    def __init__(self):
        super().__init__()
        # key: ID for metabolites Value: Common Name (the only name in this database)
        self.metaboliteCommonName = dict()
        #key: ID for pathway, Value: pathway name
        self.pathwayDictionary = dict()
        
        #key: ID for pathway, category: various categories such as cellular process, metabolic process 
        self.pathwayCategory = dict()
        
        #key: gene, value: gene mapping
        self.geneInfoDictionary = dict()
        
        #key: metabolite, value: metabolite mapping
        self.metaboliteIDDictionary = dict()
        
        #key: pathwayID, value: list of genes
        self.pathwaysWithGenesDictionary = dict()
        
        #key: pathwayId, value: list of metabolites
        self.pathwayWithMetabolitesDictionary = dict()
        
        #empty for reactome
        self.metabolitesWithSynonymsDictionary = dict()
        
        #only not empty when a catalyzed class exists 
        #empty
        self.metabolitesLinkedToGenes = dict()
        
        #key: metabolite, value: list of pathways
        self.metabolitesWithPathwaysDictionary = dict()
        
        #key: current pathway, value: list of pathways linked to this pathway
        self.pathwayOntology = dict()
        
        #stays empty for this class
        self.biofluidLocation = dict()
        
        #stays empty for this class
        self.biofluid = dict()
        
        #stays empty for this class
        self.cellularLocation = dict()
        
        #stays empty for this class
        self.cellular = dict()
        
        #stays empty
        self.exoEndoDictionary = dict()
        self.exoEndo =dict()
        #tissue location stay empty
        self.tissue = dict()
        self.tissueLocation = dict()
    def getDatabaseFile(self):
        '''
        Downloaded wikipathway file from the given url
        The url is from current(latest) version of wikipathways
        if the file name is wrong, go the the url to check if the file is updated
        '''
        url = 'http://data.wikipathways.org/current/rdf/'
        filename = 'wikipathways-20180310-rdf-wp.zip'
        path = '../misc/data/wikipathwaysRDF/'
        self.check_path(path)
        existed = os.listdir(path)
        if filename not in existed:
            print(path + filename)
            self.download_files(url+filename, path + filename)
            with zipfile.ZipFile(path+filename,'r') as zip_ref:
                zip_ref.extractall(path)
        else:
            print('Already download ...')
    
    def _getAllRDFTypes(self):
        '''
        Find all possible types form rdf files
        All namespaces for the rdf files:
        ('rdfs', rdflib.term.URIRef('http://www.w3.org/2000/01/rdf-schema#'))
        ('dc', rdflib.term.URIRef('http://purl.org/dc/elements/1.1/'))
        ('biopax', rdflib.term.URIRef('http://www.biopax.org/release/biopax-level3.owl#'))
        ('foaf', rdflib.term.URIRef('http://xmlns.com/foaf/0.1/'))
        ('dcterms', rdflib.term.URIRef('http://purl.org/dc/terms/'))
        ('wp', rdflib.term.URIRef('http://vocabularies.wikipathways.org/wp#'))
        ('void', rdflib.term.URIRef('http://rdfs.org/ns/void#'))
        ('xml', rdflib.term.URIRef('http://www.w3.org/XML/1998/namespace'))
        ('gpml', rdflib.term.URIRef('http://vocabularies.wikipathways.org/gpml#'))
        ('prov', rdflib.term.URIRef('http://www.w3.org/ns/prov#'))
        ('wprdf', rdflib.term.URIRef('http://rdf.wikipathways.org/'))
        ('rdf', rdflib.term.URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#'))
        ('skos', rdflib.term.URIRef('http://www.w3.org/2004/02/skos/core#'))
        ('freq', rdflib.term.URIRef('http://purl.org/cld/freq/'))
        ('xsd', rdflib.term.URIRef('http://www.w3.org/2001/XMLSchema#'))
        ('pav', rdflib.term.URIRef('http://purl.org/pav/'))
        '''
        rdftype = set()
        path = '../misc/data/wikipathwaysRDF/wp/Human/'
        listoffiles = os.listdir(path)
        print(listoffiles)
        listoffiles.sort()
        print('Total {} pathways in Human'.format(len(listoffiles)))
        for each in listoffiles:
            print(each)
            g = Graph()
            g.parse(path+each,format = 'n3')
            print('length of graph: {}'.format(len(g)))
            for s,p,o in g.triples((None, RDF.type,None)):
                #print('---------------------------')
                #print('Subject: {} \nPredicate: {}\nObject: {}'.format(s,p,o))
                #print('---------------------------')
                obj = o[o.find('wp#') + 3:]
                if obj == -1:
                    raise ValueError('Wrong object value from rdf')
                if obj not in rdftype:
                    rdftype.add(obj)
                #time.sleep(1)
        
        print(rdftype)
        return(rdftype)
    '''
    All possible types 
    {'Metabolite', 'TranscriptionTranslation', 'DataNode', 'Stimulation', 
    'Catalysis', 'DirectedInteraction', 'Pathway', 'Conversion', 'ComplexBinding'
    , 'Inhibition', 'Protein', 'tp://www.w3.org/2004/02/skos/core#Collection', 
    'Rna', 'Interaction', 'Binding', 'PublicationReference', 'GeneProduct', 'Complex'}        
    '''
    def getAllPathways(self):
        path = '../misc/data/wikipathwaysRDF/wp/Human/'
        listoffiles = os.listdir(path)
        print('Total {} pathways in Human'.format(len(listoffiles)))
        
        for each in listoffiles:
            this_pathway = each.replace('.ttl','')
            print('pathway id is {}'.format(this_pathway))
            g = Graph()
            g.parse(path + each,format = 'n3')
            for s,p,o in g.triples((None,DC.title,None)):
                #print('---------------------------')
                #print('Subject: {} \nPredicate: {}\nObject: {}'.format(s,p,o))
                #print('---------------------------')
                self.pathwayDictionary[this_pathway] = o
                self.pathwayCategory[this_pathway] = 'NA'
        
        print('Total pathways are {}'.format(len(self.pathwayDictionary)))
    def getIDmapping(self): 
        path = '../misc/data/wikipathwaysRDF/wp/Human/'
        listoffiles = os.listdir(path)
        print('Total {} pathways in Human'.format(len(listoffiles)))
        
        for each in listoffiles:
            this_pathway = each.replace('.ttl','')
            print('pathway id is {}'.format(this_pathway))
            g = Graph()
            g.parse(path + each,format = 'n3')
            
            analytes = []
            for s,p,o in g.triples((None,DCTERMS.identifier,None)):
                #print('---------------------------')
                #print('Subject: {} \nPredicate: {}\nObject: {}'.format(s,p,o))
                #print('---------------------------')
                analytes.append(s)
                
            print('Total analytes in this pathway: {}'\
                  .format(len(analytes)))
            time.sleep(0.5)
            pred = URIRef('http://vocabularies.wikipathways.org/wp#')
            while len(analytes)>0:
                id = analytes.pop(0)
                for s,p,o in g.triples((id,None,None)):
                    print('---------------------------')
                    print('Subject: {} \nPredicate: {}\nObject: {}'.format(s,p,o))
                    print('---------------------------')
                    time.sleep(1)
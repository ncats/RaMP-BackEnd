import sys
sys.path.append('/Users/pati13/Downloads/')
import urllib.request as RE
import time
import os
from MetabolomicsData import MetabolomicsData
import zipfile
from rdflib import URIRef,Graph
import rdflib.namespace
from rdflib.namespace import RDF, FOAF,RDFS,DC,DCTERMS
from builtins import str



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
        self.pathwaysWithMetabolitesDictionary = dict()
        
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
    def getEverything(self,writeToFile = False):
        '''
        This function pack all functions in this class together, running this function will parse all data 
        we would like to have.
        - param bool writeToFile if true, write all dictionary to misc/output/wikipathwayRDF/
        '''
        self.getDatabaseFile()
        self.getIDMapingWithPathways()
        if writeToFile:
            self.write_myself_files('wikipathwayRDF')
    def getDatabaseFile(self):
        '''
        Downloaded wikipathway file from the given url
        The url is from current(latest) version of wikipathways
        if the file name is wrong, go the the url to check if the file is updated
        '''
        url = 'http://data.wikipathways.org/current/rdf/'
        filename = 'wikipathways-20190210-rdf-wp.zip'
        path = '../misc/data/wikipathwaysRDF/'
        self.check_path(path)
        existed = os.listdir(path)
        if filename not in existed:
            #print(path + filename)
            self.download_files(url+filename, path + filename)
            with zipfile.ZipFile(path+filename,'r') as zip_ref:
                zip_ref.extractall(path)
        else:
            print('Already downloaded Wiki ...')
    
    def _getAllRDFTypes(self):
        '''
        Find all possible types form rdf files
        return:
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
        #print(listoffiles)
        listoffiles.sort()
        #print('Total {} pathways in Human'.format(len(listoffiles)))
        for each in listoffiles:
            #print(each)
            g = Graph()
            g.parse(path+each,format = 'n3')
            #print('length of graph: {}'.format(len(g)))
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
        
        #print(rdftype)
        return(rdftype)
    '''
    All possible types 
    {'Metabolite', 'TranscriptionTranslation', 'DataNode', 'Stimulation', 
    'Catalysis', 'DirectedInteraction', 'Pathway', 'Conversion', 'ComplexBinding'
    , 'Inhibition', 'Protein', 'tp://www.w3.org/2004/02/skos/core#Collection', 
    'Rna', 'Interaction', 'Binding', 'PublicationReference', 'GeneProduct', 'Complex'}        
    '''
    def getPathwayInfoFromGraph(self,g,this_pathway):
        for s,p,o in g.triples((None,DC.title,None)):
            '''
                print('---------------------------')
                print('Subject: {} \nPredicate: {}\nObject: {}'.format(s,p,o))
                print('---------------------------')
            '''
            self.pathwayDictionary[this_pathway] = o
            self.pathwayCategory[this_pathway] = 'NA'
            
    def displayRDFfile(self,second = 3):
        path = '../misc/data/wikipathwaysRDF/wp/Human/'
        listoffiles = os.listdir(path)
        #print('Total {} pathways in Human'.format(len(listoffiles)))
        
        for each in listoffiles:
            this_pathway = each.replace('.ttl','')
            #print('pathway id is {}'.format(this_pathway))
            g = Graph()
            g.parse(path + each,format = 'n3')
            for s,p,o in g.triples((None,None,None)):
                
                print('---------------------------')
                print('Subject: {} \nPredicate: {}\nObject: {}'.format(s,p,o))
                print('---------------------------')
                time.sleep(second)
        
    def getIDMapingWithPathways(self): 
        '''
        This function parse all RDF files in Human pathways from source files. Then call functions
    
        1) self.getPathwayInfoFromGraph(g, this_pathway)
        2) self.getMetabolitesIDFromGraph(g,this_pathway)
        3) self.getGenesIDFromGraph(g, this_pathway)
        4) self.getCatalyzation(g, this_pathway) // Not implemented 
        
        '''
        path = '../misc/data/wikipathwaysRDF/wp/Human/'
        self.check_path(path)
        listoffiles = os.listdir(path)
        #print('Total {} pathways in Human'.format(len(listoffiles)))
        
        
        total_files = len(listoffiles)
        i = 0
        for each in listoffiles:
            i = i + 1
            
            # get pathways id
            this_pathway = each.replace('.ttl','')
            #print('{}/{} ID:{}'.format(i,total_files,this_pathway))
            g = Graph()
            g.parse(path + each,format = 'n3')
            
            # get pathway information at first
            self.getPathwayInfoFromGraph(g, this_pathway)
            # get metabolites information at second
            self.getMetabolitesIDFromGraph(g,this_pathway)
            self.getGenesIDFromGraph(g, this_pathway)
            #self.getCatalyzation(g, this_pathway)
    def getGenesIDFromGraph(self,g,this_pathway):
        geneProduct = URIRef('http://vocabularies.wikipathways.org/wp#GeneProduct')
        type_predicate = URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
        proteins = URIRef('http://vocabularies.wikipathways.org/wp#Protein')
        possible_source = {
            'ncbigene':'Entrez',
            'ensembl':'Ensembl',
            'uniprot':'Uniprot',
            'ec-code':'Enzyme Nomenclature',
            'wikidata':'WikiData',
            'ncbiprotein':'NCBI-ProteinID'
            }
        # These source are not retrieved at this moment
        not_retrieved = ['wikipedia.en','mirbase','hgnc.symbol','ena.embl','mirbase.mature','kegg.genes','go',
                         'interpro','refseq','pfam','ecogene','chembl.compound']
        genelist = set()
        for gene in g.subjects(type_predicate,geneProduct):
            bdbLinks = {
                'Entrez':'http://vocabularies.wikipathways.org/wp#bdbEntrezGene',
                'UniProt':'http://vocabularies.wikipathways.org/wp#bdbUniprot',
                'Ensembl':'http://vocabularies.wikipathways.org/wp#bdbEnsembl',
                'WikiData':'http://vocabularies.wikipathways.org/wp#bdbWikiData'
            }
            geneMapping = {"kegg": "NA",
                         "common_name": "NA",
                         "Ensembl": "NA", 
                         "HGNC": "NA", 
                         "HPRD": "NA", 
                         "NCBI-GeneID": "NA", 
                         "NCBI-ProteinID": "NA", 
                         "OMIM": "NA", 
                         "UniProt": "NA", 
                         "Vega": "NA", 
                         "miRBase": "NA", 
                         "HMDB_protein_accession": "NA",
                         "Entrez" : "NA",
                         "Enzyme Nomenclature": "NA",
                         'WikiData':'NA'}
            genesource = gene.split('/')[-2]
            if genesource not in possible_source:
                continue
            geneid = gene.split('/')[-1]
            geneid = self.prependID(genesource, geneid)
            if genesource not in not_retrieved:
                geneMapping[possible_source[genesource]] = [geneid]
                genelist.add(geneid)
                for key,value in bdbLinks.items():
                    for links in g.objects(gene,URIRef(value)):
                        link_id = links.split('/')[-1]
                        #print("link id:", link_id)
                        #if link_id is "TP53":
                            #print("*****************************found TP53 wiki *****")
                        link_id = self.prependID(key, link_id)
                        #print("link id:", link_id)
                        # sometimes URIREF type object accidently appears in link_id var, so avoid it by checking
                        # data type here
                        genelist.add(link_id)
                        if geneMapping[key] == 'NA' and type(link_id) is str:
                            geneMapping[key] = [link_id]
                        elif type(geneMapping[key]) is list and type(link_id) is str:
                            if link_id not in geneMapping[key]:
                                geneMapping[key].append(link_id)
                #print('{}\n{}'.format(geneid,geneMapping))
                self.geneInfoDictionary[geneid] = geneMapping             
                
        for protein in g.subjects(type_predicate,proteins):
            bdbLinks = {
                'Entrez':'http://vocabularies.wikipathways.org/wp#bdbEntrezGene',
                'UniProt':'http://vocabularies.wikipathways.org/wp#bdbUniprot',
                'Ensembl':'http://vocabularies.wikipathways.org/wp#bdbEnsembl',
                'WikiData':'http://vocabularies.wikipathways.org/wp#bdbWikidata'
            }
            geneMapping = {"kegg": "NA",
                         "common_name": "NA",
                         "Ensembl": "NA", 
                         "HGNC": "NA", 
                         "HPRD": "NA", 
                         "NCBI-GeneID": "NA", 
                         "NCBI-ProteinID": "NA", 
                         "OMIM": "NA", 
                         "UniProt": "NA", 
                         "Vega": "NA", 
                         "miRBase": "NA", 
                         "HMDB_protein_accession": "NA",
                         "Entrez" : "NA",
                         "Enzyme Nomenclature": "NA",
                         'WikiData':'NA'}
            genesource = protein.split('/')[-2]
            geneid = protein.split('/')[-1]
            geneid = self.prependID(genesource, geneid)
            
            if genesource not in not_retrieved:
                geneMapping[possible_source[genesource]] = [geneid]
                genelist.add(geneid)
                for key,value in bdbLinks.items():
                    for links in g.objects(protein,URIRef(value)):
                        link_id = links.split('/')[-1]
                        link_id = self.prependID(key, link_id)
                        # sometimes URIREF type object accidently appears in link_id var, so avoid it by checking
                        # data type here
                        genelist.add(link_id)
                        if geneMapping[key] == 'NA' and type(link_id) is str:
                            geneMapping[key] = [link_id]
                        elif type(geneMapping[key]) is list and type(link_id) is str:
                            if link_id not in geneMapping[key]:
                                geneMapping[key].append(link_id)
                self.geneInfoDictionary[geneid] = geneMapping
        self.pathwaysWithGenesDictionary[this_pathway] = list(genelist)
        #print('At pathway {}, total {} genes, {} pathways with genes'\
         #     .format(this_pathway,len(self.geneInfoDictionary),len(self.pathwaysWithGenesDictionary)))
        
    def getMetabolitesIDFromGraph(self,g,this_pathway):
        type_predicate = URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
        
        metabolite_object = URIRef('http://vocabularies.wikipathways.org/wp#Metabolite')
        possible_source = {
            'kegg.compound':'kegg_id',
            'hmdb':'hmdb_id',
            'pubchem.compound':'pubchem_compound_id',
            'chebi':'chebi_id',
            'chemspider':'chemspider_id',
            'wikidata':'WikiData',
            'cas':'CAS',
            'lipidmaps':'LIPIDMAPS'
            }
        metabolite_list = set()
        for metabolites in g.subjects(type_predicate,metabolite_object):
            source = metabolites.split('/')
            source = source[len(source) - 2]
            # last items by split is retrieved
            metabolites_id = metabolites.split('/')[-1]
            #metabolites_id = self.getIDFromGraphLinks(g, metabolites)
            metabolites_id = self.prependID(source, metabolites_id)
            # predicate in RDF is defined here, this the subject/object with these predicates are extracted.
            id_mapping = {
                'chebi_id':'http://vocabularies.wikipathways.org/wp#bdbChEBI',
                'hmdb_id':'http://vocabularies.wikipathways.org/wp#bdbHmdb',
                'pubchem_compound_id':'http://vocabularies.wikipathways.org/wp#bdbPubChem',
                'WikiData':'http://vocabularies.wikipathways.org/wp#bdbWikidata',
                'chemspider_id':'http://vocabularies.wikipathways.org/wp#bdbChemspider'
            }
            # metabolites id mapping created each loop
            metaboliteMapping = {
                              "chebi_id": "NA", 
                              "drugbank_id": "NA", 
                              "drugbank_metabolite_id": "NA", 
                              "phenol_explorer_compound_id": "NA", 
                              "phenol_explorer_metabolite_id": "NA", 
                              "foodb_id": "NA", 
                              "knapsack_id": "NA", 
                              "chemspider_id": "NA",
                              "kegg_id": "NA",
                              "biocyc_id": "NA",
                              "bigg_id": "NA",
                              "wikipedia": "NA",
                              "nugowiki": "NA",
                              "metagene": "NA",
                              "metlin_id": "NA",
                              "pubchem_compound_id": "NA",
                              "het_id": "NA",
                              "hmdb_id": "NA",
                              "CAS": "NA",
                              'LIPIDMAPS':'NA',
                              'WikiData':'NA'}
            # skip pubchem.substance id at this moment
            #ttd.drug is new addition for the feb 10 2019 data
            if source not in ['pubchem.substance','drugbank','chembl.compound','kegg.drug', 'ttd.drug']:
                metaboliteMapping[possible_source[source]] = [metabolites_id]
                metabolite_list.add(metabolites_id)
                for key,value in id_mapping.items():
                    for links in g.objects(metabolites,URIRef(value)):
                        link_id = links.split('/')
                        link_id = link_id[len(link_id) - 1]
                        link_id = self.prependID(key, link_id)
                        metabolite_list.add(link_id)
                        # add id to the metabolites id mapping 
                        if metaboliteMapping[key] == 'NA' and type(link_id) is str:
                            metaboliteMapping[key] = [link_id]
                        elif type(metaboliteMapping[key]) is list and type(link_id) is str:
                            if link_id not in metaboliteMapping[key]:
                                metaboliteMapping[key].append(link_id)    
                        #print('Root {} has been linked to {}'.format(metabolites_id,link_id))
                #print(metaboliteMapping)
                self.pathwaysWithMetabolitesDictionary[this_pathway] = list(metabolite_list)
                self.metaboliteIDDictionary[metabolites_id] = metaboliteMapping
        '''
        print('At pathway {}:{}, there are {} metabolites'\
              .format(this_pathway,self.pathwayDictionary[this_pathway],len(self.metaboliteIDDictionary)))
        print('{} pathway has at least one metabolites'.format(len(self.pathwaysWithMetabolitesDictionary)))
        '''
        #print('WIKI  Total metabolites in this version of pathways (roughly):{}'.format(len(self.metaboliteIDDictionary)))
    # helper functions 
    def getIDFromGraphLinks(self,g,subject):
        '''
        Given RDF graph and subject,
        this function only looks at the object with identifier(source id) predicate
        '''
        for label in g.objects(URIRef(subject),DCTERMS.identifier):
            return label
        return 'NA'
    def getCatalyzation(self,g,this_pathway):
        '''
        This function get all catalyzation relations between metabolites and genes from the graph
        - param RDFgraph g The graph given by parsing RDF file
        - param string this_pathway the pathway ID (WP_number) given from the file.
        '''
        for s,p,o in g.triples((None,None,None)):
            if 'Interaction' not in s:
                continue
            print('---------------------------')
            print('Subject: {} \nPredicate: {}\nObject: {}'.format(s,p,o))
            print('---------------------------')
            time.sleep(1)
        
    def prependID(self,prefix,id):
        '''
        This function needs to be changed if the id_mapping dict has been changed
        This function will prepend prefix to id. Modified id based on the given prefix.
        '''
        # change id based on condition
        # Please look raw rdf file to figure out the pattern
        if prefix == 'pubchem_compound_id' or prefix == 'pubchem.compound':
            id = id.replace('CID','')
            id = 'pubchem:'+id
        elif prefix == 'chebi_id' or prefix =='chebi':
            id = id.replace('CHEBI:','')
            id = 'chebi:' + id
        elif prefix == 'hmdb_id' or prefix == 'hmdb':
            if len(id) == 9:
                id = id.replace('HMDB','HMDB00')
            elif len(id) != 9  and len(id) != 11:
                print('#### Check if hmdb id {} is correct ####'.format(id))
                time.sleep(3)
                '''
            elif len(id) == 11:
                print('HMDB id is in primary accession {}'.format(id))
                '''
            id = 'hmdb:' +id
        elif prefix == 'WikiData' or prefix =='wikidata' or prefix == 'entity':
            id = prefix.lower()+':'+id
        elif prefix == 'chemspider_id':
            id = 'chemspider:'+id
        elif prefix == 'kegg.compound' or prefix == 'kegg':
            id = 'kegg:' +id
        elif prefix == 'cas':
            id = 'CAS:' + id
        elif prefix == 'lipidmaps':
            id = 'LIPIDMAPS:' + id
        elif prefix == 'ncbigene' or prefix == 'Entrez':
            id = 'entrez:' + id
        elif prefix == 'uniprot' or prefix == 'UniProt':
            id = 'uniprot:' + id
        elif prefix == 'ec-code' or prefix == 'Enzyme Nomenclature':
            id = 'EN:' +id
        elif prefix == 'ensembl' or prefix == 'Ensembl':
            id = 'ensembl:' + id
        elif prefix == 'ncbiprotein':
            id = 'ncbiprotein:' + id
            
            
        return id
               
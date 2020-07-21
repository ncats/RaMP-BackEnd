import urllib.request as re
import time
import xml.etree.ElementTree as ET

class metaCycData():
    def __init__(self):
        # key name of organism value: number of duplicates
        self.organism = dict()
        # key name of pathway value: id of pathways
        self.pathways = dict()
        # key name of metabolites value: id of metabolites
        self.metabolites = dict()
        # key name of proteins value: id of proteins
        self.genes = dict()
    '''
    Go through biopax files to find all organisms 
    '''    
    def getDatabaseFile_level2_Organism(self):
        
        tree = ET.parse("../misc/data/MetaCyc/meta/21.1/data/biopax-level2.owl")
        root = tree.getroot()
        for child in root:
            tagtext = self.removeNameSpace(child.tag,'2')
            text = self.removeNameSpace(child.text,'2')
            #print(child.tag +"|"+ child.text)
            #print("After clearning ..." + tagtext)
            #for child2 in child.iter()
            for organism in child.iter("{http://www.biopax.org/release/biopax-level2.owl#}ORGANISM"):
                #print("Found ...")
                #print(self.removeNameSpace(organism.tag,'2'))
                for name in organism.iter("{http://www.biopax.org/release/biopax-level2.owl#}NAME"):
                    #print("Have name !")
                    nametext = name.text
                    #print("The organism's name is "+ nametext)
                    if nametext not in self.organism:
                        self.organism[nametext] = 1
                    else:
                        self.organism[nametext] +=1
                
            
            
            
    
    '''
    Replace namespace before each tags.
    '''
    def removeNameSpace(self,str,level2Or3):
        if level2Or3 == '2':
            namespace = {"{http://biocyc.org/biopax/biopax-level2#}",
                         "{http://www.w3.org/2002/07/owl#}",
                         "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}",
                         "{http://www.w3.org/2000/01/rdf-schema#}",
                         "{http://www.w3.org/2001/XMLSchema#}",
                         "{http://www.biopax.org/release/biopax-level2.owl#}"}
        elif level2Or3 == '3':
            namespace = {"{http://biocyc.org/biopax/biopax-level3#}",
                         "{http://www.w3.org/2002/07/owl#}",
                         "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}",
                         "{http://www.w3.org/2000/01/rdf-schema#}",
                         "{http://www.w3.org/2001/XMLSchema#}",
                         "{http://www.biopax.org/release/biopax-level3.owl#}"}
        for item in namespace:
            if item in str:
                str = str.replace(item,"")
                #print("cleaning ...")
                return(str)
        
        
        return(str)
    
    def getDatabaseFiles_level3(self):
        tree = ET.parse("../misc/data/MetaCyc/meta/21.1/data/biopax-level3.owl")
        root = tree.getroot()
        for child in root:
            tagtext = self.removeNameSpace(child.tag,'3')
            text = self.removeNameSpace(child.text,'3')
            #print(child.tag +"|"+ child.text)
            #print("After clearning ..." + tagtext)
            #time.sleep(3)
            #for child2 in child.iter()
            for organism in child.iter("{http://www.biopax.org/release/biopax-level3.owl#}organism"):
                #print("Found ...")
                #print(self.removeNameSpace(organism.tag,'3'))
                for name in organism.iter("{http://www.biopax.org/release/biopax-level3.owl#}name"):
                    #print("Have name !")
                    nametext = name.text
                    #print("The organism's name is "+ nametext)
                    if nametext not in self.organism:
                        self.organism[nametext] = 1
                    else:
                        self.organism[nametext] +=1    
    def outputToFile(self):
        organismName = open('../misc/output/metaCycOrganismCollection.text','wb')
        for key in self.organism:
            number = self.organism[key]
            organismName.write(key.encode('utf-8') +b'\t'+str(number).encode('utf-8') +b'\n')
            
        #print("Total organisms: "+ str(len(self.organism)))
        
    '''
    Find all pathways then add them to dictionary
    self.pathways
    '''
    def getPathways(self):
        tree = ET.parse("../misc/data/MetaCyc/meta/21.1/data/biopax-level2.owl")
        root = tree.getroot()
        
        ns = {"bp":"http://www.biopax.org/release/biopax-level2.owl#"}
        
        for child in root.findall('bp:pathway',ns):
            attrib = child.attrib
            id = "placeholder"
            for key in attrib:
                if "ID" in key:
                    id = attrib[key]
                    
            
            child2 = child.find('bp:NAME',ns)
            pathway = child2.text
            self.pathways[pathway] = id
            
    def getMetabolites(self):
        tree = ET.parse("../misc/data/MetaCyc/meta/21.1/data/biopax-level2.owl")
        root = tree.getroot()
        
        ns = {"bp":"http://www.biopax.org/release/biopax-level2.owl#"} 
        for child in root.findall('bp:physicalEntityParticipant',ns):
            #print(child.tag)
            physical = child.find('bp:PHYSICAL-ENTITY',ns)
            metabolite = physical.find('bp:smallMolecule',ns)
            protein = physical.find('bp:protein',ns)
            if metabolite is not None:
                name = metabolite.find('bp:NAME',ns)
                nametext = name.text
                attrib = metabolite.attrib
                id = "placeholder"
                for key in attrib:
                    if "ID" in key:
                        id = attrib[key]
                
                self.metabolites[nametext] = id
            
            if protein is not None:
                name = protein.find('bp:NAME',ns)
                nametext = name.text
                attrib = protein.attrib
                id = "placeholder"
                for key in attrib:
                    if "ID" in key:
                        id = attrib[key]
                
                self.genes[nametext] = id
        output1 = open('../misc/output/metaCycMetabolites.txt','wb')
        output2 = open('../misc/output/metaCycGenes.txt','wb')
        for key in self.metabolites:
            output1.write(key.encode('utf-8') +b'\t' + self.metabolites[key].encode('utf-8') +b'\n')
        
        for key in self.genes:
            output2.write(key.encode('utf-8') +b'\t' + self.genes[key].encode('utf-8') +b'\n')
            
        output1.close()
        output2.close()
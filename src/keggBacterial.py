import urllib.request as re
import time
import urllib.error 
from os import listdir
from multiprocessing import Pool
from MetabolomicsData import MetabolomicsData
import pandas as pd
class keggBacterial(MetabolomicsData):
    def __init__(self):
        # key: sub categories value: super categories
        self.SpeciesToGenus = dict()
        self.subtypeToSpecies = dict()
         
        # key: genus value: set of all pathways
        self.genusPathways = dict()
        # list of all pathway url
        self.pathwayurl = []
        # list of all compound
        self.metabolites = dict()
        # dict of pathway
        # key pathway id value pathway name without speicies
        self.pathwayDict = dict()
        self.metabolitesId = dict()
        # kye: species T ID value: Species categories
        self.keggIDToSpecies = dict()
        # total_pathways set for store all pathways
        self.total_pathway = set()
        # species pathway information
        self.keggSpeciesToPathways = dict()
    def findInt(self,str):
        try:
            str = str[str.find("(")+1:str.find(")")]
            num = int(str)
            return(num)
        except ValueError:
            return(0)
    # Total Bacterial 4387
    
    '''
    Must call in specific order 
    1. Download file with a list of 3-letters-id of bacterial
    2. Download files with name of 3-letters-id.txt from KEGG
    3. Download individual pathway id file 
    '''
    def getAllSpecies(self):
        '''
        Get all species from Kegg API.
        '''
        url = 'http://rest.kegg.jp/list/organism'
        dir = '../misc/data/kegg/keggBacterial/'
        self.check_path(dir)
        self.download_files(url, dir, file = 'allKeggSpecies.txt')
    def parseAllSpecies(self):
        '''
        Parse the species file 
        - return species {'Fungi', 'Animals', 'Bacteria', 'Plants', 'Protists', 'Archaea'}
        '''
        with open('../misc/data/kegg/keggBacterial/allKeggSpecies.txt') as f:
            content = [line.replace('\n','') for line in f if line is not '\n']
        tpOfCategories = (set(),set(),set(),set())
        for each in content:
            #print(each)
            #time.sleep(3)
            value = each.split('\t')
            categories = value[len(value) - 1].split(';')
            #print(categories)
            for i in range(0,len(categories)):
                #print(categories[i])
                #time.sleep(1)
                tpOfCategories[i].add(categories[i])
        
    def getSpeciesPathwayFile(self):
        '''
        For now, the function will collect all microbial that falls in below categories:
             {'Fungi','Bacteria','Protists', 'Archaea'}
        '''       
        target = {'Fungi','Bacteria','Protists', 'Archaea'}
        with open('../misc/data/kegg/keggBacterial/allKeggSpecies.txt') as f:
            content = [line.replace('\n','') for line in f if line is not '\n']
            for each in content:
                value = each.split('\t')
                categories = value[len(value) - 1].split(';')[1]
                if categories in target:
                    self.keggIDToSpecies[value[0]] = categories
        
        print('Total {} microbial is collected from kegg specis file.'.format(len(self.keggIDToSpecies)))
        dir = '../misc/data/kegg/keggBacterialPathway/'
        url = 'http://rest.kegg.jp/list/pathway/'
        self.check_path(dir)
        urls = [url + key for key in self.keggIDToSpecies]
        file_names = [key + '.txt' for key in self.keggIDToSpecies]
        dirs = [dir] * len(file_names)
        with Pool(30) as p:
            p.starmap(self.download_files,zip(urls,dirs,file_names))
    def getAllPathwaysFromMicrobial(self):
        '''
        After getting all files contains pathway information
        Try to find all unique pathways by looking at number part
        '''
        dir = '../misc/data/kegg/keggBacterialPathway/'
        microbial_files = listdir(dir)
        stored = set()
        for each in microbial_files:
            with open(dir+each) as f:
                content = [line.replace('\n','') for line in f if line is not '\n']
                pathways = {p.split('\t')[0][-5:] for p in content}
                self.total_pathway = self.total_pathway.union(pathways)
                species = self.keggIDToSpecies[each.replace('.txt','')]
                if species not in self.keggSpeciesToPathways:
                    self.keggSpeciesToPathways[species] = set()
                    self.keggSpeciesToPathways[species] = self.keggSpeciesToPathways[species].union(set(pathways))
                else:
                    self.keggSpeciesToPathways[species] = self.keggSpeciesToPathways[species].union(set(pathways))
        print('Total {} unique pathways are found associated to microbial.'.format(len(self.total_pathway)))
    def getAllLinkedCpdFromPathways(self):
        url = 'http://rest.kegg.jp/link/cpd/'
        dir_to_save = '../misc/data/kegg/keggBacterialCpd/'
        self.check_path(dir_to_save)
        file_names = ['map' + each + '.txt' for each in self.total_pathway]
        urls = [url +'map'+each for each in self.total_pathway]
        dir_to_save = [dir_to_save] * len(file_names)
        with Pool(3) as p:
            p.starmap(self.download_files,zip(urls,dir_to_save,file_names))
    def getAllCpdFromMapFile(self):
        dir_to_save = '../misc/data/kegg/keggBacterialCpd/'
        files_name = listdir(dir_to_save)
        setOfAllCpd = set()
        for each in files_name:
            with open(dir_to_save + each) as f:
                content = [line.replace('\n','') for line in f if line is not '\n']
                content = [line.split('\t')[1].replace('cpd:','') for line in content]
                setOfAllCpd = setOfAllCpd.union(set(content))
        df = pd.DataFrame({'cpd':list(setOfAllCpd)})
        df.to_csv('../misc/data/kegg/keggBacterialcpd.csv')
    def getDatabaseFile(self,haveFiles = True):
        pathwayFile = open('../misc/data/kegg/bacterial.txt')

        if not haveFiles:
            total = len(self.CategoiresID)
            for key in self.CategoiresID:
                pathToFile = "../misc/data/kegg/Bacterial/"
                filename = listdir(pathToFile)
                keyname = key+".txt"
                start = len(filename)
                if keyname not in filename: 
                    print("Downloading ... " + key +"\t" + str(start) +"/"+str(total))
                    start = start + 1
                    url = "http://rest.kegg.jp/list/pathway/" + key +"/"
                    try:
                        re.urlretrieve(url, "../misc/data/kegg/Bacterial/"+key+".txt") 
                    except urllib.error.HTTPError:
                        print(key + "Not Found ... ")
                        self.CategoiresIDNotFound.append(key)
                    except urllib.error.URLError:
                        print("Connection failed ...")
                    except FileNotFoundError:
                        #re.urlretrieve(url,"./misc/data/kegg/Bacterial/"+key+".txt")
                        print("Download failed for " + key)
    
    '''
    This function get all detail for each pathways based on given 3-letter-id number .txt file
    downloaded from last step.
    '''                
    def getAllPathwaysFiles(self,tot = 0):
        bacterialList = listdir("../misc/data/kegg/Bacterial")
        print("1.Total files need to be downloaded "+ str(len(bacterialList)))
        url = "http://rest.kegg.jp/get/"
        filelist = listdir("../misc/data/kegg/BacterialPathway/")
        i = 0
        for each in bacterialList:
            i = i + 1
            print("Downloading pathway in file " + str(i)+'/'+str(len(bacterialList)))
            pathwayFile = open("../misc/data/kegg/Bacterial/" + each)
            for line in pathwayFile:
                splitline = line.split("\t")
                filename = splitline[0] +".txt"
                
                filename = filename.replace("path:","")
                if filename not in filelist:
                    print("download " + filename)
                
                    re.urlretrieve(url+splitline[0],'../misc/data/kegg/BacterialPathway/'+ splitline[0].replace("path:","") + '.txt')  
            
            pathwayFile.close()      
    '''
    This function and the following functions deal with three-letters-id.txt file
    to extract all pathways to store them in dictionary.
    '''
    def parsePathwayAndProteins(self):
        pathwayFiles = "../misc/data/kegg/Bacterial/"
        filenames = listdir(pathwayFiles)
        tot = len(filenames)
        
        i = 0
        for file in filenames:
            
            print("Parsing ..."+file +":"+str(i)+"/"+str(tot))  
            i +=1
            dir = pathwayFiles + file
            self.parsePathwayFiles(dir,file.replace('.txt',''))      
                
    
    def parsePathwayFiles(self,dir,id):
        file = open(dir)
        for line in file:
            linesplit = line.split("\t")
            # print(linesplit)
            pathid = linesplit[0]
            pathid = pathid.replace('path:','')
            pathwayname = linesplit[1][:linesplit[1].find('-')-1]
            genusSpecies = linesplit[1][linesplit[1].find('-')+2:]
            genus = genusSpecies.split(" ")[0]
            species = genusSpecies.split(" ")[1]
            if(species not in self.SpeciesToGenus):
                self.SpeciesToGenus[species] = genus
            prefix = pathid[:pathid.find('0')]
            if prefix not in self.subtypeToSpecies:
                self.subtypeToSpecies[prefix] = species
            
            if genus not in self.genusPathways:
                self.genusPathways[genus] = set()
                self.genusPathways[genus].add(pathwayname)
            else:
                if pathwayname not in self.genusPathways[genus]:
                    self.genusPathways[genus].add(pathwayname)
                    '''
                else:
                    print("Duplicate:" + pathwayname)
    
                '''
    def outputToFile(self):
        keggBacterialFile = open('../misc/output/KeggBacterialPathNumbers.txt',"wb")
        keggGenusSpeciesID = open('../misc/output/KeggBacteralGenSpeId.txt','wb')
        for key in self.genusPathways:
            keggBacterialFile.write(key.encode('utf-8') +b"\t" +str(len(self.genusPathways[key])).encode('utf-8')+b"\n")
        
        for key in self.subtypeToSpecies:
            species = self.subtypeToSpecies[key]
            genus = self.SpeciesToGenus[species]
            
            keggGenusSpeciesID.write(genus.encode('utf-8') +b"\t" + species.encode('utf-8') +b'\t' +key.encode('utf-8')+b'\n')
        keggBacterialFile.close()
        keggGenusSpeciesID.close()     
        
    def urlOpenPathways(self):
        bacterialList = listdir("../misc/data/kegg/Bacterial")
        print("1.Total files need to be downloaded "+ str(len(bacterialList)))
        url = "http://rest.kegg.jp/get/"
        filelist = listdir("../misc/data/kegg/BacterialPathway/")
        i = 0
        for each in bacterialList:
            i = i + 1
            print("Downloading pathway in file " + str(i)+'/'+str(len(bacterialList)))
            pathwayFile = open("../misc/data/kegg/Bacterial/" + each)
            for line in pathwayFile:
                splitline = line.split("\t")
                filename = splitline[0] +".txt"
                if filename not in filelist:
                    filename = filename.replace("path:","")
                    self.pathwayurl.append(url+splitline[0])
            
            pathwayFile.close()
            #if i == 10:
            #    break
    '''
    Try to download all files, so slow even with multiprocessing ...
    1) self.urlopenpathways  to fill url in pathwayurl list
    2) self.download apply multiprocessing to download this each pathway files
    '''                
    def download(self,url):
        files = listdir("../misc/data/kegg/BacterialPathway/")
        
        print("Download in " + str(len(files))+"/"+str(len(self.pathwayurl)))
        splitline = url.split("/")
        id = splitline[len(splitline) - 1]
        id = id.replace("path:","")
        filename = id+".txt"
        if filename not in files:
            print("download "+ splitline[len(splitline) - 1])
            re.urlretrieve(url,"../misc/data/kegg/BacterialPathway/"+id +".txt")
        

    '''
    Get the number part of each pathway ids
    Then, the id can be used to map with metabolites or genes.
    '''    
    def getPathwayIds(self):
        dir = "../misc/data/kegg/Bacterial"
        pathwayFiles = listdir(dir)
        for each in pathwayFiles:
            file = open(dir +'/'+ each)
            for line in file:
                splitline = line.split('\t')
                idnum = splitline[0]
                idnum = idnum[idnum.find('0'):]
                pathway = splitline[1].split('-')[0]
                pathway = pathway[0:len(pathway) - 1]
                if idnum not in self.pathwayDict:
                    self.pathwayDict[idnum] = pathway
                    
    '''
    self.getmetabolitesFiles download files of mapping pathways with metabolites
    self.getMetabolitesId count how many metabolites id I have ...
    '''
    def getMetabolitesFiles(self):
        for pathwayid in self.pathwayDict:
            url = "http://rest.kegg.jp/link/cpd/map"
            url = url + pathwayid
            print(url)
            dir = '../misc/data/kegg/BacterialCompound/map' +pathwayid+".txt"
            re.urlretrieve(url,dir)            
    
    def getMetabolitesId(self):  
        dir = '../misc/data/kegg/BacterialCompound/'
        files = listdir(dir)
        for name in files:
            file = open(dir + name)
            for line in file:
                splitline = line.split('\t')
                if splitline[0] == '\n':
                    continue
                if len(splitline)>0:
                    cid = splitline[1]
                else:
                    pass
                
                cid = cid.replace('cpd:','')
                if cid not in self.metabolitesId:
                    self.metabolitesId[cid] = "Placeholder"
                
    '''
    get pathways mapping with genes KO ID for further counting.
    '''
    def getGeneMapFiles(self):
        files = listdir('../misc/data/kegg/BacterialPathway/')
        for each in self.pathwayDict:
            print('downloading ... ' + each)
            url = "http://rest.kegg.jp/link/ko/map"
            re.urlretrieve(url+each,"../misc/data/kegg/BacterialGeneLink/map" +each+".txt")
    
    def countGeneKoId(self):
        dir = '../misc/data/kegg/BacterialGeneLink/'
        koid = set()
        files = listdir(dir)
        for each in files:
            geneFile = open(dir+each)
            for line in geneFile:
                splitline = line.split('\t')
                print(splitline)
                id = splitline[1]
                id = id.replace('ko:','')
                if id not in koid:
                    koid.add(id)
            geneFile.close()
                
        
        return koid
                    
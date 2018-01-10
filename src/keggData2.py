import urllib.request as re
from os import listdir
import time

class keggData2():
    def __init__(self):
        # url list
        self.pathwayurls = []
        self.metaboliteurls = []
        self.geneurls = []
        self.genesDict = dict()
        self.metabolitesDict = dict()
    
    
    '''
    Get database files
        '''
    def getDatabaseFiles(self):
        re.urlretrieve("http://rest.kegg.jp/list/pathway/hsa", "../misc/data/kegg/hsa.txt")
        
        pathwayFilesOpen = open('../misc/data/kegg/hsa.txt')
        
        for line in pathwayFilesOpen:
            print(line)
            splitline = line.split('\t')
            pathwayid = splitline[0]
            
            url ='http://rest.kegg.jp/get/'
            self.pathwayurls.append(url + pathwayid)
    '''
    Download individual pathway files
    call after self.getDatabaseFiles
    '''
    def downloadPathways(self,urls):
        re.urlretrieve(url,'../misc/data/kegg/pathways')  
        
    
    def findGenesAndMetabolites(self):  
        files = listdir('../misc/data/kegg/pathways')
        for each in files:
            
            with open('../misc/data/kegg/pathways/' +each) as f:
                findGene = False
                for line in f:
                    leading = len(line) - len(line.lstrip())
                    # find first leading space == 0 and Gene in that
                    if leading ==0 and 'GENE' in line:
                        splitline = line.split('     ')
                        
                        findGene = True
                        gene = splitline[1].lstrip()
                        splitgene = gene.split('  ')
                        if splitgene[0] not in self.genesDict:
                            self.genesDict[splitgene[0]] = splitgene[1]
                        else:
                            print('find repeated ' + gene)
                            #time.sleep(1)
                        continue
                    if leading ==0 and 'Gene' not in line:
                        findGene = False
                    if findGene:
                        gene = line.lstrip()
                        splitgene = gene.split('  ')
                        if splitgene[0] not in self.genesDict:
                            self.genesDict[splitgene[0]] = splitgene[1]
                        else:
                            print('find repeated ' + gene)
                            #time.sleep(1)
                            
        for each in files:                    
            with open('../misc/data/kegg/pathways/' +each) as f:
                    findCompound = False
                    for line in f:
                        leading = len(line) - len(line.lstrip())
                        # find first leading space == 0 and Gene in that
                        if leading ==0 and 'COMPOUND' in line:
                            splitline = line.split('    ')
                            #print(splitline)
                            findCompound = True
                            compound = splitline[1].lstrip()
                            splitcompound = compound.split('  ')
                            if splitcompound[0] not in self.metabolitesDict:
                                #print(splitcompound)
                                #time.sleep(1)
                                self.metabolitesDict[splitcompound[0]] = splitcompound[1]
                            else:
                                print('find repeated ' + compound)
                                #time.sleep(1)
                            continue
                        if leading ==0 and 'COMPOUND' not in line:
                            findCompound = False
                        if findCompound:
                            compound = line.lstrip()
                            splitcompound = compound.split('  ')
                            if splitcompound[0] not in self.metabolitesDict:
                                if len(splitcompound) > 1:
                                    self.metabolitesDict[splitcompound[0]] = splitcompound[1]
                                #print(splitcompound)
                                #time.sleep(1)
                            else:
                                print('find repeated ' + compound)
                                #time.sleep(1)
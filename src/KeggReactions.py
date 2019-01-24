from KeggData import KeggData
import urllib.request as RE
import time
import os
from multiprocessing.dummy import Pool
import re
class keggReactions(KeggData):
    '''
    This class contains all information of kegg reaction 
    It has all method inherited from KeggData class
    The basic workflow for this class is:
    1) get all kegg human pahtways by self.getPathways() method
    2) Query API to find all linked reaction
    3) Query API to find detail for each reactions.
    '''
    def __init__(self):
        # super class is KeggData
        super().__init__()
        # key: pathway kegg ID value: Reactions linked to the pathways
        self.pathwayWithReactions = dict()
        self.Reactions = dict()
        
    def downloadPathwaysWithReactions(self,speed_up = 5):
        '''
        This function download all pathways that has reactions into 
        '../misc/kegg/pathwaysWithReactions/'
        The url to query API is 'http://rest.kegg.jp/link/rn/[pathway map id]'
        param int speed_up This parameter speed up the download process by calling multiprocess
        
        '''
        assert len(self.pathwayDictionary) != 0,\
        'You should first call self.getPathways to fill in the pathway dictionary'
        url = 'http://rest.kegg.jp/link/rn/map'
        path = '../misc/data/kegg/PathwayWithReactions/'
        # initialize all urls and files name for downloading
        # paths will be same for each urls and files
        urls,paths,files = [[],[],[]]
        self.check_path(path)
        for key in self.pathwayDictionary:
            urls.append(url+key)
            paths.append(path)
            files.append('map'+key+'rn.txt')
        #print('{} urls, {} in the paths, {} files to download'
        #      .format(len(urls),len(paths),len(files)))
        with Pool(speed_up) as p:
            p.starmap(self.download_files,zip(urls,paths,files))
        
        
    def downloadReactionsInfo(self,speed_up = 5):
        '''
        Download information for each individual reaction
        Store them at ../misc/data/kegg/Reactions/
        param int speed_up Speed download process up
        '''
        assert len(self.pathwayDictionary) != 0,\
        'You should first call self.getPathways to fill in the pathway dictionary'
        print('Start getting reactions ...')
        path = '../misc/data/kegg/PathwayWithReactions/'
        reactionsFile = os.listdir(path)
        setOfReactions = []
        setOfPathwaysWithReactions = set()
        # get all reactions from files
        for each in reactionsFile:
            with open(path+each) as f:
                reactions = [line.rstrip('\n') for line in f if line != '\n']
                reactions = list(map(lambda x:x[x.find('rn:')+3:],reactions))
                if not len(reactions) == 0 and not reactions == ['']:
                    setOfPathwaysWithReactions.add(each)
                setOfReactions.extend(reactions)
        setOfReactions = list(set(setOfReactions))
        
        print('{} pathways have reactions, total reactions are {}'
              .format(len(setOfPathwaysWithReactions),len(setOfReactions)))
        urls = list(map(lambda x:'http://rest.kegg.jp/get/'+x,setOfReactions))
        paths = ['../misc/data/kegg/Reactions/']*len(urls)
        file_names = list(map(lambda x:x+'.txt',setOfReactions))
        self.check_path(paths[0])
        with Pool(speed_up) as p:
            p.starmap(self.download_files,zip(urls,paths,file_names))
    
    def getPathwaysWithReactions(self):
        print('Get pathways and reactions ...')
        assert len(self.pathwayDictionary) != 0,\
        'You should first call self.getPathways to fill in the pathway dictionary'
        url = 'http://rest.kegg.jp/link/rn/map'
        path = '../misc/data/kegg/PathwayWithReactions/'
        # initialize all urls and files name for downloading
        # paths will be same for each urls and files
        urls,paths,files = [[],[],[]]
        self.check_path(path)
        for key in self.pathwayDictionary:
            with open(path+'map'+key+'rn.txt') as f:
                reactions = [line.rstrip('\n') for line in f if line != '\n']
                reactions = list(map(lambda x:x[x.find('rn:') + 3:],reactions))
                if len(reactions) > 0:
                    self.pathwayWithReactions[key] = reactions
        
        print('{} pathways has reaction'.format(len(self.pathwayWithReactions)))
        
    def getReactionsInfo(self):
        paths = '../misc/data/kegg/Reactions/'
        reactions = os.listdir(paths)
        # regular expression for the kegg compound id
        # C+ number part
        pat = re.compile(r'(C\d+)')
        
        for react in reactions:
            with open(paths + react,'r') as f:
                reaction = {'substrates':[],
                            'products':[]}
                rnid = react.replace('.txt','')
                #print('Reaction {}'.format(rnid))
                for line in f:
                    # get rid of new line character
                    line = line.rstrip('\n')
                    if 'NAME' in line:
                        #print(line)
                        time.sleep(1)
                    elif 'EQUATION' in line:
                        # split to substrates and products
                        eqn = line.split('<=>')
                        substrates = pat.findall(eqn[0])
                        products = pat.findall(eqn[1])
                        reaction['substrates'] = substrates
                        reaction['products'] = products
                        
                        #print("{} = {}".format(substrates,products))
                        
                        time.sleep(1)
        


if __name__ =="__main__":
    kegg = keggReactions()
    kegg.getPathways()
    kegg.getPathwaysWithReactions()
    print(kegg.pathwayWithReactions['00010'])
    kegg.getReactionsInfo()
        
        
    
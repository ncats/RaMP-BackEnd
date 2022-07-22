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
                
        self.isTransport = 0        
                
        # compounds ids for left side, populate using rhea2ec tsv file.
        self.left_comp_ids = []
        
        # compound ids for right side
        self.right_comp_ids = []
        
        # compounds ids for left side, populate using rhea2ec tsv file.
        self.left_comps = dict()
        
        # compound ids for right side
        self.right_comps = dict()
                
    def getBasicRecordString(self):
        ec = self.ec
        if ec is None:
           ec = "" 
        dir = self.direction
        if dir is None:
            dir = ""
        
        s = self.rhea_id + "\t" + str(self.status) + "\t" + str(self.isTransport) + "\t" +self.direction + "\t" + self.rhea_label + "\t" + self.rhea_equation + "\t" + self.rhea_html_eq + "\t" + ec + "\n" 
       
        return s
        
    def getRheaIdToCompMappingString(self):
        s = ""
        for cid in self.left_comp_ids:
            s = s + self.rhea_id + "\t" + cid + "\t0\n" 
        for cid in self.right_comp_ids:
            s = s + self.rhea_id + "\t" + cid + "\t1\n"             
        return s
    
    def getRheaIdToUniprotMappingString(self):    
        s = ""
        for pid in self.proteins:
            s = s + self.rhea_id + "\t" + pid + "\t0\n"             
        return s

    def getCompoundToProteinString(self):
        s = ""
        for pid in self.proteins:
            for cid in self.left_comp_ids:
                s = s + self.rhea_id + "\t" + pid + "\t" + cid + "\t0\n"
            for cid in self.right_comp_ids:
                s = s + self.rhea_id + "\t" + pid + "\t" + cid + "\t1\n"
        return s
                

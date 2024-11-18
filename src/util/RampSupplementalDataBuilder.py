'''
Created on Aug 2, 2023

@author: braistedjc
'''
from sqlalchemy import create_engine
from src.util.SimilarityMatrix import SimilarityMatrix, SimilarityMatrix_metabolite, SimilarityMatrix_gene


class RampSupplementalDataBuilder(object):
    '''
    classdocs
    '''
    analytesMatrix = None
    metabolitesMatrix = None
    genesMatrix = None

    def __init__(self, sqliteCreds=None):
        '''
        Constructor
        '''
        
        # a MySQL RaMP db_properties file, or an SQLite DB file 
        self.credInfo = sqliteCreds

        self.conn = self.createSQLiteEngine(self.credInfo).connect()

    def getPathwaysWithSameAnalytes(self):
        if self.analytesMatrix is None:
            self.initialize_similarity_matrices()
        return self.analytesMatrix.getDuplicates()

    
    def createSQLiteEngine(self, sqliteFile=None):
        engine = create_engine('sqlite:///'+sqliteFile, echo=False)
        return engine

    def getMergedSimilarityMatrix(self):
        compress = True
        if self.analytesMatrix is None:
            self.initialize_similarity_matrices()

        metaboliteCounts = self.metabolitesMatrix.getCounts()
        geneCounts = self.genesMatrix.getCounts()

        analyteSparseTuples = self.analytesMatrix.getSerializedSparseTuples(compress)
        metaboliteSparseTuples = self.metabolitesMatrix.getSerializedSparseTuples(compress)
        geneSparseTuples = self.genesMatrix.getSerializedSparseTuples(compress)

        allKeys = sorted(set(metaboliteCounts) | set(geneCounts))

        merged_dict = {
            key: {
                "metabolite_count": metaboliteCounts.get(key),
                "gene_count": geneCounts.get(key),
                "analyte_blob": analyteSparseTuples.get(key),
                "metabolite_blob": metaboliteSparseTuples.get(key),
                "gene_blob": geneSparseTuples.get(key)
            }
            for key in allKeys
        }
        return merged_dict

    def initialize_similarity_matrices(self):
        self.analytesMatrix = SimilarityMatrix(connection=self.conn)
        self.metabolitesMatrix = SimilarityMatrix_metabolite(connection=self.conn)
        self.genesMatrix = SimilarityMatrix_gene(connection=self.conn)

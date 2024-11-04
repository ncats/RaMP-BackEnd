import zlib
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict

import numpy as np
from bitarray import bitarray
from sqlalchemy.engine import Connection


@dataclass
class SimilarityMatrix:
    connection: Connection
    analyte_filter: str = None
    include_SMPDB: bool = False
    include_PFOCR: bool = True
    analyte_count_cutoff: int = 10
    valid_pathways: List[str] = field(default_factory=list)
    analyte_index: Dict[str, int] = field(default_factory=dict)
    pathway_vectors: Dict[str, bitarray] = None


    def validPathwayList(self):
        if len(self.valid_pathways) == 0:
            rows = self._getValidPathways()
            self.valid_pathways = [row['pathwayRampId'] for row in rows]
        return self.valid_pathways

    def _getValidPathways(self):
        return self.connection.execute(self.allPathwaysQuery()).all()

    def getCounts(self):
        vector_dictionary = self._getPathwayVectors()
        count_dictionary = {}
        for pathway, vector in vector_dictionary.items():
            count_dictionary[pathway] = vector.count()
        return count_dictionary

    def getDuplicates(self):
        vector_dictionary = self._getPathwayVectors()

        hash_map = defaultdict(list)

        for pathway, bitarr in vector_dictionary.items():
            bit_hash = bitarr.tobytes()
            hash_map[bit_hash].append(pathway)

        duplicates = {k: v for k, v in hash_map.items() if len(v) > 1}

        pairs = []
        for bit_hash, pathway_list in duplicates.items():
            for index1, pathway1 in enumerate(pathway_list):
                for index2, pathway2 in enumerate(pathway_list):
                    if index2 <= index1:
                        continue
                    pairs.append((pathway1, pathway2))

        return pairs

    def getSimilarityMatrix(self):
        vector_dictionary = self._getPathwayVectors()
        pathways = self.validPathwayList()
        matrix = np.zeros(shape=(len(pathways), len(pathways)))
        for left_index, left in enumerate(pathways):
            for right_index, right in enumerate(pathways):
                if right_index <= left_index:
                    continue
                left_vector = vector_dictionary[left]
                right_vector = vector_dictionary[right]
                dist = self._calculate_jaccard(left_vector, right_vector)
                matrix[left_index, right_index] = dist
                matrix[right_index, left_index] = dist
        return matrix

    def getSparseSimilarityTuples(self):
        vector_dictionary = self._getPathwayVectors()
        pathways = self.validPathwayList()
        sparse_tuples = {}
        for left_index, left in enumerate(pathways):
            row_tuples = [[left_index, '-']]
            last_index = left_index
            for right_index, right in enumerate(pathways):
                if right_index <= left_index:
                    continue
                left_vector = vector_dictionary[left]
                right_vector = vector_dictionary[right]
                dist = self._calculate_jaccard(left_vector, right_vector)
                if dist > 0:
                    index_diff = right_index - last_index
                    row_tuples.append([index_diff, dist])
                    last_index = right_index
            sparse_tuples[left] = row_tuples
        return sparse_tuples

    def getSerializedSparseTuples(self, compress):
        sparse_tuples = self.getSparseSimilarityTuples()
        return {
            pathwayId: self._serialize_one_row_of_tuples(tuples, compress) for pathwayId, tuples in sparse_tuples.items()
        }

    def _serialize_one_row_of_tuples(self, simTuple, compress = False):
        row_pieces = []
        for pair in simTuple:
            pair_str = f"{pair[0]},{pair[1]}"
            row_pieces.append(pair_str)
        if len(row_pieces) == 0:
            return None
        if compress:
            return zlib.compress("|".join(row_pieces).encode())
        return "|".join(row_pieces)



    def _getPathwayVectors(self):
        if self.pathway_vectors is None:
            self.pathway_vectors = {}
            associations = self._getAssociations()
            analyteCount = self.getAnalyteCount()
            for row in associations:
                pathwayRampId = row[0]
                analyte_ids = row[1].split(',')
                vector = bitarray(analyteCount)
                vector.setall(0)
                self.pathway_vectors[pathwayRampId] = vector
                for analyte_id in analyte_ids:
                    vector[self.getAnalyteIndex(analyte_id)] = 1
        return self.pathway_vectors


    def _calculate_jaccard(self, bitarray1, bitarray2):
        intersection = bitarray1 & bitarray2
        intersection_count = intersection.count()
        if intersection_count == 0:
            return 0

        union = bitarray1 | bitarray2
        union_count = union.count()
        if union_count == 0:
            return 0

        return round(1000 * intersection_count / union_count)

    def _getAssociations(self):
        results = self.connection.execute(self.allAssociationsQuery()).all()
        return results


    def allPathwaysQuery(self):
        return f"""select distinct ap.pathwayRampId 
                    from analytehaspathway ap, pathway p
                    where {self._pathwayFilter()}
                    and ap.pathwayRampId = p.pathwayRampId
                    {self._analyteFilter()}
                    GROUP BY ap.pathwayRampId
                    HAVING count(distinct ap.rampId) >= {self.analyte_count_cutoff}
                    order by ap.pathwayRampId"""
    def allAnalytesQuery(self):
        return f"""select distinct rampID 
                    from analytehaspathway ap, pathway p 
                    where ap.pathwayRampId in ({self.quotedPathwayList()})
                    {self._analyteFilter()}
                    and ap.pathwayRampId = p.pathwayRampId
                    order by ap.rampID"""
    def allAssociationsQuery(self):
        return f"""select ap.pathwayRampId, group_concat(ap.rampID) 
                    from analytehaspathway ap, pathway p 
                    where ap.pathwayRampId in ({self.quotedPathwayList()})
                    {self._analyteFilter()}
                    and ap.pathwayRampId = p.pathwayRampId
                    group by ap.pathwayRampId
                    order by ap.pathwayRampId"""

    def _pathwayFilter(self):
        clauses = []
        if not self.include_SMPDB:
            clauses.append("p.type != 'hmdb'")
        if not self.include_PFOCR:
            clauses.append("p.type != 'pfocr'")
        return " AND ".join(clauses)

    def _analyteFilter(self):
        if self.analyte_filter is None:
            return ""
        return "AND ap.rampId like '" + self.analyte_filter + "'"




    def quotedPathwayList(self):
        pathways = self.validPathwayList()
        quotedPathwayRampIds = [f"'{pathway}'" for pathway in pathways]
        return ', '.join(quotedPathwayRampIds)

    def getAnalyteCount(self):
        analyte_index = self._getInternalAnalyteIndex()
        return len(analyte_index)

    def getAnalyteIndex(self, analyte_id):
        analyte_index = self._getInternalAnalyteIndex()
        return analyte_index[analyte_id]

    def _getInternalAnalyteIndex(self):
        if len(self.analyte_index) == 0:
            self._initializeAnalyteIndex()
        return self.analyte_index

    def _initializeAnalyteIndex(self):
        analyte_result = self.connection.execute(self.allAnalytesQuery()).all()
        self.analyte_index = {analyte[0]: idx for idx, analyte in enumerate(analyte_result)}


@dataclass
class SimilarityMatrix_metabolite(SimilarityMatrix):
    analyte_filter: str = "RAMP_C%"
    analyte_count_cutoff: int = 5


@dataclass
class SimilarityMatrix_gene(SimilarityMatrix):
    analyte_filter: str = "RAMP_G%"
    analyte_count_cutoff: int = 5

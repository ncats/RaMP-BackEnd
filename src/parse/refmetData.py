import csv
import os
from os.path import exists
from typing import List, Dict

from parse.MetabolomicsData import MetabolomicsData

class refmetData(MetabolomicsData):
    def __init__(self, resConfig):
        self.config = resConfig
        self.pathwayDictionary: Dict[str, str] = {}
        self.metabolitesWithPathwaysDictionary: Dict[str, List[str]] = {}
        self.metaboliteCommonName: Dict[str, str] = {}
        self.pathwaysWithGenesDictionary: Dict[str, List[str]] = {}
        self.pathwayCategory: Dict[str, str] = {}
        self.metaboliteIDDictionary: Dict[str, Dict] = {}
        self.geneInfoDictionary: Dict[str, Dict] = {}

    def getEverything(self):

        met_conf = self.config.getConfig("refmet_met")

        file = met_conf.sourceFileName
        url = met_conf.sourceURL
        local_directory = met_conf.localDir

        if not exists(local_directory):
            os.mkdir(local_directory)

        if not exists(local_directory + file):
            self.download_files(url, local_directory, file)

    def idHasCrossRefs(self, refmet_row):
        id_types = ['pubchem_cid', 'chebi_id', 'hmdb_id', 'lipidmaps_id', 'kegg_id']
        for id_field in id_types:
            if refmet_row[id_field].strip() != '':
                return True
        return False

    def processRefMet(self):
        met_conf = self.config.getConfig("refmet_met")
        local_directory = met_conf.localDir
        file = met_conf.sourceFileName

        with open(local_directory + file, mode='r') as csvfile:
            csvreader = csv.DictReader(csvfile, skipinitialspace=True)
            for row in csvreader:
                if not self.idHasCrossRefs(row):
                    continue
                refmet_id = row['refmet_id'].strip()

                id_types = ['pubchem', 'chebi', 'hmdb', 'LIPIDMAPS', 'kegg']
                column_names = ['pubchem_cid', 'chebi_id', 'hmdb_id', 'lipidmaps_id', 'kegg_id']
                id_mapping = {'refmet': "refmet:" + refmet_id}
                for id_type, column_name in zip(id_types, column_names):
                    equivalent_id = row[column_name].strip()
                    if equivalent_id != "":
                        id_mapping[id_type] = f"{id_type}:{equivalent_id}"
                self.metaboliteIDDictionary["refmet:" + refmet_id] = id_mapping
                self.metaboliteCommonName["refmet:" + refmet_id] = row['refmet_name'].strip()

        self.write_myself_files(database = "refmet")
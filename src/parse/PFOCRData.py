import os
from dataclasses import dataclass, field
from os.path import exists
from typing import List, Dict

from parse.MetabolomicsData import MetabolomicsData

class PFOCRData(MetabolomicsData):
    def __init__(self, resConfig):
        self.config = resConfig
        self.pathwayDictionary: Dict[str, str] = {}
        self.metabolitesWithPathwaysDictionary: Dict[str, List[str]] = {}
        self.pathwaysWithGenesDictionary: Dict[str, List[str]] = {}
        self.pathwayCategory: Dict[str, str] = {}


    def getEverything(self):

        met_conf = self.config.getConfig("pfocr_met")
        gene_conf = self.config.getConfig("pfocr_gene")

        for conf in [met_conf, gene_conf]:
            file = conf.sourceFileName
            url = conf.sourceURL
            local_directory = conf.localDir

            self.download_files(url, local_directory, file)
            self.checkFileFormat(local_directory + file)

    def checkFileFormat(self, file_path):
        with open(file_path, 'r') as gmt_file:
            lines = gmt_file.readlines()
            for index, line in enumerate(lines):
                if not line.startswith('PMC'):
                    for bad_line in lines[index-1:index+2]:
                        print(bad_line)
                    raise Exception(f'{file_path} needs fixing. Lines should start with PMC')

    def ensure_id_is_unique(self, id, title) -> str:
        if id not in self.pathwayDictionary:
            return id
        existing_title = self.pathwayDictionary[id]
        if existing_title == title:
            return id

        pieces = id.split('.')
        id_to_use = pieces[0]
        count = 0
        if len(pieces) > 1:
            count = int(pieces[1]) + 1

        return self.ensure_id_is_unique(f"{id_to_use}.{count}", title)

    def processPathways(self):
        for conf_string in ['pfocr_met', 'pfocr_gene']:
            conf = self.config.getConfig(conf_string)
            file = conf.sourceFileName
            local_directory = conf.localDir

            with open(local_directory + file, 'r') as gmt_file:
                for line in gmt_file:
                    parts = line.split("\t")
                    figure_id =  parts[0].strip()
                    title = parts[1].strip()
                    id_list = [id.strip() for id in parts[2:]]

                    figure_id = self.ensure_id_is_unique(figure_id, title)
                    if figure_id in self.pathwayDictionary:
                        old_title = self.pathwayDictionary[figure_id]
                        assert old_title == title

                    self.pathwayDictionary[figure_id] = title
                    self.pathwayCategory[figure_id] = "NA"

                    if conf_string.endswith('met'):
                        for id in id_list:
                            if id in self.metabolitesWithPathwaysDictionary:
                                self.metabolitesWithPathwaysDictionary[id].append(figure_id)
                            else:
                                self.metabolitesWithPathwaysDictionary[id] = [figure_id]

                    if conf_string.endswith('gene'):
                        if figure_id in self.pathwaysWithGenesDictionary:
                            raise Exception('i dont think this happens')
                        prefixed_id_list = [f"entrez:{id}" for id in id_list]
                        self.pathwaysWithGenesDictionary[figure_id] = prefixed_id_list

        self.write_myself_files('pfocr')


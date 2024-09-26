from random import Random
from typing import List

from src.rampConfig.RampConfig import RampConfig
from src.util.EntityBuilder import EntityBuilder

resourceConfFile = "../config/external_resource_config.txt"
resourceConf = RampConfig()
resourceConf.loadConfig(resourceConfFile)

builder = EntityBuilder(resourceConf)
builder.loadMetaboList()
builder.addMetaboliteCommonName()
builder.addMetaboliteSynonyms()

best_name_counts = {}

class NameList:
    names: List[str] = []
    uniqueness_threshold = 1 # must be unique

    name_count: dict  # Class-level dict
    def __init__(self):
        if not hasattr(NameList, 'name_count'):
            NameList.name_count = {}
        self.names = []

    def add_name(self, new_name):
        if new_name not in self.names:
            self.names.append(new_name)
            NameList.add_to_global_list(new_name)

    @staticmethod
    def add_to_global_list(name):
        if name in NameList.name_count:
            NameList.name_count[name] += 1
        else:
            NameList.name_count[name] = 1

    def get_name(self):
        if len(self.names) == 0:
            return None
        self.names.sort(key=len)
        most_unique_name = None
        best_ncompounds = None
        for name in self.names:
            ncompounds = NameList.name_count[name]
            if ncompounds <= self.uniqueness_threshold:
                return name
            if best_ncompounds is None or ncompounds < best_ncompounds:
                best_ncompounds = ncompounds
                most_unique_name = name

        print(f'there is no unique name in this list, but {most_unique_name} has {best_ncompounds} and that is best we can do')
        if most_unique_name is not None:
            return most_unique_name
        return self.names[0]

    @staticmethod
    def get_names_from_metabolite(met):
        name_list = NameList()
        for source, id_dict in met.commonNameDict.items():
            for id_source, val in id_dict.items():
                name_list.add_name(val)
        return name_list

    @staticmethod
    def get_names_and_synonyms_from_metabolite(met):
        name_list = NameList()
        for source, id_dict in met.commonNameDict.items():
            for id_source, val in id_dict.items():
                name_list.add_name(val)
        for source, ids in met.synonymDict.items():
            for val in ids:
                name_list.add_name(val)
        return name_list

    @staticmethod
    def get_names_by_id_type(met):
        name_dicts = {}
        source_info = met.getSourceData()
        for row in source_info:
            source = get_source_name(row.IDtype)
            name = row.commonName
            if name is None or name == 'None':
                continue
            if source in name_dicts:
                source_list = name_dicts[source]
            else:
                source_list = []
                name_dicts[source] = source_list
            if name not in source_list:
                source_list.append(name)

        for source in name_dicts.keys():
            name_dicts[source].sort(key=len)

        return name_dicts

def get_source_name(val):
    if val == 'rhea-comp':
        return 'rhea'
    if val == 'LIPIDMAPS':
        return 'lipidmaps'
    if val == 'wikidata':
        return 'wiki'
    return val

name_dict = {}
name_syn_dict = {}
common_name_dict = {}
source_name_dicts = {}
all_sources = []

count = 0

met_list = builder.metaboliteList.getUniqueMetabolites()
met_list = sorted(met_list, key=lambda x: x.rampId)

for met in met_list:
    # count += 1
    # if count > 10:
    #     continue
    name_dict[met.rampId] = NameList.get_names_from_metabolite(met)
    name_syn_dict[met.rampId] = NameList.get_names_and_synonyms_from_metabolite(met)
    common_name_dict[met.rampId] = met.get_most_common_value()
    source_name_dicts[met.rampId] = NameList.get_names_by_id_type(met)
    for source in source_name_dicts[met.rampId]:
        if source not in all_sources:
            all_sources.append(source)

all_sources.sort()

html_out = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>How can we find the best Metabolite name?</title>
    <style>
    
        .badge {
            border-radius: 50%;
            width: 20px;
            height: 20px;
            display: inline-block;
            margin-right: 5px;
            margin-left: 10px;
        }
        .badge-yellow {
            background-color: orange;
        }
        .badge-green {
            background-color: lightgreen;
        }
        .badge-blue {
            background-color: lightblue;
        }
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f4f4f9;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 20px;
        }

        th, td {
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
            column-width: 50em;
        }

        th {
            background-color: #333;
            color: white;
        }

        tr:nth-child(even) {
            background-color: #f2f2f2;
        }

        tr:nth-child(odd) {
            background-color: #ffffff;
        }

        tr:hover {
            background-color: #ddd;
        }

        .center {
            text-align: center;
        }
    </style>
</head>
<body>
    <h1>How can we find the best Metabolite name?</h1>
    
    <p> Review the random sample of metabolites below.</p> 
    <p>You will find all the metabolite names and synonyms, sorted by the length of the name.</p>
    <p>The number to the left of the name is the count of metabolites in RaMP that have that same name. (1) means the name is unique to one metabolite.</p> 

    <h3>Shortest Unique Name</h3>
    <p> This name is labeled with a <b class='badge badge-yellow'></b>. It is the shortest common name (from the 'source' table) that is unique in the database. If there are no unique names, it is the closest to being unique.</p>
    
    <h3>Shortest Unique Name or Synonym</h3>
    <p> This name is labeled with a <b class='badge badge-green'></b>. Considering both names and synonyms (both the 'source' and 'analytesynonym' tables). </p>

    <h3>Most common common name</h3>
    <p> This name is labeled with a <b class='badge badge-blue'></b>. It is the most common common name in the 'source' table across all data sources.</p>    
    
    <h3>Sources</h3>
    <p>The IDs divided by the data source</p>
    
    <table>
        <thead>
            <tr>
                <th>ramp ID</th>
                <th>All Names and Synonyms</th>
                <th>IDs by Source</th>
            </tr>
        </thead>
        <tbody>
            {{body_placeholder}}
        </tbody>
    </table>
</body>
</html>
"""

data_html = []

r = Random(1)
random_keys = sorted(name_dict.keys(), key=lambda x: r.random())[:1000]

def get_name_html(name, best_1, best_2, best_3):
    name_count = NameList.name_count[name]
    html = f'<b>({name_count}) - {name}</b>'
    if name == best_1:
        html = html + f'<b class="badge badge-yellow"></b>'
    if name == best_2:
        html = html + f'<b class="badge badge-green"></b>'
    if name == best_3:
        html = html + f'<b class="badge badge-blue"></b>'
    return html


for rampId in random_keys:
    name_list = name_dict[rampId]
    name_syn_list = name_syn_dict[rampId]
    best_name_or_synonym = name_syn_list.get_name()
    best_name = name_list.get_name()
    commonest_name = common_name_dict[rampId]
    source_dicts = source_name_dicts[rampId]
    source_data = []
    for source in all_sources:
        if source in source_dicts:
            table_html = f"""
            <table>
            <tr>
            <th>{source}</th>
            </tr>"""
            source_list = source_dicts[source]
            table_html = table_html + "".join(f"<tr><td>{get_name_html(id, best_name, best_name_or_synonym, commonest_name)}</td></tr>" for id in source_list)
            table_html = table_html + "</table>"
            source_data.append(table_html)

    row = f"""<tr>
    <td>{rampId}</td>
    <td>{"<br />".join([get_name_html(name, best_name, best_name_or_synonym, commonest_name) for name in name_syn_list.names])}</td>
    <td>{"".join(source_data)}</td>
    </tr>"""
    data_html.append(row)

with open("metabolite_name_report.html", 'w') as file:
    html_out = html_out.replace("{{body_placeholder}}", "".join(data_html))
    file.write(html_out)

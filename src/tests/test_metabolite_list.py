import sys

sys.path.append('../src')
import pytest
from unittest.mock import patch
from src.rampEntity.Metabolite import Metabolite
from src.rampEntity.MetaboliteList import MetaboliteList

no_merges = [
    {
        "name": "r1",
        "ids": ['idr1','idr1_1'],
        "keys": ["1", "2"]
    },
    {
        "name": "r2",
        "ids": ['idr2'],
        "keys": ["3", "4"]
    },
    {
        "name": "r3",
        "ids": ['idr3'],
        "keys": ["5", "6"]
    },
    {
        "name": "r4",
        "ids": ['idr4'],
        "keys": ["7", "8"]
    }
]

partial_merges = no_merges + [
    {
        "name": "r5",
        "ids": ['idr5'],
        "keys": ["1", "4"]
    }
]

all_merges = partial_merges + [
    {
        "name": "r6",
        "ids": ['idr6'],
        "keys": ["5", "8"]
    },
    {
        "name": "r7",
        "ids": ['idr7'],
        "keys": ["2", "7"]
    }
]


@pytest.fixture
def metabolite_list():
    return MetaboliteList()

def initialize_metabolite_list(metabolite_list, test_data: list):
    for entry in test_data:
        m = Metabolite()
        m.rampId = entry["name"]
        for id in entry["ids"]:
            m.addId(id, 'source')
            metabolite_list.addMetaboliteByAltId(id, m)

def mock_get_inchi_keys(self, metabolite):
    test_entry = next((item for item in all_merges if item['name'] == metabolite.rampId), None)
    return set(test_entry["keys"])


@patch.object(MetaboliteList, 'get_inchi_keys', mock_get_inchi_keys)
def test_multi_level_mergers(metabolite_list):
    initialize_metabolite_list(metabolite_list, all_merges)
    metabolite_list.collapseMetsOnInchiKeyPrefix()
    mets = metabolite_list.getUniqueMetabolites()
    assert len(mets) == 1

    valid_id_list = (['idr1', 'idr1_1', 'idr2', 'idr3', 'idr4', 'idr5', 'idr6', 'idr7'])
    valid_id_list.sort()
    valid_id_string = "-".join(valid_id_list)

    id_list = mets[0].idList
    id_list.sort()
    id_list_string = "-".join(id_list)

    assert valid_id_string == id_list_string

@patch.object(MetaboliteList, 'get_inchi_keys', mock_get_inchi_keys)
def test_no_mergers(metabolite_list):
    initialize_metabolite_list(metabolite_list, no_merges)
    metabolite_list.collapseMetsOnInchiKeyPrefix()
    mets = metabolite_list.getUniqueMetabolites()
    assert len(mets) == 4
    assert mets[0].rampId != mets[1].rampId
    assert mets[1].rampId != mets[2].rampId
    assert mets[2].rampId != mets[3].rampId
    assert mets[3].rampId != mets[0].rampId

@patch.object(MetaboliteList, 'get_inchi_keys', mock_get_inchi_keys)
def test_mixed_merges(metabolite_list):
    initialize_metabolite_list(metabolite_list, partial_merges)
    metabolite_list.collapseMetsOnInchiKeyPrefix()
    mets = metabolite_list.getUniqueMetabolites()
    assert len(mets) == 3

    valid_id_list_1 = (['idr3'])
    valid_id_list_1.sort()
    valid_id_string_1 = "-".join(valid_id_list_1)

    valid_id_list_2 = ['idr1', 'idr1_1', 'idr2', 'idr5']
    valid_id_list_2.sort()
    valid_id_string_2 = "-".join(valid_id_list_2)

    valid_id_list_3 = ['idr4']
    valid_id_list_3.sort()
    valid_id_string_3 = "-".join(valid_id_list_3)

    for met in mets:
        id_list = met.idList
        id_list.sort()
        id_list_string = "-".join(id_list)
        assert id_list_string in [valid_id_string_1, valid_id_string_2, valid_id_string_3]

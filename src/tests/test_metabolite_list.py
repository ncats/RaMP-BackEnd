import sys

sys.path.append('../src')
import pytest
from unittest.mock import patch
from src.rampEntity.Metabolite import Metabolite
from src.rampEntity.MetaboliteList import MetaboliteList

no_merges = [
    {
        "name": "r1",
        "keys": ["1", "2"]
    },
    {
        "name": "r2",
        "keys": ["3", "4"]
    },
    {
        "name": "r3",
        "keys": ["5", "6"]
    },
    {
        "name": "r4",
        "keys": ["7", "8"]
    }
]

partial_merges = no_merges + [
    {
        "name": "r5",
        "keys": ["1", "4"]
    }
]

all_merges = partial_merges + [
    {
        "name": "r6",
        "keys": ["5", "8"]
    },
    {
        "name": "r7",
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
        metabolite_list.addMetaboliteByAltId(f"source:{m.rampId}", m)


def mock_get_inchi_keys(self, metabolite):
    test_entry = next((item for item in all_merges if item['name'] == metabolite.rampId), None)
    return set(test_entry["keys"])


@patch.object(MetaboliteList, 'get_inchi_keys', mock_get_inchi_keys)
def test_multi_level_mergers(metabolite_list):
    initialize_metabolite_list(metabolite_list, all_merges)
    metabolite_list.collapseMetsOnInchiKeyPrefix()
    mets = metabolite_list.getAllMetabolites()
    assert len(mets) == 7
    assert mets[0].rampId == mets[1].rampId
    assert mets[1].rampId == mets[2].rampId
    assert mets[2].rampId == mets[3].rampId
    assert mets[3].rampId == mets[4].rampId
    assert mets[4].rampId == mets[5].rampId
    assert mets[5].rampId == mets[6].rampId

@patch.object(MetaboliteList, 'get_inchi_keys', mock_get_inchi_keys)
def test_no_mergers(metabolite_list):
    initialize_metabolite_list(metabolite_list, no_merges)
    metabolite_list.collapseMetsOnInchiKeyPrefix()
    mets = metabolite_list.getAllMetabolites()
    assert len(mets) == 4
    assert mets[0].rampId != mets[1].rampId
    assert mets[1].rampId != mets[2].rampId
    assert mets[2].rampId != mets[3].rampId
    assert mets[3].rampId != mets[0].rampId

@patch.object(MetaboliteList, 'get_inchi_keys', mock_get_inchi_keys)
def test_mixed_merges(metabolite_list):
    initialize_metabolite_list(metabolite_list, partial_merges)
    metabolite_list.collapseMetsOnInchiKeyPrefix()
    mets = metabolite_list.getAllMetabolites()
    assert len(mets) == 5
    assert mets[0].rampId == mets[1].rampId
    assert mets[0].rampId == mets[4].rampId
    assert mets[1].rampId != mets[2].rampId
    assert mets[2].rampId != mets[3].rampId
    assert mets[3].rampId != mets[0].rampId

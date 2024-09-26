import sys

from src.rampEntity.Gene import Gene
from src.rampEntity.Molecule import Molecule
from src.rampEntity.Ontology import Ontology

sys.path.append('../src')
import pytest
from unittest.mock import patch
from src.rampEntity.Metabolite import Metabolite
from src.rampEntity.MetaboliteList import MetaboliteList

no_merges = [
    {
        "name": "r1",
        "ids": ['idr1','idr1_1'],
        "chemprops": ['a'],
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
        "chemprops": ['b'],
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
        if "chemprops" in entry:
            for chemprop in entry['chemprops']:
                mol = Molecule()
                mol.id = chemprop
                mol.source = "test"
                mol.names = ["name"]
                m.addChemProps(mol)
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

@patch.object(MetaboliteList, 'get_inchi_keys', mock_get_inchi_keys)
def test_chemprops_merge(metabolite_list):
    initialize_metabolite_list(metabolite_list, partial_merges)
    metabolite_list.collapseMetsOnInchiKeyPrefix()
    mets = metabolite_list.getUniqueMetabolites()
    merged_mets = [m for m in mets if "idr1" in m.idList]
    assert len(merged_mets) == 1
    merged_met = merged_mets[0]
    assert 'test' in merged_met.chemPropsMolecules
    chemprops = merged_met.chemPropsMolecules['test']
    assert len(chemprops) == 2
    ids = [chem.id for chem in chemprops.values()]
    assert 'a' in ids
    assert 'b' in ids

def test_subsume_pathways():
    met1 = Metabolite()
    met1.addPathway('p1', 'reactome')
    met1.addPathway('p2', 'kegg')

    met2 = Metabolite()
    met2.addPathway('p3', 'wiki')
    met2.addPathway('p1', 'reactome')
    met2.addPathway('p4', 'reactome')

    met1.subsumeMetabolite(met2)
    assert 'reactome' in met1.pathways
    assert 'kegg' in met1.pathways
    assert 'wiki' in met1.pathways

    assert len(met1.pathways['reactome']) == 2
    assert len(met1.pathways['wiki']) == 1
    assert len(met1.pathways['kegg']) == 1

def test_subsume_genes():
    g1 = Gene()
    g1.rampId = 'g1'
    g2 = Gene()
    g2.rampId = 'g2'
    g3 = Gene()
    g3.rampId = 'g3'

    met1 = Metabolite()
    met1.addAssociatedGene(g1)
    met1.addAssociatedGene(g2)

    met2 = Metabolite()
    met2.addAssociatedGene(g2)
    met2.addAssociatedGene(g3)

    met1.subsumeMetabolite(met2)
    assert len(met1.associatedGenes) == 3

def test_subsume_ontologies():
    ont1 = Ontology()
    ont1.ontolChild = 'child1'
    ont1.ontolParent = 'parent1'
    ont1.ontolRampId = "o1"

    ont2 = Ontology()
    ont2.ontolChild = 'child2'
    ont2.ontolParent = 'parent2'
    ont2.ontolRampId = "o2"

    ont3 = Ontology()
    ont3.ontolChild = 'child3'
    ont3.ontolParent = 'parent3'
    ont3.ontolRampId = "o3"

    met1 = Metabolite()
    met1.addOntologyTerm(ont1)
    met1.addOntologyTerm(ont2)

    met2 = Metabolite()
    met2.addOntologyTerm(ont2)
    met2.addOntologyTerm(ont3)

    met1.subsumeMetabolite(met2)
    assert len(met1.ontologyTerms) == 3

def test_subsume_classes():
    met1 = Metabolite()
    met1.addMetClass('hmdb', "1a", "2a", "3a")
    met1.addMetClass('hmdb', "1a", "2a", "3b")
    met1.addMetClass('hmdb', "1a", "2b", "3a")
    met1.addMetClass('hmdb', "1a", "2b", "3b")
    met1.addMetClass('hmdb', "1b", "2a", "3a")
    met1.addMetClass('hmdb', "1b", "2a", "3b")
    met1.addMetClass('hmdb', "1b", "2b", "3a")
    met1.addMetClass('hmdb', "1b", "2b", "3b")
    met1.addMetClass('lipidmaps','0', '0', '0')

    met2 = Metabolite()
    met2.addMetClass('hmdb', "1a", "2a", "3a")
    met2.addMetClass('hmdb', "1a", "2a", "3c")
    met2.addMetClass('hmdb', "1a", "2c", "3a")
    met2.addMetClass('hmdb', "1a", "2c", "3c")
    met2.addMetClass('hmdb', "1c", "2a", "3a")
    met2.addMetClass('hmdb', "1c", "2a", "3c")
    met2.addMetClass('hmdb', "1c", "2c", "3a")
    met2.addMetClass('hmdb', "1c", "2c", "3c")
    met2.addMetClass('wiki','0', '0', '0')

    met1.subsumeMetabolite(met2)

    assert 'hmdb' in met1.metClasses
    assert 'lipidmaps' in met1.metClasses
    assert 'wiki' in met1.metClasses

    hmdbClasses = met1.metClasses['hmdb']
    assert ''.join(hmdbClasses.keys()) == '1a1b1c'

    assert ''.join(hmdbClasses['1a']) == '2a2b2c'
    assert ''.join(hmdbClasses['1a']['2a']) == '3a3b3c'
    assert ''.join(hmdbClasses['1a']['2b']) == '3a3b'
    assert ''.join(hmdbClasses['1a']['2c']) == '3a3c'

    assert ''.join(hmdbClasses['1b']) == '2a2b'
    assert ''.join(hmdbClasses['1b']['2a']) == '3a3b'
    assert ''.join(hmdbClasses['1b']['2b']) == '3a3b'

    assert ''.join(hmdbClasses['1c']) == '2a2c'
    assert ''.join(hmdbClasses['1c']['2a']) == '3a3c'
    assert ''.join(hmdbClasses['1c']['2c']) == '3a3c'


def test_subsume_hmdb_status():
    met1 = Metabolite()
    met2 = Metabolite()
    met2.setPriorityHMDBStatus('predicted')

    met1.subsumeMetabolite(met2)
    assert met1.hmdbStatus == 'predicted'

    met3 = Metabolite()
    met1.subsumeMetabolite(met3) # can't overwrite with None
    assert met1.hmdbStatus == 'predicted'

    met4 = Metabolite()
    met4.setPriorityHMDBStatus('expected')
    met1.subsumeMetabolite(met4)
    assert met1.hmdbStatus == 'expected'

    met5 = Metabolite()
    met5.setPriorityHMDBStatus('detected')
    met1.subsumeMetabolite(met5)
    assert met1.hmdbStatus == 'detected'

    met6 = Metabolite()
    met6.setPriorityHMDBStatus('quantified')
    met1.subsumeMetabolite(met6)
    assert met1.hmdbStatus == 'quantified'

    met7 = Metabolite()
    met7.setPriorityHMDBStatus('detected')
    met1.subsumeMetabolite(met1)  # no stepping down
    assert met1.hmdbStatus == 'quantified'
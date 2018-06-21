from spada.methods import get_switches
from spada.network import gene_network

import os
import pickle
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

g = get_switches.GetSwitches(dataPath + 'annotation.pkl')
g.run(dataPath + 'switches')

gn = g._genes
txs = g._txs

def test_readSwitches():

	# number of switches is correct
	assert len(gn.nodes()["ENSG00.5"]["switches"]) == 2
	assert len(gn.nodes()["ENSG03.3"]["switches"]) == 2
	assert len(gn.nodes()["ENSG05.2"]["switches"]) == 1
	assert len(gn.nodes()["ENSG09.6"]["switches"]) == 3

	# check the switches are the right ones
	assert [ x for x in gn.nodes()["ENSG00.5"]["switches"] if x.ctrl == "ENST01.2" and x.case == "ENST02.2" and x.samples == {"A","B","C"} ]
	assert [ x for x in gn.nodes()["ENSG00.5"]["switches"] if x.ctrl == "ENST02.2" and x.case == "ENST01.2" and x.samples == {"D"} ]
	assert [ x for x in gn.nodes()["ENSG03.3"]["switches"] if x.ctrl == "ENST06.1" and x.case == "ENST05.1" and x.samples == {"A","B"} ]
	assert [ x for x in gn.nodes()["ENSG03.3"]["switches"] if x.ctrl == "ENST05.1" and x.case == "ENST06.1" and x.samples == {"C","D"} ]
	assert [ x for x in gn.nodes()["ENSG05.2"]["switches"] if x.ctrl == "ENST08.1" and x.case == "ENST09.1" and x.samples == {"A","B"} ]
	assert [ x for x in gn.nodes()["ENSG09.6"]["switches"] if x.ctrl == "ENST15.1" and x.case == "ENST13.5" and x.samples == {"A","B","C","D"} ]
	assert [ x for x in gn.nodes()["ENSG09.6"]["switches"] if x.ctrl == "ENST16.2" and x.case == "ENST13.5" and x.samples == {"E","F","G"} ]
	assert [ x for x in gn.nodes()["ENSG09.6"]["switches"] if x.ctrl == "ENST16.2" and x.case == "ENST14.5" and x.samples == {"H","I"} ]

def test_add_node():

	old_size = len(gn.nodes())

	assert 'kk' not in gn.nodes()
	assert gn.add_node(full_name = 'kk')

	new_size = len(gn.nodes())

	assert 'kk' in gn.nodes()
	assert old_size == (new_size - 1)
	assert gn.add_node('kk')

def test_add_edge():

	assert not gn.add_edge(full_name1 = 'MADEUP1', full_name2 = 'MADEUP2')
	assert not gn.add_edge(full_name1 = 'kk', full_name2 = 'MADEUP2')
	assert not gn.add_edge(full_name1 = 'MADEUP1', full_name2 = 'kk')
	assert gn.add_edge(full_name1 = 'ENSG03.3', full_name2 = 'kk')

def test_NotImplementedError():

	bad_network = gene_network.GeneNetwork('bad')

	with pytest.raises(NotImplementedError):
		bad_network.nameFilter()

def test_sampleSwitches():

	with pytest.raises(ValueError):
		[ x for x in gn.sampleSwitches(txs) ]

	assert len([ x for x in gn.sampleSwitches(txs, 3) ]) == 3
	assert len([ x for x in gn.sampleSwitches(txs, 4) ]) == 4

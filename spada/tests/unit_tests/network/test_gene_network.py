from spada.network import gene_network

import os
import pickle
import pytest

def test_readSwitches():

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../../data"

	gn = pickle.load(open(dataPath + "/genes.pkl", "rb"))
	gn.createLogger()
	txs = pickle.load(open(dataPath + "/transcripts.pkl", "rb"))
	txs.createLogger()
	gn.flushSwitches()
	gn.readSwitches(dataPath + "/switches", txs)

	# number of switches is correct
	assert len(gn.nodes()["ENSG00.5"]["switches"]) == 2
	assert len(gn.nodes()["ENSG03.3"]["switches"]) == 2
	assert len(gn.nodes()["ENSG05.2"]["switches"]) == 1
	assert len(gn.nodes()["ENSG09.6"]["switches"]) == 3

	# check the switches are the right ones
	assert [ x for x in gn.nodes()["ENSG00.5"]["switches"] if x.nTx == "ENST01.2" and x.tTx == "ENST02.2" and x.samples == {"A","B","C"} ]
	assert [ x for x in gn.nodes()["ENSG00.5"]["switches"] if x.nTx == "ENST02.2" and x.tTx == "ENST01.2" and x.samples == {"D"} ]
	assert [ x for x in gn.nodes()["ENSG03.3"]["switches"] if x.nTx == "ENST06.1" and x.tTx == "ENST05.1" and x.samples == {"A","B"} ]
	assert [ x for x in gn.nodes()["ENSG03.3"]["switches"] if x.nTx == "ENST05.1" and x.tTx == "ENST06.1" and x.samples == {"C","D"} ]
	assert [ x for x in gn.nodes()["ENSG05.2"]["switches"] if x.nTx == "ENST08.1" and x.tTx == "ENST09.1" and x.samples == {"A","B"} ]
	assert [ x for x in gn.nodes()["ENSG09.6"]["switches"] if x.nTx == "ENST15.1" and x.tTx == "ENST13.5" and x.samples == {"A","B","C","D"} ]
	assert [ x for x in gn.nodes()["ENSG09.6"]["switches"] if x.nTx == "ENST16.2" and x.tTx == "ENST13.5" and x.samples == {"E","F","G"} ]
	assert [ x for x in gn.nodes()["ENSG09.6"]["switches"] if x.nTx == "ENST16.2" and x.tTx == "ENST14.5" and x.samples == {"H","I"} ]

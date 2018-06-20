from spada.io import io
from spada.methods import create_network, get_switches

import pickle
import pytest
import os

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_printSwitches():

	c = create_network.CreateNetwork("test", "gencode", new = False)
	g = get_switches.GetSwitches(c._genes, c._txs)
	g.run(dataPath + 'switches')

	io.printSwitches(g._genes, g._txs)

	switches = [ x for x in io.readTable("switches_spada.tsv") ]

	assert len(switches) == 8

	os.remove("switches_spada.tsv")

def test_getGene2Tx():

	genes = pickle.load(open(dataPath + 'genes.pkl', 'rb'))
	txs = pickle.load(open(dataPath + 'transcripts.pkl', 'rb'))
	gene2tx = io.getGene2Tx(txs)

	assert len(genes.nodes()) == len(gene2tx)
	assert len(txs.nodes()) == len([ t for g,T in gene2tx.items() for t in T ])

def test_readSamples():

	EXPR = open(dataPath + 'expression', "r")
	ids = io.readSamples(EXPR)

	assert ids == ['1','2','3']

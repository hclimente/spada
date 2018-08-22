from spada.io import io
from spada.methods import get_switches, method

import pickle
import pytest
import os

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_printSwitches():

	g = get_switches.GetSwitches(dataPath + 'annotation.pklz')
	g.run(dataPath + 'switches')

	io.printSwitches(g._genes, g._txs)

	switches = [ x for x in io.readTable("switches_spada.tsv") ]

	assert len(switches) == 8

	os.remove("switches_spada.tsv")

def test_getGene2Tx():

	m = method.Method('test_getGene2Tx', dataPath + 'annotation.pklz')
	gene2tx = io.getGene2Tx(m._txs)

	assert len(m._genes.nodes()) == len(gene2tx)
	assert len([ t for t,i in m._txs.transcripts() if i['canonical'] ]) == len([ t for g,T in gene2tx.items() for t in T ])

def test_readSamples():

	EXPR = open(dataPath + 'expression', "r")
	ids = io.readSamples(EXPR)

	assert ids == ['1','2','3']

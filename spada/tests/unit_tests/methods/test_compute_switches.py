from spada.methods import compute_switches

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

c = compute_switches.ComputeSwitches(dataPath + "genes.pkl",
									 dataPath + "transcripts.pkl",
									 dataPath + 'expression',
									 dataPath + 'expression')

def test_run():

	# switches are flushed
	assert not [ x for x,i in c._genes.nodes(data=True) if i['switches'] ]

	c.run(0.1)

	assert os.stat('switches_spada.tsv').st_size > 0
	os.remove('switches_spada.tsv')
	assert [ x for x,i in c._genes.nodes(data=True) if i['switches'] ]

def test_getGene2Tx():

	gene2tx = c.getGene2Tx()

	assert len(c._genes.nodes()) == len(gene2tx)
	assert len(c._txs.nodes()) == len([ t for g,T in gene2tx.items() for t in T ])

def test_readSamples():

	EXPR = open(dataPath + 'expression', "r")
	ids = c.readSamples(EXPR)

	assert ids == ['1','2','3']

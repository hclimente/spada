#!/usr/bin/env python

from spada.methods import get_switches

import os
import pytest

def test_init():

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../data/"


	g = get_switches.GetSwitches(dataPath + "genes.pkl",
								 dataPath + "transcripts.pkl")
	g.run(dataPath + "switches")

	# number of switches is correct
	assert len(g._genes.nodes()["ENSG00000223972.5"]["switches"]) == 2
	assert len(g._genes.nodes()["ENSG00000243485.3"]["switches"]) == 2
	assert len(g._genes.nodes()["ENSG00000237613.2"]["switches"]) == 1
	assert len(g._genes.nodes()["ENSG00000238009.6"]["switches"]) == 3

	# check the switches are the right ones
	assert [ x for x in g._genes.nodes()["ENSG00000223972.5"]["switches"] if x.nTx == "ENST00000456328.2" and x.tTx == "ENST00000450305.2" and x.samples == ["A","B","C"] ]
	assert [ x for x in g._genes.nodes()["ENSG00000223972.5"]["switches"] if x.nTx == "ENST00000450305.2" and x.tTx == "ENST00000456328.2" and x.samples == ["D"] ]
	assert [ x for x in g._genes.nodes()["ENSG00000243485.3"]["switches"] if x.nTx == "ENST00000469289.1" and x.tTx == "ENST00000473358.1" and x.samples == ["A","B"] ]
	assert [ x for x in g._genes.nodes()["ENSG00000243485.3"]["switches"] if x.nTx == "ENST00000473358.1" and x.tTx == "ENST00000469289.1" and x.samples == ["C","D"] ]
	assert [ x for x in g._genes.nodes()["ENSG00000237613.2"]["switches"] if x.nTx == "ENST00000417324.1" and x.tTx == "ENST00000461467.1" and x.samples == ["A","B"] ]
	assert [ x for x in g._genes.nodes()["ENSG00000238009.6"]["switches"] if x.nTx == "ENST00000610542.1" and x.tTx == "ENST00000466430.5" and x.samples == ["A","B","C","D"] ]
	assert [ x for x in g._genes.nodes()["ENSG00000238009.6"]["switches"] if x.nTx == "ENST00000453576.2" and x.tTx == "ENST00000466430.5" and x.samples == ["E","F","G"] ]
	assert [ x for x in g._genes.nodes()["ENSG00000238009.6"]["switches"] if x.nTx == "ENST00000453576.2" and x.tTx == "ENST00000477740.5" and x.samples == ["H","I"] ]

	# switches have at most 1 candidate, some have noise
	assert len([ x for x in g._genes.nodes()["ENSG00000223972.5"]["switches"] if x.isCandidate ]) == 1
	assert len([ x for x in g._genes.nodes()["ENSG00000223972.5"]["switches"] if x.isNoise ]) == 1
	assert len([ x for x in g._genes.nodes()["ENSG00000243485.3"]["switches"] if x.isCandidate ]) == 0
	assert len([ x for x in g._genes.nodes()["ENSG00000243485.3"]["switches"] if x.isNoise ]) == 2
	assert len([ x for x in g._genes.nodes()["ENSG00000237613.2"]["switches"] if x.isCandidate ]) == 1
	assert len([ x for x in g._genes.nodes()["ENSG00000238009.6"]["switches"] if x.isCandidate ]) == 1
	assert len([ x for x in g._genes.nodes()["ENSG00000238009.6"]["switches"] if x.isNoise ]) == 1

	# candidate is the expected switch
	assert g._genes.getSwitch("ENSG00000223972.5","ENST00000456328.2","ENST00000450305.2").isCandidate
	assert g._genes.getSwitch("ENSG00000237613.2","ENST00000417324.1","ENST00000461467.1").isCandidate
	assert g._genes.getSwitch("ENSG00000238009.6","ENST00000453576.2","ENST00000466430.5").isCandidate

	# noise are the expected switch
	assert g._genes.getSwitch("ENSG00000223972.5","ENST00000450305.2","ENST00000456328.2").isNoise
	assert g._genes.getSwitch("ENSG00000237613.2","ENST00000417324.1","ENST00000461467.1").isNoise
	assert g._genes.getSwitch("ENSG00000238009.6","ENST00000453576.2","ENST00000477740.5").isNoise

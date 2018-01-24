#!/usr/bin/env python

from spada.methods import get_switches

import os
import pytest

def test_init():

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../data/"


	g = get_switches.GetSwitches(dataPath + "genes.pkl",
								 dataPath + "transcripts.pkl")
	g._genes.clean()
	g.run(dataPath + "switches")

	# number of switches is correct
	assert len(g._genes.nodes()["ENSG00.5"]["switches"]) == 2
	assert len(g._genes.nodes()["ENSG03.3"]["switches"]) == 2
	assert len(g._genes.nodes()["ENSG05.2"]["switches"]) == 1
	assert len(g._genes.nodes()["ENSG09.6"]["switches"]) == 3

	# check the switches are the right ones
	assert [ x for x in g._genes.nodes()["ENSG00.5"]["switches"] if x.nTx == "ENST01.2" and x.tTx == "ENST02.2" and x.samples == ["A","B","C"] ]
	assert [ x for x in g._genes.nodes()["ENSG00.5"]["switches"] if x.nTx == "ENST02.2" and x.tTx == "ENST01.2" and x.samples == ["D"] ]
	assert [ x for x in g._genes.nodes()["ENSG03.3"]["switches"] if x.nTx == "ENST06.1" and x.tTx == "ENST05.1" and x.samples == ["A","B"] ]
	assert [ x for x in g._genes.nodes()["ENSG03.3"]["switches"] if x.nTx == "ENST05.1" and x.tTx == "ENST06.1" and x.samples == ["C","D"] ]
	assert [ x for x in g._genes.nodes()["ENSG05.2"]["switches"] if x.nTx == "ENST08.1" and x.tTx == "ENST09.1" and x.samples == ["A","B"] ]
	assert [ x for x in g._genes.nodes()["ENSG09.6"]["switches"] if x.nTx == "ENST15.1" and x.tTx == "ENST13.5" and x.samples == ["A","B","C","D"] ]
	assert [ x for x in g._genes.nodes()["ENSG09.6"]["switches"] if x.nTx == "ENST16.2" and x.tTx == "ENST13.5" and x.samples == ["E","F","G"] ]
	assert [ x for x in g._genes.nodes()["ENSG09.6"]["switches"] if x.nTx == "ENST16.2" and x.tTx == "ENST14.5" and x.samples == ["H","I"] ]

	# switches have at most 1 candidate, some have noise
	assert len([ x for x in g._genes.nodes()["ENSG00.5"]["switches"] if x.isCandidate ]) == 1
	assert len([ x for x in g._genes.nodes()["ENSG00.5"]["switches"] if x.isNoise ]) == 1
	assert len([ x for x in g._genes.nodes()["ENSG03.3"]["switches"] if x.isCandidate ]) == 0
	assert len([ x for x in g._genes.nodes()["ENSG03.3"]["switches"] if x.isNoise ]) == 2
	assert len([ x for x in g._genes.nodes()["ENSG05.2"]["switches"] if x.isCandidate ]) == 1
	assert len([ x for x in g._genes.nodes()["ENSG09.6"]["switches"] if x.isCandidate ]) == 1
	assert len([ x for x in g._genes.nodes()["ENSG09.6"]["switches"] if x.isNoise ]) == 1

	# candidate is the expected switch
	assert g._genes.getSwitch("ENSG00.5","ENST01.2","ENST02.2").isCandidate
	assert g._genes.getSwitch("ENSG05.2","ENST08.1","ENST09.1").isCandidate
	assert g._genes.getSwitch("ENSG09.6","ENST16.2","ENST13.5").isCandidate

	# noise are the expected switch
	assert g._genes.getSwitch("ENSG00.5","ENST02.2","ENST01.2").isNoise
	assert g._genes.getSwitch("ENSG05.2","ENST08.1","ENST09.1").isNoise
	assert g._genes.getSwitch("ENSG09.6","ENST16.2","ENST14.5").isNoise

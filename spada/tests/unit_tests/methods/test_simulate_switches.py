#!/usr/bin/env python

from spada.methods import simulate_switches

import itertools
import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

sim = simulate_switches.SimulateSwitches(dataPath + "genes.pkl",
										 dataPath + "transcripts.pkl")

def test_run():

	sim.run(dataPath + "expression", dataPath + "expression_case", 'random', 0.1)
	os.remove("switches_simulated_random.tsv")
	sim._genes.flushSwitches()

	sim.run(dataPath + "expression", dataPath + "expression_case", 'fix_expressed', 0.1)
	os.remove("switches_simulated_fix_expressed.tsv")
	sim._genes.flushSwitches()

	sim.run(dataPath + "expression", dataPath + "expression_case", 'fix_main', 0.1)
	os.remove("switches_simulated_fix_main.tsv")

def test_sampleTranscripts_fixControl():

	cases = ['B', 'C', 'D']
	switches = sim.sampleTranscripts_fixControl('A', cases)

	assert len(switches) == 3
	for ctrl,case in switches:
		assert ctrl == 'A'
		assert case in cases

def test_sampleTranscripts_random():

	txs = ['A', 'B', 'C', 'D']
	switches = sim.sampleTranscripts_random(txs)

	allswitches = set([ x for x in itertools.combinations(txs, 2) ])

	assert len(switches) == 5
	assert len(set(switches) & allswitches) == 5

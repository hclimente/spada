#!/usr/bin/env python

from spada.methods import simulate_switches

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_sampleTranscripts_fixNormal():

	sim = simulate_switches.SimulateSwitches(dataPath + "genes.pkl",
											 dataPath + "transcripts.pkl")
	switches = sim.sampleTranscripts_fixNormal()

def test_sampleTranscripts_random():

	sim = simulate_switches.SimulateSwitches(dataPath + "genes.pkl",
											 dataPath + "transcripts.pkl")
	switches = sim.sampleTranscripts_random()

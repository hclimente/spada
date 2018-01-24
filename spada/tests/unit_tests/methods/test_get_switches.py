#!/usr/bin/env python

from spada.methods import get_switches

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_run():

	g = get_switches.GetSwitches(dataPath + "genes.pkl",
								 dataPath + "transcripts.pkl")
	g.run(dataPath + "switches")

	assert os.stat("genes.pkl").st_size > 0
	#os.remove("genes.pkl")

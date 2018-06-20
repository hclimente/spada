from spada.methods import compute_switches

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

c = compute_switches.ComputeSwitches(dataPath + "genes.pkl",
									 dataPath + "transcripts.pkl")

def test_run():

	# switches are flushed
	assert not [ x for x,i in c._genes.nodes(data=True) if i['switches'] ]

	c.run(dataPath + 'expression', dataPath + 'expression_case', 0.1)

	assert os.stat('switches_spada.tsv').st_size > 0
	os.remove('switches_spada.tsv')
	assert [ x for x,i in c._genes.nodes(data=True) if i['switches'] ]

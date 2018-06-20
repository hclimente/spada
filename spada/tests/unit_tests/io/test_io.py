from spada.io import io
from spada.methods import create_network, get_switches

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

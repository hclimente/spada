from spada.io import io
from spada.methods import create_network

import pytest
import os

def test_printSwitches():

	c = create_network.CreateNetwork("test", "gencode", new = False)

	io.printSwitches(c._genes, c._txs)

	switches = [ x for x in io.readTable("switches_spada.tsv") ]

	assert len(switches) == 8

	os.remove("switches_spada.tsv")

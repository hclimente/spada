#!/usr/bin/env python

from spada.methods import get_switches

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_run():

	g = get_switches.GetSwitches(dataPath + 'annotation.pklz')

	g.run(dataPath + "switches")

	assert [ x for x in g._genes.switches(g._txs) ]

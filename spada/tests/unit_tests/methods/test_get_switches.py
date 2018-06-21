#!/usr/bin/env python

from spada.methods import get_switches

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_run():

	g = get_switches.GetSwitches(dataPath + 'annotation.pkl')

	g.run(dataPath + "switches")

	assert os.stat('annotation.pkl').st_size > 0
	os.remove('annotation.pkl')

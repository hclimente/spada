from spada.methods import summary

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_init():

	s = summary.Summary(dataPath + "genes.pkl",
						dataPath + "transcripts.pkl")

def test_run():

	s = summary.Summary(dataPath + "genes.pkl",
						dataPath + "transcripts.pkl")
	s.run()

	os.remove('isoform_length.tsv')
	os.remove('exons.tsv')
	os.remove('exons_new.tsv')
	os.remove('structural_features.tsv')
	os.remove("switches_spada.tsv")

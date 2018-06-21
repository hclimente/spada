from spada.io import io
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
	s.run(dataPath + "expression", dataPath + "expression_case")

	os.remove('proteome_features.tsv')
	os.remove('isoform_length.tsv')
	os.remove('exons.tsv')
	os.remove('exons_new.tsv')
	os.remove('structural_features.tsv')
	os.remove("switches_spada.tsv")

def test_proteomeStatistics():

	s = summary.Summary(dataPath + "genes.pkl",
						dataPath + "transcripts.pkl")
	s.proteomeStatistics(dataPath + "expression", dataPath + "expression_case")
	proteome = [ x for x in io.readTable("proteome_features.tsv") ]

	assert len(proteome) == 4 # ENSG00.5 isnot complete, as we don't have info for JKLM-.3 and DE_FG_HI.2

	assert len([ x for x in proteome if x['Feature_type'] == 'Prosite' ]) == 3
	assert len([ x for x in proteome if x['Feature_type'] == 'Pfam' ]) == 1

	# ENST08.1 goes first in the file, and both have same median
	assert len([ x for x in proteome if x['Transcript'] == 'ENST08.1' ]) == 2
	assert len([ x for x in proteome if x['Transcript'] == 'ENST09.1' ]) == 0

	assert len([ x for x in proteome if x['Transcript'] == 'ENST18.3' ]) == 1

	os.remove('proteome_features.tsv')

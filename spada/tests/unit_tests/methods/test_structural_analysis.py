#!/usr/bin/env python

from spada.methods import structural_analysis
from spada.biological_entities import switch

import os
import shutil
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_init():

	s = structural_analysis.StructuralAnalysis(dataPath + "genes.pkl",
								 			   dataPath + "transcripts.pkl")

	assert os.stat("structural_analysis/pfam_analysis.tsv").st_size == 0
	assert os.stat("structural_analysis/prosite_analysis.tsv").st_size == 0
	assert os.stat("structural_analysis/idr_analysis.tsv").st_size == 0
	shutil.rmtree("structural_analysis")

	s = structural_analysis.StructuralAnalysis(dataPath + "genes.pkl",
									 		   dataPath + "transcripts.pkl",
											   isRandom = True)

	assert os.stat("structural_analysis/pfam_analysis_random.tsv").st_size == 0
	assert os.stat("structural_analysis/prosite_analysis_random.tsv").st_size == 0
	assert os.stat("structural_analysis/idr_analysis_random.tsv").st_size == 0
	shutil.rmtree("structural_analysis")

def test_run():

	s = structural_analysis.StructuralAnalysis(dataPath + "genes.pkl",
								 			   dataPath + "transcripts.pkl")
	s.run()

	assert os.stat("structural_analysis/pfam_analysis.tsv").st_size > 0
	assert os.stat("structural_analysis/prosite_analysis.tsv").st_size > 0
	assert os.stat("structural_analysis/idr_analysis.tsv").st_size > 0
	shutil.rmtree("structural_analysis")

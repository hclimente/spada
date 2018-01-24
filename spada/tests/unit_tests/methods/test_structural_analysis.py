#!/usr/bin/env python

from spada.methods import structural_analysis
from spada import utils

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_run():

	s = structural_analysis.StructuralAnalysis(dataPath + "genes.pkl",
								 			   dataPath + "transcripts.pkl")
	s.run()

	# pfams
	pfam = [ x for x in utils.readTable("pfam_analysis.tsv") ]
	assert len(pfam) == 8
	assert [ x for x in pfam if x[1] == "ENST01.2" and x[2] == "ENST02.2" and x[3] == "Nothing" and x[4] == "D1" ]
	assert [ x for x in pfam if x[1] == "ENST01.2" and x[2] == "ENST02.2" and x[3] == "Lost_in_tumor" and x[4] == "D1" ]
	assert [ x for x in pfam if x[1] == "ENST01.2" and x[2] == "ENST02.2" and x[3] == "Gained_in_tumor" and x[4] == "D4" ]
	assert [ x for x in pfam if x[1] == "ENST01.2" and x[2] == "ENST02.2" and x[3] == "Gained_in_tumor" and x[4] == "D2" ]

	assert [ x for x in pfam if x[1] == "ENST02.2" and x[2] == "ENST01.2" and x[3] == "Nothing" and x[4] == "D1" ]
	assert [ x for x in pfam if x[1] == "ENST02.2" and x[2] == "ENST01.2" and x[3] == "Gained_in_tumor" and x[4] == "D1" ]
	assert [ x for x in pfam if x[1] == "ENST02.2" and x[2] == "ENST01.2" and x[3] == "Lost_in_tumor" and x[4] == "D4" ]
	assert [ x for x in pfam if x[1] == "ENST02.2" and x[2] == "ENST01.2" and x[3] == "Lost_in_tumor" and x[4] == "D2" ]

	# prosite
	prosite = [ x for x in utils.readTable("prosite_analysis.tsv") ]
	assert len(prosite) == 3
	assert [ x for x in prosite if x[1] == "ENST08.1" and x[2] == "ENST09.1" and x[3] == "Nothing" and x[4] == "P1" ]
	assert [ x for x in prosite if x[1] == "ENST08.1" and x[2] == "ENST09.1" and x[3] == "Lost_in_tumor" and x[4] == "P1" ]
	assert [ x for x in prosite if x[1] == "ENST08.1" and x[2] == "ENST09.1" and x[3] == "Gained_in_tumor" and x[4] == "P2" ]

	# idr
	idr = [ x for x in utils.readTable("idr_analysis.tsv") ]
	assert len(idr) == 6
	assert [ x for x in idr if x[1] == "ENST16.2" and x[2] == "ENST14.5" and x[3] == "Lost_in_tumor" and x[4] == "ABcd" ]
	assert len([ x for x in idr if x[1] == "ENST16.2" and x[2] == "ENST14.5"]) == 1

	os.remove("pfam_analysis.tsv")
	os.remove("prosite_analysis.tsv")
	os.remove("idr_analysis.tsv")

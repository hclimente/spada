from spada.methods import structural_analysis
from spada.io import io

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_featureAnalysis():

	s = structural_analysis.StructuralAnalysis(dataPath + "genes.pkl",
								 			   dataPath + "transcripts.pkl")
	s.featureAnalysis()

	# pfams
	pfam = [ x for x in io.readTable("pfam_analysis.tsv") ]
	assert len(pfam) == 8
	assert [ x for x in pfam if x['Normal_transcript'] == "ENST01.2" and x['Tumor_transcript'] == "ENST02.2" and x['What'] == "Nothing" and x['Feature'] == "D1" ]
	assert [ x for x in pfam if x['Normal_transcript'] == "ENST01.2" and x['Tumor_transcript'] == "ENST02.2" and x['What'] == "Lost_in_tumor" and x['Feature'] == "D1" ]
	assert [ x for x in pfam if x['Normal_transcript'] == "ENST01.2" and x['Tumor_transcript'] == "ENST02.2" and x['What'] == "Gained_in_tumor" and x['Feature'] == "D4" ]
	assert [ x for x in pfam if x['Normal_transcript'] == "ENST01.2" and x['Tumor_transcript'] == "ENST02.2" and x['What'] == "Gained_in_tumor" and x['Feature'] == "D2" ]

	assert [ x for x in pfam if x['Normal_transcript'] == "ENST02.2" and x['Tumor_transcript'] == "ENST01.2" and x['What'] == "Nothing" and x['Feature'] == "D1" ]
	assert [ x for x in pfam if x['Normal_transcript'] == "ENST02.2" and x['Tumor_transcript'] == "ENST01.2" and x['What'] == "Gained_in_tumor" and x['Feature'] == "D1" ]
	assert [ x for x in pfam if x['Normal_transcript'] == "ENST02.2" and x['Tumor_transcript'] == "ENST01.2" and x['What'] == "Lost_in_tumor" and x['Feature'] == "D4" ]
	assert [ x for x in pfam if x['Normal_transcript'] == "ENST02.2" and x['Tumor_transcript'] == "ENST01.2" and x['What'] == "Lost_in_tumor" and x['Feature'] == "D2" ]

	# prosite
	prosite = [ x for x in io.readTable("prosite_analysis.tsv") ]
	assert len(prosite) == 3
	assert [ x for x in prosite if x['Normal_transcript'] == "ENST08.1" and x['Tumor_transcript'] == "ENST09.1" and x['What'] == "Nothing" and x['Feature'] == "P1" ]
	assert [ x for x in prosite if x['Normal_transcript'] == "ENST08.1" and x['Tumor_transcript'] == "ENST09.1" and x['What'] == "Lost_in_tumor" and x['Feature'] == "P1" ]
	assert [ x for x in prosite if x['Normal_transcript'] == "ENST08.1" and x['Tumor_transcript'] == "ENST09.1" and x['What'] == "Gained_in_tumor" and x['Feature'] == "P2" ]

	# idr
	idr = [ x for x in io.readTable("idr_analysis.tsv") ]
	assert len(idr) == 6
	assert [ x for x in idr if x['Normal_transcript'] == "ENST16.2" and x['Tumor_transcript'] == "ENST14.5" and x['What'] == "Lost_in_tumor" and x['Sequence'] == "ABcd" ]
	assert len([ x for x in idr if x['Normal_transcript'] == "ENST16.2" and x['Tumor_transcript'] == "ENST14.5"]) == 1

	os.remove("pfam_analysis.tsv")
	os.remove("prosite_analysis.tsv")
	os.remove("idr_analysis.tsv")

def test_ppiAnalysis():

	s = structural_analysis.StructuralAnalysis(dataPath + "genes.pkl",
								 			   dataPath + "transcripts.pkl")
	s.ppiAnalysis()

	ddi = [ x for x in io.readTable("ppi_analysis.tsv") ]
	assert len(ddi) == 6
	assert len([ x for x in ddi if x['What'] == "Unaffected" ]) == 2
	assert len([ x for x in ddi if x['What'] == "Affected" ]) == 2
	assert len([ x for x in ddi if x['What'] == "Lost_in_tumor" ]) == 1
	assert len([ x for x in ddi if x['What'] == "Gained_in_tumor" ]) == 1

	os.remove("ppi_analysis.tsv")

def test_analyzeDDIs():

	s = structural_analysis.StructuralAnalysis(dataPath + "genes.pkl",
								 			   dataPath + "transcripts.pkl")

	thisSwitch = s._genes.nodes()["ENSG00.5"]["switches"][0]
	ddiChanges = s.analyzeDDIs(thisSwitch)

	assert len(ddiChanges) == 3
	assert ddiChanges['ENST20.1']['what'] == 'Gained_in_tumor'
	assert ddiChanges['ENST20.1']['nDDIs'] == set()
	assert ddiChanges['ENST20.1']['tDDIs'] == {'D2@D2','D4@D2'}
	assert ddiChanges['ENST20.1']['bothDDIs'] == set()

	assert ddiChanges['ENST01.2']['what'] == 'Unaffected'
	assert ddiChanges['ENST01.2']['nDDIs'] == set()
	assert ddiChanges['ENST01.2']['tDDIs'] == set()
	assert ddiChanges['ENST01.2']['bothDDIs'] == {'D1@D1'}

	assert ddiChanges['ENST02.2']['what'] == 'Affected'
	assert ddiChanges['ENST02.2']['nDDIs'] == set()
	assert ddiChanges['ENST02.2']['tDDIs'] == {'D2@D2','D2@D4','D4@D2'}
	assert ddiChanges['ENST02.2']['bothDDIs'] == {'D1@D1'}

	thisSwitch = s._genes.nodes()["ENSG00.5"]["switches"][1]
	ddiChanges = s.analyzeDDIs(thisSwitch)

	assert len(ddiChanges) == 3
	assert ddiChanges['ENST20.1']['what'] == 'Lost_in_tumor'
	assert ddiChanges['ENST20.1']['nDDIs'] == {'D2@D2','D4@D2'}
	assert ddiChanges['ENST20.1']['tDDIs'] == set()
	assert ddiChanges['ENST20.1']['bothDDIs'] == set()

	assert ddiChanges['ENST01.2']['what'] == 'Unaffected'
	assert ddiChanges['ENST01.2']['nDDIs'] == set()
	assert ddiChanges['ENST01.2']['tDDIs'] == set()
	assert ddiChanges['ENST01.2']['bothDDIs'] == {'D1@D1'}

	assert ddiChanges['ENST02.2']['what'] == 'Affected'
	assert ddiChanges['ENST02.2']['nDDIs'] == {'D2@D2','D2@D4','D4@D2'}
	assert ddiChanges['ENST02.2']['tDDIs'] == set()
	assert ddiChanges['ENST02.2']['bothDDIs'] == {'D1@D1'}

def test_analyzeDDIs():

	s = structural_analysis.StructuralAnalysis(dataPath + "genes.pkl",
								 			   dataPath + "transcripts.pkl")
	s.proteomeStatistics()
	proteome = [ x for x in io.readTable("proteome_features.tsv") ]

	assert len(proteome) == 6
	assert len([ x for x in proteome if x['Feature_type'] == 'Pfam' ]) == 3
	assert len([ x for x in proteome if x['Transcript'] == 'ENST01.2' ]) == 2
	assert len([ x for x in proteome if x['Transcript'] == 'ENST02.2' ]) == 0

	assert len([ x for x in proteome if x['Feature_type'] == 'Prosite' ]) == 3
	assert len([ x for x in proteome if x['Transcript'] == 'ENST18.3' ]) == 1

	os.remove('proteome_features.tsv')

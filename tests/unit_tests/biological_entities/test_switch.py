#!/usr/bin/env python

from spada.biological_entities import switch

import math
import pytest

a1Info = {	"gene_id":			"A",
			"exons": 			[[1,31],[40,60],[70,101]],
		  	"txCoords": 		[1,100],
		  	"CDS":				[21,82],
		  	"strand":			"+",
		  	"chr":				1,
			"start_codon":		21,
			"stop_codon":		83,
			"proteinSequence":	"ABCDEFGHIJNOPQR",
			"Pfam":				{"D1": [[1,3],[5,7]]},
			"Prosite":			{"P1": [[7,9]], "P2": [[8,9]]},
			"IDR": 				{"EFGHIJNOPQ": [[5, 14]]} }
a2Info = {	"gene_id":			"A",
			"exons": 			[[1,31],[49,101],[120,151]],
		  	"txCoords": 		[1,150],
		  	"CDS":				[21,124],
		  	"strand":			"+",
		  	"chr":				1,
			"start_codon":		21,
			"stop_codon":		125,
			"proteinSequence":	"ABCDGHIJKLMNOPKRSTUVWXY",
			"Pfam":				{"D1": [[1,3]], "D2": [[7,9]]},
			"Prosite":			{"P1": [[4,6],[8,10]]},
			"IDR":				{"GHIJKLMNO": [[5, 13]]} }

def test_init():

	s = switch.IsoformSwitch("A", "B", [1,2,3,4])

	assert s.ctrl == "A"
	assert s.case == "B"
	assert s.samples == [1,2,3,4]
	assert s.ctrlTranscript == None
	assert s.caseTranscript == None
	assert s.ctrlIsoform == None
	assert s.caseIsoform == None

	s = switch.IsoformSwitch("Y", "Z", ["P1","P2","P3"])

	assert s.ctrl == "Y"
	assert s.case == "Z"
	assert s.samples == ["P1","P2","P3"]
	assert s.ctrlTranscript == None
	assert s.caseTranscript == None
	assert s.ctrlIsoform == None
	assert s.caseIsoform == None

def test_addTxInfo():

	s = switch.IsoformSwitch("A", "B", [1,2,3,4])
	s.addTxInfo(a1Info, a2Info)

	# properties of transcripts are kept
	assert s.ctrlTranscript != None
	assert s.ctrlTranscript.name == "A"
	assert len(s.ctrlTranscript.cds) == 3 * len(a1Info["proteinSequence"])
	assert min(s.ctrlTranscript.cds) == a1Info["CDS"][0]
	assert max(s.ctrlTranscript.cds) == a1Info["CDS"][1]

	assert s.caseTranscript != None
	assert s.caseTranscript.name == "B"
	assert len(s.caseTranscript.cds) == 3 * len(a2Info["proteinSequence"])
	assert min(s.caseTranscript.cds) == a2Info["CDS"][0]
	assert max(s.caseTranscript.cds) == a2Info["CDS"][1]

	# properties of proteins are kept
	assert s.ctrlIsoform != None
	assert s.ctrlIsoform.tx == "A"
	assert s.ctrlIsoform.seq == a1Info["proteinSequence"]
	assert "".join( x.res for x in s.ctrlIsoform.structure ) == a1Info["proteinSequence"]
	assert s.ctrlIsoform._pfam == a1Info["Pfam"]
	assert s.ctrlIsoform._prosite == a1Info["Prosite"]
	assert s.ctrlIsoform._idr == a1Info["IDR"]

	assert s.caseIsoform != None
	assert s.caseIsoform.tx == "B"
	assert s.caseIsoform.seq == a2Info["proteinSequence"]
	assert "".join( x.res for x in s.caseIsoform.structure ) == a2Info["proteinSequence"]
	assert s.caseIsoform._pfam == a2Info["Pfam"]
	assert s.caseIsoform._prosite == a2Info["Prosite"]
	assert s.caseIsoform._idr == a2Info["IDR"]

def test_analyzeDomains():

	thisSwitch = switch.IsoformSwitch("A1", "A2", [1,2,3])
	thisSwitch.addTxInfo(a1Info, a2Info)

	pfams = thisSwitch.analyzeDomains("Pfam")
	assert len(pfams) == 3
	assert len([ x for x in pfams if x["feature"] == "D1"]) == 2
	assert len([ x for x in pfams if x["feature"] == "D2"]) == 1
	assert math.isnan([ x['nStart'] for x in pfams if x["feature"] == "D2"][0])
	assert math.isnan([ x['nEnd'] for x in pfams if x["feature"] == "D2"][0])
	assert [ (x['tStart'],x['tEnd']) for x in pfams if x["feature"] == "D2"][0] == (7,9)

	prosites = thisSwitch.analyzeDomains("Prosite")
	assert len(prosites) == 3
	assert len([ x for x in prosites if x["feature"] == "P1"]) == 2
	assert len([ x for x in prosites if x["feature"] == "P2"]) == 1
	assert math.isnan([ x['tStart'] for x in prosites if x["feature"] == "P2"][0])
	assert math.isnan([ x['tEnd'] for x in prosites if x["feature"] == "P2"][0])
	assert [ (x['nStart'],x['nEnd']) for x in prosites if x["feature"] == "P2"][0] == (8,9)

def test_analyzeIDR():

	thisSwitch = switch.IsoformSwitch("A1", "A2", [1,2,3])
	thisSwitch.addTxInfo(a1Info, a2Info)

	idrs = thisSwitch.analyzeIDR(0.2)
	assert len(idrs) == 2
	assert len([ x for x in idrs if x["feature"] == "EFGhijnopq" ]) == 1
	assert [ (x["start"], x["end"]) for x in idrs if x["feature"] == "EFGhijnopq" ] == [(5,14)]
	assert [ (x["m"], x["M"],x["J"]) for x in idrs if x["feature"] == "EFGhijnopq" ] == [(1, 3/10, 3/10)]
	assert len([ x for x in idrs if x["feature"] == "ghijKLMno" ]) == 1
	assert [ (x["start"], x["end"]) for x in idrs if x["feature"] == "ghijKLMno" ] == [(5,13)]
	assert [ (x["m"], x["M"],x["J"]) for x in idrs if x["feature"] == "ghijKLMno" ] == [(1, 3/9, 3/9)]

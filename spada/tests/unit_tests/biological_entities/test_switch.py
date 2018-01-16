#!/usr/bin/env python

from spada.biological_entities import switch

import pytest

a1Info = {	"gene_id":			"A",
			"exonStructure": 	[[1,31],[40,61],[70,101]],
		  	"txCoords": 		[1,100],
		  	"cdsCoords":		[21,87],
		  	"strand":			"+",
		  	"chr":				1,
			"proteinSequence":	"ABCDEFGHIJNOPQR",
			"Pfam":				{"D1": [[1,3],[5,7]]},
			"Prosite":			{"P1": [[7,9]], "P2": [[8,9]]},
			"IDR": 				{"EFGHIJNOPQ": [[5, 14]]} }
a2Info = {	"gene_id":			"A",
			"exonStructure": 	[[1,31],[49,101],[120,151]],
		  	"txCoords": 		[1,150],
		  	"cdsCoords":		[21,130],
		  	"strand":			"+",
		  	"chr":				1,
			"proteinSequence":	"ABCDGHIJKLMNOPKRSTUVWXY",
			"Pfam":				{"D1": [[1,3]], "D2": [[7,9]]},
			"Prosite":			{"P1": [[4,6],[8,10]]},
			"IDR":				{"GHIJKLMNO": [[5, 13]]} }

def test_analyzeDomains():

	thisSwitch = switch.IsoformSwitch("A1", "A2", [1,2,3])
	thisSwitch.addTxInfo(a1Info, a2Info)

	pfams = thisSwitch.analyzeDomains("Pfam")
	assert len(pfams) == 3
	assert len([ x for x in pfams if x["feature"] == "D1"]) == 2
	assert len([ x for x in pfams if x["feature"] == "D2"]) == 1

	prosites = thisSwitch.analyzeDomains("Prosite")
	assert len(prosites) == 3
	assert len([ x for x in prosites if x["feature"] == "P1"]) == 2
	assert len([ x for x in prosites if x["feature"] == "P2"]) == 1


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

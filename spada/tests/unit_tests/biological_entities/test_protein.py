#!/usr/bin/env python

from spada.biological_entities import protein

import pytest

info = {	"gene_id":			"A",
			"exonStructure": 	[[1,30],[70,100],[120,150],[170,190]],
			"txCoords": 		[1,190],
			"cdsCoords":		[21,172],
		  	"strand":			"+",
		  	"chr":				1,
			"proteinSequence":	"EFGHIJKABCDEFGHIJKLMNOPKR",
			"Pfam":				{"D1": [[1,3]], "D2": [[7,9]]},
			"Prosite":			{"P1": [[4,6],[7,9]]},
			"IDR":				{"ABCDEFG": [[8, 14]]} }

def test_init():

	p = protein.Protein("X", info)

	# basic information
	assert p._tx == "X"
	assert p._gene == info["gene_id"]
	assert p._sequence == info["proteinSequence"]
	assert p._pfam == info["Pfam"]
	assert p._prosite == info["Prosite"]
	assert p._idr == info["IDR"]

	# explore structure
	structure = [ x for x in p.structure_ordered ]
	assert len(structure) == len(info["proteinSequence"])
	assert structure[0].genomicPosition == 21
	assert structure[24].genomicPosition == 170
	assert "".join([ x.res for x in structure ]) == info["proteinSequence"]
	assert all( (y.num - x.num) == 1 for x,y in zip(structure,structure[1:]) )
	assert all( y.genomicPosition > x.genomicPosition for x,y in zip(structure,structure[1:]) )

	infominus = info
	infominus["strand"] = "-"
	pminus = protein.Protein("X", infominus)

	# explore structure if strand is minus
	structureminus = [ x for x in pminus.structure_ordered ]
	assert len(structureminus) == len(info["proteinSequence"])
	assert structureminus[0].genomicPosition == 172
	assert structureminus[24].genomicPosition == 23
	assert "".join([ x.res for x in structureminus ]) == info["proteinSequence"]
	assert all( (y.num - x.num) == 1 for x,y in zip(structureminus,structureminus[1:]) )
	assert all( y.genomicPosition < x.genomicPosition for x,y in zip(structureminus,structureminus[1:]) )

	# explore features
	assert sum([ "D1" in x._features for x in structure ]) == 3
	assert "D1" in structure[0]._features
	assert "D1" in structure[1]._features
	assert "D1" in structure[2]._features
	assert sum([ "D2" in x._features for x in structure ]) == 3
	assert "D2" in structure[6]._features
	assert "D2" in structure[7]._features
	assert "D2" in structure[8]._features
	assert sum([ "P1" in x._features for x in structure ]) == 6
	assert "P1" in structure[3]._features
	assert "P1" in structure[4]._features
	assert "P1" in structure[5]._features
	assert "P1" in structure[6]._features
	assert "P1" in structure[7]._features
	assert "P1" in structure[8]._features

def test_getFeature():

	p = protein.Protein("X", info)

	pfam1 = p.getFeature("Pfam", "D1")
	assert len(pfam1) == 1
	assert len(pfam1[0]) == 3
	assert "".join([ x.res for x in pfam1[0] ]) == "EFG"

	pfam2 = p.getFeature("Pfam", "D2")
	assert len(pfam2) == 1
	assert len(pfam2[0]) == 3
	assert "".join([ x.res for x in pfam2[0] ]) == "KAB"

	prosite = p.getFeature("Prosite", "P1")
	assert len(prosite) == 2
	assert len(prosite[0]) == 3
	assert len(prosite[1]) == 3
	assert "".join([ x.res for x in prosite[0] ]) == "HIJ"
	assert "".join([ x.res for x in prosite[1] ]) == "KAB"

	idrs = p.getFeature("IDR", "ABCDEFG")
	assert len(idrs) == 1
	assert len(idrs[0]) == 7
	assert "".join([ x.res for x in idrs[0] ]) == "ABCDEFG"

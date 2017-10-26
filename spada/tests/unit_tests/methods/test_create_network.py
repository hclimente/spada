#!/usr/bin/env python

from spada.methods import create_network
from spada import utils

import pandas as pd
import pytest
import os

def test_createNetworks():

	c = create_network.CreateNetwork("test", "gencode")

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../../data"
	gtf = dataPath + "/gtf"
	c.createNetworks(gtf)

	assert len(c._genes.nodes()) == 14
	assert len(c._txs.nodes()) == 21
	assert len(c._txs.nodes()["ENST00000450305.2"]["exonStructure"]) == 6
	assert c._txs.nodes()["ENST00000450305.2"]["strand"] == "+"
	assert c._txs.nodes()["ENST00000450305.2"]["chr"] == "chr1"
	assert not c._txs.nodes()["ENST00000450305.2"]["cdsCoords"]
	assert len(c._txs.nodes()["ENST00000494149.2"]["exonStructure"]) == 1
	assert c._txs.nodes()["ENST00000494149.2"]["strand"] == "-"
	assert c._txs.nodes()["ENST00000494149.2"]["chr"] == "chr1"
	assert not c._txs.nodes()["ENST00000494149.2"]["cdsCoords"]
	assert len(c._txs.nodes()["ENST00000335137.3"]["exonStructure"]) == 1
	assert c._txs.nodes()["ENST00000335137.3"]["strand"] == "+"
	assert c._txs.nodes()["ENST00000335137.3"]["chr"] == "chr1"
	assert c._txs.nodes()["ENST00000335137.3"]["cdsCoords"] == [69091,70005]

def test_readExpression():

    expression = pd.DataFrame({'sample1': [1., 2., 3.],
                               'sample2': [4., 5., 6.],
                               'sample3': [7., 8., 9.]},
                               index = pd.Series(["A", "B", "C"]))

    c = create_network.CreateNetwork("test", "gencode")

    assert c.readExpression(expression) == {"A": 4, "B": 5, "C": 6}

def test_calculatePSI():

    expression = pd.DataFrame({'sample1': [1., 2., 3.],
                               'sample2': [4., 5., 6.],
                               'sample3': [7., 8., 9.]},
                               index = pd.Series(["A", "B", "C"]))

    c = create_network.CreateNetwork("test", "gencode")

    c._txs.add_node("A", "1")
    c._txs.add_node("B", "1")
    c._txs.add_node("C", "2")

    assert c.calculatePSI(expression) == {"A": 1/3, "B": 2/3, "C": 1}

def test_isExpressed():

    expression = pd.DataFrame({'sample1': [1., 2., 3.],
                               'sample2': [4., 5., 6.],
                               'sample3': [7., 8., 9.]},
                               index = pd.Series(["A", "B", "C"]))

    c = create_network.CreateNetwork("test", "gencode")

    c._txs.add_node("A", "1")
    c._txs.add_node("B", "1")
    c._txs.add_node("C", "2")

    median = c.readExpression(expression)

    assert c.isExpressed(median, 3) == {"1": set(["A","B"]), "2": set("C")}
    assert c.isExpressed(median, 4.5) == {"1": set("B"), "2": set("C")}
    assert c.isExpressed(median, 7) == { }

def test_getInteractions():

	c = create_network.CreateNetwork("test", "gencode")
	c._genes.add_node(gene_id = "A", gene_symbol = "HSPB8")
	c._genes.add_node(gene_id = "B", gene_symbol = "HSPB7")
	c._genes.add_node(gene_id = "C", gene_symbol = "RPIA")
	c._genes.add_node(gene_id = "D", gene_symbol = "WDYHV1")
	c._genes.add_node(gene_id = "E", gene_symbol = "KRT15")
	c._genes.add_node(gene_id = "F", gene_symbol = "KRT20")
	c._genes.add_node(gene_id = "G", gene_symbol = "TPM1")
	c._genes.add_node(gene_id = "H", gene_symbol = "KXD1")

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../../data"
	mitab = dataPath + "/ppis"
	c.getInteractions(mitab)

	assert c._genes._net.has_edge("A", "B")
	assert c._genes._net.has_edge("C", "D")
	assert c._genes._net.has_edge("E", "F")
	assert c._genes._net.has_edge("G", "H")
	assert len(c._genes._net.edges()) == 4

def test_readDrivers():

	c = create_network.CreateNetwork("test", "gencode")

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../../data"
	drivers = dataPath + "/drivers"
	drivers, specificDrivers = c.readDrivers(drivers)

	assert len(drivers) == 3
	assert len(specificDrivers) == 2
	assert "CICP27" in drivers
	assert "RP11-34P13.3" in specificDrivers
	assert "WASH7P" in specificDrivers

def test_getIsoformSequences():

	c = create_network.CreateNetwork("test", "gencode")
	c._txs.add_node("ENST00000335137.3", "1")
	c._txs.add_node("ENST00000417324.1", "1")
	c._txs.add_node("ENST00000442987.3", "2")
	c._txs.add_node("ENST00000450305.2", "2")
	c._txs.add_node("ENST00000453576.2", "3")
	c._txs.add_node("ENST00000456328.2", "3")
	c._txs.add_node("ENST00000461467.1", "3")
	c._txs.add_node("ENST00000466430.5", "4")

	assert not c._txs.nodes()["ENST00000335137.3"]["proteinSequence"]
	assert not c._txs.nodes()["ENST00000417324.1"]["proteinSequence"]
	assert not c._txs.nodes()["ENST00000442987.3"]["proteinSequence"]
	assert not c._txs.nodes()["ENST00000450305.2"]["proteinSequence"]
	assert not c._txs.nodes()["ENST00000453576.2"]["proteinSequence"]
	assert not c._txs.nodes()["ENST00000456328.2"]["proteinSequence"]
	assert not c._txs.nodes()["ENST00000461467.1"]["proteinSequence"]
	assert not c._txs.nodes()["ENST00000466430.5"]["proteinSequence"]

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../../data"
	fasta = dataPath + "/fasta"
	c.getIsoformSequences(fasta)

	assert c._txs.nodes()["ENST00000335137.3"]["proteinSequence"] == "ASDFAFAFA"
	assert c._txs.nodes()["ENST00000417324.1"]["proteinSequence"] == "ASDASDASD"
	assert c._txs.nodes()["ENST00000442987.3"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST00000450305.2"]["proteinSequence"] == "ASFASFAS"
	assert c._txs.nodes()["ENST00000453576.2"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST00000456328.2"]["proteinSequence"] == "ASFSAFASFS"
	assert c._txs.nodes()["ENST00000461467.1"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST00000466430.5"]["proteinSequence"] == "ASFASFASASFASF"

def test_getIsoformFeatures():

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../../data"
	features = dataPath + "/features"

	c = create_network.CreateNetwork("test", "gencode")
	c._txs.add_node("ENST00000335137.3", "1")
	c._txs.add_node("ENST00000417324.1", "1")
	c._txs.add_node("ENST00000442987.3", "2")

	assert not c._txs.nodes()["ENST00000335137.3"]["Pfam"]
	assert not c._txs.nodes()["ENST00000417324.1"]["IDR"]
	assert not c._txs.nodes()["ENST00000442987.3"]["Prosite"]

	c.getIsoformFeatures(features)

	assert c._txs.nodes()["ENST00000335137.3"]["Pfam"]["xyz"] == {(3,6), (40,93)}
	assert c._txs.nodes()["ENST00000335137.3"]["Pfam"]["kile"] == {(40,93)}
	assert len(c._txs.nodes()["ENST00000335137.3"]["Pfam"]) == 2
	assert c._txs.nodes()["ENST00000417324.1"]["IDR"]["ASDASD"] == {(4,13)}
	assert c._txs.nodes()["ENST00000442987.3"]["Prosite"]["aasd"] == {(23,123)}

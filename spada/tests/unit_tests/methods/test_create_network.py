#!/usr/bin/env python

from spada.methods import create_network
from spada.utils import SpadaError

import networkx as nx
import os
import pandas as pd
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_createNetworks():

	c = create_network.CreateNetwork("test", "gencode")
	c.createNetworks(dataPath + "gtf")

	assert len(c._genes.nodes()) == 14
	assert len(c._txs.nodes()) == 21
	assert len(c._txs.nodes()["ENST02.2"]["exonStructure"]) == 6
	assert c._txs.nodes()["ENST02.2"]["strand"] == "+"
	assert c._txs.nodes()["ENST02.2"]["chr"] == "chr1"
	assert c._txs.nodes()["ENST02.2"]["cdsCoords"] == [13000, 13027]
	assert len(c._txs.nodes()["ENST19.2"]["exonStructure"]) == 1
	assert c._txs.nodes()["ENST19.2"]["strand"] == "-"
	assert c._txs.nodes()["ENST19.2"]["chr"] == "chr1"
	assert not c._txs.nodes()["ENST19.2"]["cdsCoords"]
	assert len(c._txs.nodes()["ENST12.3"]["exonStructure"]) == 1
	assert c._txs.nodes()["ENST12.3"]["strand"] == "+"
	assert c._txs.nodes()["ENST12.3"]["chr"] == "chr1"
	assert c._txs.nodes()["ENST12.3"]["cdsCoords"] == [69091,70021]

	assert os.stat("genes.pkl").st_size > 0
	assert os.stat("transcripts.pkl").st_size > 0
	os.remove("genes.pkl")
	os.remove("transcripts.pkl")

	with pytest.raises(SpadaError):
		c = create_network.CreateNetwork("test", "gencode", new = False)
		c.createNetworks(dataPath + "gtf")

def test_measureExpression():

	with pytest.raises(SpadaError):
		c = create_network.CreateNetwork("test", "gencode")
		c.measureExpression(None, -3.3, "N")

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

	c.getInteractions(dataPath + "ppis")

	assert c._genes._net.has_edge("A", "B")
	assert c._genes._net.has_edge("C", "D")
	assert c._genes._net.has_edge("E", "F")
	assert c._genes._net.has_edge("G", "H")
	assert len(c._genes._net.edges()) == 4

	with pytest.raises(SpadaError):
		c = create_network.CreateNetwork("test", "gencode")
		c.getInteractions(None)

def test_getDomainInteractions():

	c = create_network.CreateNetwork("test", "gencode")

	c._genes.add_node(gene_id = "A", gene_symbol = "")
	c._genes.add_node(gene_id = "B", gene_symbol = "")
	c._genes.add_node(gene_id = "C", gene_symbol = "")
	c._genes.add_edge(gene_id1 = "A", gene_id2 = "B")
	c._genes.add_edge(gene_id1 = "A", gene_id2 = "C")

	c._txs.add_node("A1", "A")
	c._txs.update_node("A1", "Pfam", (0,0), "D1")
	c._txs.update_node("A1", "Pfam", (0,0), "D2")
	c._txs.update_node("A1", "Pfam", (0,0), "D6")
	c._txs.add_node("A2", "A")
	c._txs.update_node("A2", "Pfam", (0,0), "D1")
	c._txs.update_node("A2", "Pfam", (0,0), "D3")
	c._txs.add_node("B1", "B")
	c._txs.update_node("B1", "Pfam", (0,0), "D4")
	c._txs.add_node("C1", "C")
	c._txs.update_node("C1", "Pfam", (0,0), "D5")

	c.getDomainInteractions(dataPath + "ddis")

	assert c._txs._net.has_edge("A1", "B1")
	assert c._txs._net["A1"]["B1"]["ddi"] == {frozenset({"D2","D4"}), frozenset({"D6","D4"})}
	assert c._txs._net.has_edge("A2", "C1")
	assert c._txs._net["A2"]["C1"]["ddi"] == {frozenset({"D3","D5"})}
	assert not c._txs._net.has_edge("A1", "A2")
	assert not c._txs._net.has_edge("B1", "C1")

	with pytest.raises(SpadaError):
		c = create_network.CreateNetwork("test", "gencode")
		c.getDomainInteractions(None)

def test_readDrivers():

	c = create_network.CreateNetwork("test", "gencode")
	drivers, specificDrivers = c.readDrivers(dataPath + "drivers")

	assert len(drivers) == 3
	assert len(specificDrivers) == 2
	assert "GeneL" in drivers
	assert "GeneD" in specificDrivers
	assert "GeneB" in specificDrivers

def test_getIsoformSequences():

	c = create_network.CreateNetwork("test", "gencode")
	c._txs.add_node("ENST12.3", "1")
	c._txs.add_node("ENST08.1", "1")
	c._txs.add_node("ENST18.3", "2")
	c._txs.add_node("ENST02.2", "2")
	c._txs.add_node("ENST16.2", "3")
	c._txs.add_node("ENST01.2", "3")
	c._txs.add_node("ENST09.1", "3")
	c._txs.add_node("ENST13.5", "4")

	assert not c._txs.nodes()["ENST12.3"]["proteinSequence"]
	assert not c._txs.nodes()["ENST08.1"]["proteinSequence"]
	assert not c._txs.nodes()["ENST18.3"]["proteinSequence"]
	assert not c._txs.nodes()["ENST02.2"]["proteinSequence"]
	assert not c._txs.nodes()["ENST16.2"]["proteinSequence"]
	assert not c._txs.nodes()["ENST01.2"]["proteinSequence"]
	assert not c._txs.nodes()["ENST09.1"]["proteinSequence"]
	assert not c._txs.nodes()["ENST13.5"]["proteinSequence"]

	c.getIsoformSequences(dataPath + "fasta")

	assert c._txs.nodes()["ENST12.3"]["proteinSequence"] == "ASDFAFAFA"
	assert c._txs.nodes()["ENST08.1"]["proteinSequence"] == "ASDASDASD"
	assert c._txs.nodes()["ENST18.3"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST02.2"]["proteinSequence"] == "ASFASFAS"
	assert c._txs.nodes()["ENST16.2"]["proteinSequence"] == "ABCDEFGHI"
	assert c._txs.nodes()["ENST01.2"]["proteinSequence"] == "ASFSAFASFS"
	assert c._txs.nodes()["ENST09.1"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST13.5"]["proteinSequence"] == "ASFASFASASFASF"

	with pytest.raises(SpadaError):
		c = create_network.CreateNetwork("test", "gencode")
		c.getIsoformSequences(None)

def test_getIsoformFeatures():

	c = create_network.CreateNetwork("test", "gencode")
	c._txs.add_node("ENST01.2", "1")
	c._txs.add_node("ENST02.2", "1")
	c._txs.add_node("ENST20.1", "1")
	c._txs.add_node("ENST08.1", "1")
	c._txs.add_node("ENST18.3", "2")

	assert not c._txs.nodes()["ENST20.1"]["Pfam"]
	assert not c._txs.nodes()["ENST08.1"]["IDR"]
	assert not c._txs.nodes()["ENST18.3"]["Prosite"]

	c.getIsoformFeatures(dataPath + "features")

	assert c._txs.nodes()["ENST01.2"]["Pfam"]["D1"] == {(1,2), (3,4)}
	assert len(c._txs.nodes()["ENST01.2"]["Pfam"]) == 1
	assert c._txs.nodes()["ENST02.2"]["Pfam"]["D1"] == {(1,2)}
	assert c._txs.nodes()["ENST02.2"]["Pfam"]["D2"] == {(3,4)}
	assert c._txs.nodes()["ENST02.2"]["Pfam"]["D4"] == {(5,9)}
	assert len(c._txs.nodes()["ENST02.2"]["Pfam"]) == 3
	assert c._txs.nodes()["ENST20.1"]["Pfam"]["D2"] == {(3,6)}
	assert len(c._txs.nodes()["ENST20.1"]["Pfam"]) == 1
	assert c._txs.nodes()["ENST08.1"]["IDR"]["I1"] == {(4,13)}
	assert c._txs.nodes()["ENST18.3"]["Prosite"]["P1"] == {(23,123)}

	with pytest.raises(SpadaError):
		c = create_network.CreateNetwork("test", "gencode")
		c.getIsoformFeatures(None)

def test_CreateNetwork():

	with pytest.raises(SpadaError):
		c = create_network.CreateNetwork("test", "madeup")

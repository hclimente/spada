#!/usr/bin/env python

from spada.methods import create_network
from spada.io.error import SpadaError

import networkx as nx
import pytest
import os

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_init():

	c = create_network.CreateNetwork("test", "gencode")

	assert type(c._genes).__name__ == 'GENCODEGeneNetwork'
	assert type(c._txs).__name__ == 'GENCODETranscriptNetwork'

	c = create_network.CreateNetwork("test", "ucsc")
	assert type(c._genes).__name__ == 'UCSCGeneNetwork'
	assert type(c._txs).__name__ == 'UCSCTranscriptNetwork'

	with pytest.raises(SpadaError):
		create_network.CreateNetwork("test", "Unexistant")

def test_createNetworks():

	c = create_network.CreateNetwork("test", "gencode")
	c.createNetworks(dataPath + "gtf")

	assert len(c._genes.nodes()) == 16
	assert len(c._txs.nodes()) == 24
	assert len(c._txs.nodes()["ENST02.2"]["exons"]) == 6
	assert c._txs.nodes()["ENST02.2"]["strand"] == "+"
	assert c._txs.nodes()["ENST02.2"]["chr"] == "chr1"
	assert c._txs.nodes()["ENST02.2"]["CDS"] == [13000, 13023]
	assert len(c._txs.nodes()["ENST19.2"]["exons"]) == 1
	assert c._txs.nodes()["ENST19.2"]["strand"] == "-"
	assert c._txs.nodes()["ENST19.2"]["chr"] == "chr1"
	assert not c._txs.nodes()["ENST19.2"]["CDS"]
	assert not c._txs.nodes()["ENST19.2"]["main"]
	assert len(c._txs.nodes()["ENST12.3"]["exons"]) == 1
	assert c._txs.nodes()["ENST12.3"]["strand"] == "+"
	assert c._txs.nodes()["ENST12.3"]["chr"] == "chr1"
	assert c._txs.nodes()["ENST12.3"]["CDS"] == [69978, 70004]
	assert c._txs.nodes()["ENST12.3"]["main"]
	assert len([ x for x,i in c._txs.nodes(data=True) if i['main'] ]) == 3
	assert "ENSG999" in c._genes.nodes()
	assert "ENSG999" in c._genes.nodes()
	assert "ENST99.1" in c._txs.nodes()
	assert "ENST99.2" in c._txs.nodes()
	assert "ENST99.3" in c._txs.nodes()

	# test rejected genes, transcript and CDS
	assert "bad_status" not in c._genes.nodes()
	assert "bad_accepted_tag" not in c._txs.nodes()
	assert "bad_status" not in c._txs.nodes()
	assert "rejected_status" not in c._txs.nodes()
	assert "test" in c._genes.nodes()
	assert "test.1" in c._txs.nodes()
	assert c._txs.nodes()["test.1"]["CDS"] == None
	assert "ENST99.4" not in c._txs.nodes()
	assert "ENST99.5" not in c._txs.nodes()
	assert "ENST99.NA" not in c._txs.nodes()

	with pytest.raises(SpadaError):
		c = create_network.CreateNetwork("test", dataPath + 'annotation.pklz', new = False)
		c.createNetworks(dataPath + "gtf")

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

def test_gecaseIsoformSequences():

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

	c.gecaseIsoformSequences(dataPath + "fasta")

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
		c.gecaseIsoformSequences(None)

def test_gecaseIsoformFeatures():

	c = create_network.CreateNetwork("test", "gencode")
	c._txs.add_node("ENST01.2", "1")
	c._txs.add_node("ENST02.2", "1")
	c._txs.add_node("ENST20.1", "1")
	c._txs.add_node("ENST08.1", "1")
	c._txs.add_node("ENST18.3", "2")

	assert not c._txs.nodes()["ENST20.1"]["Pfam"]
	assert not c._txs.nodes()["ENST08.1"]["IDR"]
	assert not c._txs.nodes()["ENST18.3"]["Prosite"]

	c.gecaseIsoformFeatures(dataPath + "features")

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
		c.gecaseIsoformFeatures(None)

def test_addAberrant():

	c = create_network.CreateNetwork("test", "gencode")
	c.addAberrant(dataPath + 'aberrant')

	assert c._txs.nodes()["ABC.1"]["gene_id"] == "ENSG08.4"
	assert c._txs.nodes()["DE_FG_HI.2"]["gene_id"] == "ENSG00.5"
	assert c._txs.nodes()["JKLM-.3"]["gene_id"] == "ENSG00.5"

def test_CreateNetwork():

	with pytest.raises(SpadaError):
		c = create_network.CreateNetwork("test", "madeup")

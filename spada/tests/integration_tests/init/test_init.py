from spada.methods import create_network

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_init():

	c = create_network.CreateNetwork("test", "gencode")
	c.run(dataPath + "gtf",
		  dataPath + "expression",
		  dataPath + "expression",
		  -3.3,
		  dataPath + "fasta",
		  dataPath + "mitab",
		  dataPath + "ddis",
		  dataPath + "drivers",
		  dataPath + "features",
		  dataPath + "aberrant")

	# gtf
	assert len(c._genes.nodes()) == 14
	assert len(c._txs.nodes()) == 24 # 21 txs + 3 aberrant
	assert len(c._txs.nodes()["ENST02.2"]["exonStructure"]) == 6
	assert c._txs.nodes()["ENST02.2"]["strand"] == "+"
	assert c._txs.nodes()["ENST02.2"]["chr"] == "chr1"
	assert c._txs.nodes()["ENST02.2"]["cdsCoords"] == [13000, 13023]
	assert len(c._txs.nodes()["ENST19.2"]["exonStructure"]) == 1
	assert c._txs.nodes()["ENST19.2"]["strand"] == "-"
	assert c._txs.nodes()["ENST19.2"]["chr"] == "chr1"
	assert not c._txs.nodes()["ENST19.2"]["cdsCoords"]
	assert len(c._txs.nodes()["ENST12.3"]["exonStructure"]) == 1
	assert c._txs.nodes()["ENST12.3"]["strand"] == "+"
	assert c._txs.nodes()["ENST12.3"]["chr"] == "chr1"
	assert c._txs.nodes()["ENST12.3"]["cdsCoords"] == [69978, 70004]

	# aberrant
	assert c._txs.nodes()["ABC.1"]["gene_id"] == "ENSG08.4"
	assert c._txs.nodes()["DE_FG_HI.2"]["gene_id"] == "ENSG00.5"
	assert c._txs.nodes()["JKLM-.3"]["gene_id"] == "ENSG00.5"

	# interactions
	assert c._genes._net.number_of_edges() == 2
	assert c._genes._net.has_edge("ENSG00.5", "ENSG00.5")
	assert c._genes._net.has_edge("ENSG13.1", "ENSG00.5")

	# ddis
	assert c._txs._net.number_of_edges() == 4
	assert c._txs._net.has_edge("ENST01.2", "ENST01.2")
	assert c._txs._net.has_edge("ENST01.2", "ENST02.2")
	assert c._txs._net.has_edge("ENST02.2", "ENST02.2")
	assert c._txs._net.has_edge("ENST20.1", "ENST02.2")
	assert not c._txs._net.has_edge("ENST20.1", "ENST01.2")
	assert c._txs._net["ENST01.2"]["ENST01.2"]["ddi"] == {frozenset({"D1"})}
	assert c._txs._net["ENST01.2"]["ENST02.2"]["ddi"] == {frozenset({"D1"})}
	assert c._txs._net["ENST02.2"]["ENST02.2"]["ddi"] == {frozenset({"D1"}), frozenset({"D2"}), frozenset({'D4', 'D2'})}
	assert c._txs._net["ENST20.1"]["ENST02.2"]["ddi"] == {frozenset({"D2","D4"}), frozenset({"D2"})}

	# driver
	assert not c._genes.nodes()["ENSG08.4"]["driver"] and not c._genes.nodes()["ENSG08.4"]["specificDriver"]
	assert not c._genes.nodes()["ENSG09.6"]["driver"] and not c._genes.nodes()["ENSG09.6"]["specificDriver"]
	assert c._genes.nodes()["ENSG11.3"]["driver"] and not c._genes.nodes()["ENSG11.3"]["specificDriver"]
	assert c._genes.nodes()["ENSG03.3"]["driver"] and c._genes.nodes()["ENSG03.3"]["specificDriver"]
	assert c._genes.nodes()["ENSG01.1"]["driver"] and c._genes.nodes()["ENSG01.1"]["specificDriver"]

	# fasta
	assert c._txs.nodes()["ENST12.3"]["proteinSequence"] == "ASDFAFAFA"
	assert c._txs.nodes()["ENST08.1"]["proteinSequence"] == "ASDASDASD"
	assert c._txs.nodes()["ENST18.3"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST02.2"]["proteinSequence"] == "ASFASFAS"
	assert c._txs.nodes()["ENST16.2"]["proteinSequence"] == "ABCDEFGHI"
	assert c._txs.nodes()["ENST01.2"]["proteinSequence"] == "ASFSAFASFS"
	assert c._txs.nodes()["ENST09.1"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST13.5"]["proteinSequence"] == "ASFASFASASFASF"

	# features
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

from spada.methods import create_network

import os
import pytest

def test_init():

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../../data"

	c = create_network.CreateNetwork("test", "gencode")
	c.run("{}/gtf".format(dataPath),
		  "{}/expression".format(dataPath),
		  "{}/expression".format(dataPath),
		  -3.3,
		  "{}/fasta".format(dataPath),
		  "{}/mitab".format(dataPath),
		  "{}/ddis".format(dataPath),
		  "{}/drivers".format(dataPath),
		  "{}/features".format(dataPath))

	# gtf
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

	# interactions
	assert len(c._genes._net.edges()) == 0

	# ddis
	assert len(c._txs._net.edges()) == 0

	# driver
	assert not c._genes.nodes()["ENSG00000186092.4"]["driver"] and not c._genes.nodes()["ENSG00000186092.4"]["specificDriver"]
	assert not c._genes.nodes()["ENSG00000238009.6"]["driver"] and not c._genes.nodes()["ENSG00000238009.6"]["specificDriver"]
	assert c._genes.nodes()["ENSG00000233750.3"]["driver"] and not c._genes.nodes()["ENSG00000233750.3"]["specificDriver"]
	assert c._genes.nodes()["ENSG00000243485.3"]["driver"] and c._genes.nodes()["ENSG00000243485.3"]["specificDriver"]
	assert c._genes.nodes()["ENSG00000227232.5"]["driver"] and c._genes.nodes()["ENSG00000227232.5"]["specificDriver"]

	# fasta
	assert c._txs.nodes()["ENST00000335137.3"]["proteinSequence"] == "ASDFAFAFA"
	assert c._txs.nodes()["ENST00000417324.1"]["proteinSequence"] == "ASDASDASD"
	assert c._txs.nodes()["ENST00000442987.3"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST00000450305.2"]["proteinSequence"] == "ASFASFAS"
	assert c._txs.nodes()["ENST00000453576.2"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST00000456328.2"]["proteinSequence"] == "ASFSAFASFS"
	assert c._txs.nodes()["ENST00000461467.1"]["proteinSequence"] == "ASFASFASF"
	assert c._txs.nodes()["ENST00000466430.5"]["proteinSequence"] == "ASFASFASASFASF"

	# features
	assert c._txs.nodes()["ENST00000335137.3"]["Pfam"]["xyz"] == {(3,6), (40,93)}
	assert c._txs.nodes()["ENST00000335137.3"]["Pfam"]["kile"] == {(40,93)}
	assert len(c._txs.nodes()["ENST00000335137.3"]["Pfam"]) == 2
	assert c._txs.nodes()["ENST00000417324.1"]["IDR"]["ASDASD"] == {(4,13)}
	assert c._txs.nodes()["ENST00000442987.3"]["Prosite"]["aasd"] == {(23,123)}

from spada.network import transcript_network

import os
import gzip,pickle
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

gn,txs = pickle.load(gzip.open(dataPath + "annotation.pklz", "rb"))
txs.createLogger()

def test_add_node():

	assert 'kk' not in txs.nodes()
	assert txs.add_node('kk', 'PEO')
	assert txs.nodes()['kk']['gene_id'] == 'PEO'
	assert txs.add_node('kk', 'PEO')

def test_transcripts():

	assert len([ x for x in txs.transcripts() ]) > len([ x for x in txs.transcripts(onlyMain = True) ])

def test_NotImplementedError():

	bad_network = transcript_network.TranscriptNetwork('bad')

	with pytest.raises(NotImplementedError):
		bad_network.acceptCDS()

	with pytest.raises(NotImplementedError):
		bad_network.genenameFilter()

	with pytest.raises(NotImplementedError):
		bad_network.isMain()

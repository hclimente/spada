from spada.methods import method
from spada.network import gencode_transcript_network

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"
m = method.Method('test', dataPath + 'annotation.pklz')
txs = m._txs

def test_init():

	new = gencode_transcript_network.GENCODETranscriptNetwork('aName')

	assert new._name == 'aName'

def test_genenameFilter():

	assert txs.genenameFilter(full_name = 'id') == ('id', None)
	assert txs.genenameFilter(gene_symbol = 'id', gene_id = 'id') == (None, None)
	assert txs.genenameFilter(gene_symbol = 'locus', gene_id = 'locus') == (None, None)

def test_txFilter():

    assert txs.txFilter('made_up_tx') == None
    assert txs.txFilter('ENST16.2') == 'ENST16.2'
    assert txs.txFilter('ENST16') == 'ENST16.2'
    assert txs.txFilter('ENST17.1') == 'ENST17.1'
    assert txs.txFilter('ENST17') == 'ENST17.1'
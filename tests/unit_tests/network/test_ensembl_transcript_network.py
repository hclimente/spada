from spada.methods import method
from spada.network import ensembl_transcript_network

import os
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

txs = ensembl_transcript_network.ENSEMBLTranscriptNetwork('aName')
txs.skip_filter = False

def test_init():

	assert txs._name == 'aName'

def test_genenameFilter():

	assert txs.genenameFilter(full_name = 'id') == ('id', None)
	assert txs.genenameFilter(gene_symbol = 'id', gene_id = 'id') == (None, None)
	assert txs.genenameFilter(gene_symbol = 'locus', gene_id = 'locus') == (None, None)

def test_txFilter():

    assert txs.txFilter('made_up_tx') == None
    assert txs.txFilter('ENST16.2') == 'ENST16.2'
    txs.add_node('ENST16.2', 'test')
    txs.add_node('ENST17.1', 'test')
    assert txs.txFilter('ENST16.2') == 'ENST16.2'
    assert txs.txFilter('ENST16') == 'ENST16.2'
    assert txs.txFilter('ENST17.1') == 'ENST17.1'
    assert txs.txFilter('ENST17') == 'ENST17.1'
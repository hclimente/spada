from spada.network import ucsc_transcript_network

import os
import pickle
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data"

txs = ucsc_transcript_network.UCSCTranscriptNetwork('aName')

def test_init():

	assert txs._name == 'aName'

def test_genenameFilter():

	assert txs.genenameFilter(full_name = 'symbol|id') == ('id', 'symbol')
	assert txs.genenameFilter(gene_symbol = 'symbol', gene_id = 'id') == ('id', 'symbol')
	assert txs.genenameFilter(gene_symbol = 'locus', gene_id = 'locus') == (None, None)

from spada.network import ucsc_gene_network

import os
import pickle
import pytest

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data"

genes = ucsc_gene_network.UCSCGeneNetwork('aName')

def test_init():

	assert genes._name == 'aName'

def test_genenameFilter():

	assert genes.nameFilter(full_name = 'symbol|id') == ('id', 'symbol')
	assert genes.nameFilter(gene_symbol = 'symbol', gene_id = 'id') == ('id', 'symbol')
	assert genes.nameFilter(gene_symbol = 'locus', gene_id = 'locus') == (None, None)

	genes.add_node(gene_id = 'id1', gene_symbol = 'symbolA')

	assert genes.nameFilter(gene_symbol = 'symbolA') == ('id1', 'symbolA')


def test_add_node():

	assert not genes.add_node(gene_symbol = 'locus', gene_id = 'locus')

def test_update_node():

	assert not genes.update_node('fakeKey', 0, full_name = 'locus')

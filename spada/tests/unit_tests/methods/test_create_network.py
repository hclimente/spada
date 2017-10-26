#!/usr/bin/env python

from spada.methods import create_network
from spada import utils

import pandas as pd
import pytest
import os

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

	scriptPath = os.path.realpath(__file__)
	dataPath = os.path.dirname(scriptPath) + "/../../data"
	mitab = dataPath + "/ppis"

	c = create_network.CreateNetwork("test", "gencode")
	c._genes.add_node(gene_id = "A", gene_symbol = "HSPB8")
	c._genes.add_node(gene_id = "B", gene_symbol = "HSPB7")
	c._genes.add_node(gene_id = "C", gene_symbol = "RPIA")
	c._genes.add_node(gene_id = "D", gene_symbol = "WDYHV1")
	c._genes.add_node(gene_id = "E", gene_symbol = "KRT15")
	c._genes.add_node(gene_id = "F", gene_symbol = "KRT20")
	c._genes.add_node(gene_id = "G", gene_symbol = "TPM1")
	c._genes.add_node(gene_id = "H", gene_symbol = "KXD1")
	c.getInteractions(mitab)

	assert c._genes._net.has_edge("A", "B")
	assert c._genes._net.has_edge("C", "D")
	assert c._genes._net.has_edge("E", "F")
	assert c._genes._net.has_edge("G", "H")

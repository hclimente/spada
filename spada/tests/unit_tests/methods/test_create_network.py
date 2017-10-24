#!/usr/bin/env python

from spada.methods import create_network

import pandas as pd
import pytest

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

    assert c.calculatePSI(expression) == {"A": 4/9, "B": 5/9, "C": 1}

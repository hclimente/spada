# SPADA

[![Build Status](https://travis-ci.org/hclimente/spada.svg?branch=master)](https://travis-ci.org/hclimente/spada)
[![codecov](https://codecov.io/gh/hclimente/spada/branch/master/graph/badge.svg)](https://codecov.io/gh/hclimente/spada)

**WARNING:** *SPADA is still a work-in-progress. Take this README as a statement of intent.*

SPADA (Splicing-led Protein Alterations Discovered Agilely) is a tool to study the functional impact of alternative splicing changes between two conditions. The alterations have to be represented as isoform switches i.e. when one condition is best represented by one isoform, and the second condition by another one. Then, SPADA is able to predict which protein features are changing between both isoforms, and their impact on the protein-protein interaction network. A SPADA-based analysis of isoform switches found in the TCGA dataset was published in [Climente-González *et al.* (2017)](http://www.cell.com/cell-reports/abstract/S2211-1247(17)31104-X).

To start using SPADA, simply install it with

``` python
pip install spada
```

and run it in a toy dataset with

``` python
spada.py function --annotation gencode --switches spada.example
```

For information about how to run SPADA on your dataset, please read the [documentation](docs/index.md).

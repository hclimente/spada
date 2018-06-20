[![Build Status](https://travis-ci.org/hclimente/spada.svg?branch=master)](https://travis-ci.org/hclimente/spada)
[![codecov](https://codecov.io/gh/hclimente/spada/branch/master/graph/badge.svg)](https://codecov.io/gh/hclimente/spada)
[![PyPI version](https://badge.fury.io/py/spada.svg)](https://badge.fury.io/py/spada)

## spada

`spada` (Splicing-led Protein Alterations Discovered Agilely) is a tool to study the functional impact of alternative splicing changes between two conditions. The alterations have to be represented as isoform switches i.e. when one condition is best represented by one isoform, and the second condition by another one. Then, `spada` is able to predict which protein features are changing between both isoforms, and their impact on the protein-protein interaction network.

To start using `spada`, simply install it with

``` bash
pip install spada
```

### Documentation

#### Citation

If you use spada in a scientific publication, we would appreciate citations:

Climente-González, H., Porta-Pardo, E., Godzik, A., and Eyras, E. (2017). [The Functional Impact of Alternative Splicing in Cancer.]((http://www.cell.com/cell-reports/abstract/S2211-1247(17)31104-X)) Cell Rep. 20, 2215–2226.

#### Related projects

* [smartas](https://github.com/hclimente/smartas): Analysis of `spada` applied results applied to TCGA data.

[![Build Status](https://travis-ci.org/hclimente/spada.svg?branch=master)](https://travis-ci.org/hclimente/spada)
[![codecov](https://codecov.io/gh/hclimente/spada/branch/master/graph/badge.svg)](https://codecov.io/gh/hclimente/spada)
[![PyPI version](https://badge.fury.io/py/spada.svg)](https://badge.fury.io/py/spada)

## spada

spada (Splicing-led Protein Alterations Discovered Agilely) is a tool to study the functional impact of alternative splicing changes in case-control studies. The alterations are represented as isoform switches i.e. the situation when cases are best represented by one isoform, while controls are best represented by a different one. Then, spada is able to predict which protein features are changing between both isoforms, and their impact on the protein-protein interaction network.

To start using spada, simply install it with

``` bash
pip install spada
```

This package has two main functionalities:

- Calculating isoform switches from case-control transcript expression data:
  ```
  spada switches --expression-control expression_ctrl.tsv \
    --expression-case expression_case.tsv --minimum-expression 0.1
  ```
- Predicting the effects of splicing on the protein structure and the interactome of a set of switches:
  ```
  spada function --switches switches.tsv
  ```

For detailed information on how to use spada, please visit the [wiki](https://github.com/hclimente/spada/wiki).

### Documentation

#### Citation

If you use spada in a scientific publication, we would appreciate citations:

> Climente-González, H., Porta-Pardo, E., Godzik, A., and Eyras, E. (2017). [The Functional Impact of Alternative Splicing in Cancer.](http://www.cell.com/cell-reports/abstract/S2211-1247(17)31104-X) *Cell Reports 20*, 2215–2226.

#### Related projects

* [smartas](https://github.com/hclimente/smartas): we applied spada to transcriptomics data from the TCGA. This repository contains the code for the analyses and the main conclusions.

SmartAS
=======

The working directory must have the following structure:

* Data: contains the input data:

  * Sailfish output, in k-mer folders.
  * Human Interactome Project (HI_2012_PRE.tsv , downloaded from http://interactome.dfci.harvard.edu/H_sapiens/index.php)

* Pipeline: contains the code.

Dependencies:

* R v3.0.0
* Ensembl Perl API
* BioPerl 1.2.3
* Cytoscape v3.0.2, with the following Apps:
  * GeneMANIA (Cytoscape App) + Homo sapiens dataset
  * Jython Scripting Engine v2.5.2
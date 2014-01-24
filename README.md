SmartAS
=======

SmartAS is a pipeline oriented to finding interesting isoform switches between two conditions, with as many replicates as desired.

The working directory must have the following structure:

* Data: contains the input data:

	* Sailfish output, in k-mer folders.

* Pipeline: contains the code.

Dependencies:

* R v3.0.0
* Ensembl Perl API
* BioPerl 1.2.3
* Python 2.7 with the following modules:
	* Biana v1.3.1
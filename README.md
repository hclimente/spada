SmartAS
=======

SmartAS is a pipeline oriented to finding interesting isoform switches between two conditions, with as many replicates as desired.

The working directory must have the following structure:

* Data: contains the input data:

	* Sailfish output, in k-mer folders named as x-kmer-length. Inside, the sailfish output in the following format:
		x-kmer-length/<cell><compartment><replicate>_<kmer length>
	* GENCODE information:
		* annotation.gtf (gencode.v19.annotation.filter_pc.gtf)
		* proteins.fa (gencode.v19.pc_translations.fa)
		* transcripts.fa (gencode.v19.pc_transcripts.fa)

* Pipeline: contains the code.

Dependencies:

* R v3.0.0
* Python 2.7 with the following modules:
	* Biana v1.3.1
	* RPy2 2.2.2, with R 2.11.1 (this version is the "R" command that must be called from the shell).
	* Other libraries: subprocess, sys, getopt, shutil, os
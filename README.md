SmartAS
=======

SmartAS is a pipeline oriented to finding interesting isoform switches between two conditions, with as many replicates as desired.

The working directory must have the following structure:

* Data: contains the input data:
	* Databases
		* Using gene symbol
			* Druggable targets (dgidb_export_all_drivers_bygene_results.tsv)
			* Driver information (cancer_networks_SuppTables_v7_S7.csv)
			* AS Driver information (cancer_networks_SuppTables_v7_S6.csv)
			* Cancer drivers from COSMIC
			* compilationTable.tsv -> Endre/01_tcga_genome_annotation_full_noQuotes.txt
		* Using entrez id
			* Broad Institute gene sets:
				* Canonical pathways (c2.cp.v4.0.entrez.gmt)
				* GO biological process (c5.bp.v4.0.entrez.gmt)
				* Oncogenic signatures (c6.all.v4.0.entrez.gmt)
				* Immunologic signatures (c7.all.v4.0.entrez.gmt)
		* Using UniProt
			* Interactome3D

	* Precalculated files per protein
		* InterPro annotation
		* GPS phosphorylations
		* ProSite motifs
		* Fasta file containing as id the transcript id, gene fullname (entrez+symbol), the UniProt and the loops mapped with iLoops; the sequence is the annotated protein sequence (UnifiedFasta_iLoops13.fa)

	* Annotation files
		* TCGA transcript annotation file (knownGene.txt)
		* GTF with the transcript annotation

* testResults: output folder for SmartAS results.
* Pipeline: contains the code (optional).

Dependencies:

* R v3.0.0
* Python 2.7 with the following modules:
	* Biana v1.3.1
	* SBI libraries
	* GUILD
	* RPy2 2.2.2, with R 2.11.1 (this version is the "R" command that must be called from the shell).
	* Some other libraries and functions:
		* abc
		* argparse
		* base64
		* Bio.pairwise2
		* collections.Counter
		* cPickle
		* fisher
		* fnmatch
		* logging
		* loxun
		* math
		* networkx
		* numpy
		* os
		* pandas
		* re
		* scipy
		* SOAPpy
		* subprocess
		* sys
		* tarfile
		* time
		* transcript
		* xml.etree.ElementTree
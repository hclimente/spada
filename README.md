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

# Import your data #
In order to import and analyze your externally calculated switches, you need the following things:
* A list of your switches
* A parent SmartAS experiment to copy the networks from, with the same notation
* A list of expressed transcripts

The list of switches is a tab-delimited table with a header and a line per switch containing the gene name, the tumor isoform and the normal isoform. Its path is specified using the argument external-switches e.g. external-switches=my_switches.tsv. The table has the following format:

```
#!text

# One line header e.g. Gene	Tumor_isoform	NormalIsoform
Gene1	tumor_isoform	normal_isoform
Gene2	tumor_isoform	normal_isoform

```
The parent SmartAS experiment, identified by a tag, needs to be located under the same annotation directory in the working directory. It must me specified using the argument parent-tag e.g. parent-tag=brca. The new experiment will copy both the gene (geneNetwork.pkl) and the transcript (txNetwork.pkl) networks from its parent.
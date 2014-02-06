SmartAS
=======

SmartAS is a pipeline oriented to finding interesting isoform switches between two conditions, with as many replicates as desired.

The working directory must have the following structure:

* Data: contains the input data:

	├── Databases
	│   ├── compilationTable.tsv -> Endre/01_tcga_genome_annotation_full_noQuotes.txt
	│   ├── Endre
	│   │   ├── 01_tcga_genome_annotation_full_noQuotes.txt
	│   │   ├── 01_tcga_genome_annotation_full.txt
	│   │   ├── 02_baltz_gene_hgnc.txt
	│   │   ├── 03_castello_gene_hgnc.txt
	│   │   ├── 04_kwon_gene_hgnc.txt
	│   │   ├── 05_gonzalez_gene_hgnc.txt
	│   │   ├── 06_brosseau_gene_hgnc.txt
	│   │   ├── 07_vogelstein_gene_hgnc.txt
	│   │   ├── 08_han_gene_hgnc.txt
	│   │   ├── 09_juan_pr_list.txt
	│   │   ├── 10_juan_ap_list.txt
	│   │   ├── 11_biomart_gene_hgnc.txt
	│   │   └── 12_cosmic_gene_hgnc.txt
	│   └── Intogen.tsv
	├── GENCODE
	│   ├── annotation.gtf -> /home/hector/Data/GENCODE19/gencode.v19.annotation.gtf
	│   ├── Filtered
	│   │   ├── 10C1_20.filtered.sf
	│   │   ├── 10C2_20.filtered.sf
	│   │   ├── 7C1_20.filtered.sf
	│   │   └── 7C2_20.filtered.sf
	│   ├── proteins.fa -> /home/hector/Data/GENCODE19/gencode.v19.pc_translations.fa
	│   ├── Rawdata
	│   │   ├── 20-kmer-length -> /home/hector/Data/Sailfish/
	│   │   ├── 25-kmer-length -> /projects/rg/TCGA/users/gael/sailfish/25-kmer-length/
	│   │   └── 30-kmer-length -> /projects/rg/TCGA/users/gael/sailfish/30-kmer-length/
	│   └── transcripts.fa -> /home/hector/Data/GENCODE19/gencode.v19.pc_transcripts.fa
	└── TCGA
	    └── Rawdata -> /projects/rg/TCGA/pipeline/run6/

* Pipeline: contains the code.

Dependencies:

* R v3.0.0
* Python 2.7 with the following modules:
	* Biana v1.3.1
	* RPy2 2.2.2, with R 2.11.1 (this version is the "R" command that must be called from the shell).
	* Other libraries: subprocess, sys, getopt, shutil, os
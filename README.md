SmartAS
=======

SmartAS is a pipeline oriented to finding interesting isoform switches between two conditions, with as many replicates as desired.

The working directory must have the following structure:

* Data: contains the input data:
	* Databases: information from public Databases
		* c2.cp.v4.0.entrez.gmt
		* c5.bp.v4.0.entrez.gmt
		* c6.all.v4.0.entrez.gmt
		* compilationTable.tsv -> Endre/01_tcga_genome_annotation_full_noQuotes.txt
		* ids_all.lst
		* ids_varspl.lst
		* Interactome3D
		* Intogen.tsv
		* Uniprot.fasta
		* uniprot_sprot.fasta
		* uniprot_sprot_hs.fasta
		* uniprot_sprot_varsplic.fasta
		* uniprot_sprot_varsplic_hs.fasta
	* GENCODE
		* annotation.gtf -> /home/hector/Data/GENCODE19/gencode.v19.annotation.gtf
		* Filtered
		* proteins.fa -> /home/hector/Data/GENCODE19/gencode.v19.pc_translations.fa
		* Rawdata
		* transcripts.fa -> /home/hector/Data/GENCODE19/gencode.v19.pc_transcripts.fa
		* UnifiedFasta.fa
	* Input
	* TCGA
	    * annotation.gtf -> /projects/rg/TCGA/pipeline/run10/tcga_exon_annotation_full.gff
	    * External
	    * Filtered
	    * knownGene.txt
	    * proteins.fa -> /projects/rg/TCGA/users/hector/fastaFiles/proteins.fa
	    * Rawdata -> /projects/rg/TCGA/pipeline/run10
	    * specificDrivers -> /projects/rg/TCGA/users/hector/specificDrivers/
	    * UnifiedFasta_iLoops13.fa -> /projects/rg/TCGA/users/hector/fastaFiles/UnifiedFasta_iLoops13.fa
	    * UnifiedFasta_iLoops13_loopFamilies.txt
	    * UnifiedFasta_iLoops13_noLoops.li
	    * UnifiedFasta_iLoops_devel.fa -> /projects/rg/TCGA/users/hector/fastaFiles/UnifiedFasta_iLoops_devel.fa
	    * UnifiedFasta_iLoops_devel_loopFamilies.txt
	    * UnifiedFasta_iLoops_devel_noLoops.li

* Results: output folder for SmartAS results.
* InterPro: output folder for InterPro motif analysis on isoforms.
* iLoops: output folder for iLoops results in protein-protein interactions.
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
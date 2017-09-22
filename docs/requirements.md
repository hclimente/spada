* annotation
  * knownGene.txt

* driverness
  * [Mutational_drivers_per_tumor_type.tsv](http://www.intogen.org/downloads)
  * [cancer_networks_SuppTables_v7_S7.csv]()
  * data/ucsc/drivers_ucsc_notation.txt
  * data/ucsc/asdrivers_ucsc_notation.txt

* drugabble: data/Databases/dgidb_export_all_drivers_bygene_results.tsv

* function
  * [BIOGRID-MV-Physical-3.4.152.tab2.txt](https://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.152/BIOGRID-MV-Physical-3.4.152.tab2.zip)
  * Eduard intx: interactions_found_more_than_three_times.txt & Switched_interactions_consensus.txt.
  * interpro/
  * anchor/
  * iupred/
  * prosite/
  * iLoops: sequences.uniprot.loops.fa & sequences.uniprot.loops.fa

* genesets: Genesets from the Molecular Signatures Database
  * [H: hallmark gene sets](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/h.all.v6.0.symbols.gmt)
  * [C2: Canonical pathways](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.v6.0.symbols.gmt)
  * C5: GO gene sets
    * [Biological process](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c5.bp.v6.0.symbols.gmt)
    * [Cellular compartment](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c5.cc.v6.0.symbols.gmt)
    * [Molecular function](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c5.mf.v6.0.symbols.gmt)
  * [C6: Oncogenic signatures](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c6.all.v6.0.symbols.gmt)
* seq
  * _exon_mutation-functional-count_full.txt
  * \*_gene_read_[tumor/paired]-filtered.txt
  * \*_iso_psi_[tumor/paired]-filtered.txt
  * \*_iso_tpm_[tumor/paired]-filtered.txt
  * WGS_505_TCGA

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

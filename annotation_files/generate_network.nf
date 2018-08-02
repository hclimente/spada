#!/usr/bin/env nextflow

// Help message
helpMessage = """
Generate annotation files for spada based on GENCODE annotation.
Usage:
  ./generate_network.nf --db gencode --v 28 [ --genome = GRCh37 ]

PARAMETERS
- db                   Database (gencode, ensembl).
- v                    Version of the database.
- genome               (Optional, default GRCh38) Version of the genome (GRCh38 or GRCh37).
- out                  (Optional, default .) Output directory.
- spada_dir            (Optional, default ..) Directory of spada.

OUTPUT
- \${db}_v\${v}.pklz   spada annotation file.
"""

if (params.help || params.v == null){
    log.info helpMessage
    exit 0
}

params.out = '.'
params.spada_dir = file('..')
params.domine = "$params.spada_dir/INTERACTION.txt"
params.genome = 'GRCh38'

gencode2ensembl = [ 28: 92, 27: 91, 26: 89, 25: 87, 24: 84,
                    23: 82, 22: 80, 21: 78, 20: 76, 19: 75,
                    18: 73, 17: 72, 16: 71, 15: 70, 14: 69,
                    13: 68, 12: 67 ]

if ( params.db == 'gencode' | params.db == 'ensembl' ) {
  if ( params.db == 'gencode' ) ENSEMBL_VERSION = gencode2ensembl[params.v]
  else if ( params.db == 'ensembl' ) ENSEMBL_VERSION = params.v

  feature_dbs = (ENSEMBL_VERSION > 88)? ['pfam','scanprosite'] : ['pfam','prosite']
}

GENOME_RELEASE = params.genome

// DOWNLOAD DB FILES
////////////////////////////////////////
process get_gtf {

  output:
    file 'gtf' into gtf

  script:
  if (params.db == 'gencode') template 'gencode/download_gtf.sh'
  else if (params.db == 'ensembl') template 'ensembl/download_gtf.sh'

}

process get_fasta {

  output:
    file 'fasta' into fasta, fasta_features, fasta_idr

  script:
  if (params.db == 'gencode') template 'gencode/download_fasta.sh'
  else if (params.db == 'ensembl') template 'ensembl/download_fasta.py'

}

process get_features {

  input:
    val FASTA from fasta_features
    each DB from feature_dbs

  output:
    file '*_features.tsv' into structured_features

  script:
  if (params.db == 'gencode' | params.db == 'ensembl' ) {
    if (ENSEMBL_VERSION > 78) template 'ensembl/download_features.py'
    else template 'computation/interpro.sh'
  }
}

process get_ppi {

  output:
    file 'mitab' into mitab

  """
  wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.mitab.zip
  unzip BIOGRID-MV-Physical-LATEST.mitab.zip 
  mv *.mitab.txt mitab
  """

}

domine = file(params.domine)

process get_ddi {

  input:
    file domine

  output:
    file "ddi.tsv" into ddi

  """
  # 3did
  wget https://3did.irbbarcelona.org/download/current/3did_flat.gz
  gunzip -c 3did_flat.gz >3did_flat

  grep ID 3did_flat | sed -E 's/^[^(]+//' | sed 's/[()]//g' | sed 's/@Pfam//g' | sed -E 's/\\.[0-9]+//g' >>ddi.tmp

  # DOMINE
  sed -E 's/\\|/\t/g' $domine | cut -f1,2 >>ddi.tmp

  echo -e "domain1\tdomain2" >ddi.tsv
  sort ddi.tmp | uniq >>ddi.tsv
  """

}

// GET ISOFORM IDRs
////////////////////////////////////////
process run_iupred {

  input:
    file(FASTA) from fasta_idr.splitFasta(file: true, by: 100)

  output:
    file "*.tsv" into iupred

  script:
  template 'computation/iupred.sh'

}

process collect_iupred {

  executor 'local'

  input:
    file "idr*.tsv" from iupred.collect()

  output:
    file "idr.tsv" into idr

  """
  cat idr* >idr.tsv
  """

}

process get_features {

  input:
    file '*tsv' from structured_features
    file idr

  output:
    file 'features' into features

  """
  echo -e "Transcript\tFeature_type\tFeature_id\tStart\tEnd" >features
  cat *tsv >>features
  """

}

// RUN SPADA
////////////////////////////////////////
process create_spada_annotation {

  publishDir "$params.out", overwrite: true, mode: "copy"
  time '4d'

  input:
    file gtf
    file mitab
    file ddi
    file fasta
    file features

  output:
    file "${params.db}_v${params.v}.pklz"

  """
  spada init --name ${params.db}_v${params.v} --new --gtf $gtf --annotation ${params.db} \
--ppi $mitab --ddi $ddi --seq $fasta --features $features
  mv annotation.pklz ${params.db}_v${params.v}.pklz
  """

}

#!/usr/bin/env nextflow

// Help message
helpMessage = """
Generate annotation files for spada based on GENCODE annotation.
Usage:
  ./get_gencode_annotation.nf --v 28 [ --genome = GRCh38 ]

PARAMETERS
- v                    Version of GENCODE (from 12 to 28).
- genome               (Optional, default GRCh38) Version of the genome (GRCh38 or GRCh37).
- out                  (Optional, default .) Output directory.
- spada_dir            (Optional, default ..) Directory of spada.

OUTPUT
- gtf, transcripts.fa, *_features.tsv
"""

if (params.help || params.v == null){
  log.info helpMessage
  exit 0
}

params.out = '.'
params.spada_dir = file('..')

// DOWNLOAD DB FILES
////////////////////////////////////////
gencode2ensembl = [ '28': '92', '27': '91', '26': '89', '25': '87', '24': '84',
                    '23': '82', '22': '80', '21': '78', '20': '76', '19': '75',
                    '18': '73', '17': '72', '16': '71', '15': '70', '14': '69',
                    '13': '68', '12': '67' ]

v = params.v
v_ens = gencode2ensembl["$v"]

params.genome = 'GRCh38'
tag1 = (params.genome == 'GRCh37') ? '/GRCh37_mapping' : ''
tag2 = (params.genome == 'GRCh37') ? 'lift37' : ''

// DOWNLOAD GTF
////////////////////////////////////////
process get_gtf {

  publishDir "$params.out", overwrite: true, mode: "copy"

  output:
    file "*gtf"

  """
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$v$tag1/gencode.v${v}${tag2}.annotation.gtf.gz
  gunzip -c *gtf.gz >gtf
  """

}

// DOWNLOAD FASTA
////////////////////////////////////////
process get_fasta {

  publishDir "$params.out", overwrite: true, mode: "copy"

  output:
    file 'fasta'

  """
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$v$tag1/gencode.v${v}${tag2}.pc_translations.fa.gz
  gunzip -c *fa.gz | sed -E 's/[^>|]+\\|//' | sed -E 's/\\|.+//' >fasta
  """

}

// DOWNLOAD ENSEMBL FEATURES
////////////////////////////////////////
get_ensembl = file("$params.spada_dir/annotation_files/get_ensembl_annotation.nf")

process get_features {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
   file get_ensembl

  output:
    file "*_features.tsv"

  """
  nextflow run $get_ensembl --v $v_ens --genome $params.genome
  """

}
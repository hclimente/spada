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

// DOWNLOAD DB FILES
////////////////////////////////////////
if (params.db == 'gencode') {
    get_annotation = file("$params.spada_dir/annotation_files/get_gencode_annotation.nf")
    ensembl_complete = ''
} else if (params.db == 'ensembl') {
    get_annotation = file("$params.spada_dir/annotation_files/get_ensembl_annotation.nf")
    ensembl_complete = '--complete'
}

process get_annotation_files {

    input:
        file get_annotation

    output:
        file 'gtf' into gtf
        file 'fasta' into fasta_idr, fasta
        file '*_features.tsv' into structured_features

    """
    nextflow run $get_annotation --v $params.v --genome $params.genome --spada_dir $params.spada_dir $ensembl_complete -profile cluster
    """

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

  // executor 'local'
  errorStrategy 'retry'
  maxRetries 3

  input:
    file "protein.fa" from fasta_idr.splitFasta(file: true)

  output:
    file "idr.tsv" into iupred

  """
  iupred2a.py protein.fa long | grep -v '#' | awk '\$3 > 0.5' | sed 's/\\t/ /g' >iupred.out
  tx=`head -n1 protein.fa | sed 's/>//'`

  touch idr.tsv
  start=""
  end=""
  seq=""
  while read p; do
    pos=`echo \$p | cut -f1 -d' '`
    res=`echo \$p | cut -f2 -d' '`

    if [ "\$start" == "" ]; then
      start="\$pos"
    elif [ "\$pos" -ne "\$((end + 1))" ]; then
      if [ 5 -lt "\$((end - start))" ]; then
        echo -e "\$tx\tIDR\t\$seq\t\$start\t\$end" >>idr.tsv
      fi
      seq="\$res"
      start="\$pos"
      end="\$pos"
    else
      end="\$pos"
      seq="\$seq\$res"
    fi
  done <iupred.out

  if [ 5 -lt "\$((end - start))" ]; then
    echo -e "\$tx\tIDR\t\$seq\t\$start\t\$end" >>idr.tsv
  fi
  """

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

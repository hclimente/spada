#!/usr/bin/env nextflow

// Help message
helpMessage = """
Generate annotation files for spada based on GENCODE annotation.
Usage:
  ./get_gencode_annotation.nf --v 28 [ --genome = GRCh38 ]

PARAMETERS
- v                    Version of GENCODE (from 12 to 28).
- genome               (Optional.) Version of the genome (GRCh38 or GRCh37).
                       By default, it uses GRCh38.

OUTPUT
- gencode_vXX.pklz     spada annotation file.
"""

// Show help when needed
if (params.help || params.v == null){
    log.info helpMessage
    exit 0
}

params.out = '.'
params.domine = 'INTERACTION.txt'

// DOWNLOAD DB FILES
////////////////////////////////////////
gencode2ensembl = [ '28': '92', '27': '91', '26': '89', '25': '87', '24': '84',
                    '23': '82', '22': '80', '21': '78', '20': '76', '19': '75',
                    '18': '73', '17': '72', '16': '71', '15': '70', '14': '69',
                    '13': '68', '12': '67' ]

ensemblUrls = [ '92': 'http://apr2018.archive.ensembl.org/biomart', '91': 'http://dec2017.archive.ensembl.org/biomart',
                '90': 'http://aug2017.archive.ensembl.org/biomart', '89': 'http://may2017.archive.ensembl.org/biomart',
                '88': 'http://mar2017.archive.ensembl.org/biomart', '87': 'http://dec2016.archive.ensembl.org/biomart',
                '86': 'http://oct2016.archive.ensembl.org/biomart', '85': 'http://jul2016.archive.ensembl.org/biomart',
                '84': 'http://mar2016.archive.ensembl.org/biomart', '83': 'http://dec2015.archive.ensembl.org/biomart',
                '82': 'http://sep2015.archive.ensembl.org/biomart', '81': 'http://jul2015.archive.ensembl.org/biomart',
                '80': 'http://may2015.archive.ensembl.org/biomart', '79': 'http://mar2015.archive.ensembl.org/biomart',
                '78': 'http://dec2014.archive.ensembl.org/biomart', '77': 'http://oct2014.archive.ensembl.org/biomart',
                '76': 'http://aug2014.archive.ensembl.org/biomart', '75': 'http://feb2014.archive.ensembl.org/biomart',
                '74': 'http://dec2013.archive.ensembl.org/biomart', '67': 'http://may2012.archive.ensembl.org/biomart' ]

v = params.v
v_ens = gencode2ensembl["$v"]

params.genome = 'GRCh38'
tag1 = (params.genome == 'GRCh37') ? '/GRCh37_mapping' : ''
tag2 = (params.genome == 'GRCh37') ? 'lift37' : ''

process get_gtf {

  output:
    file '*gtf' into gtf

  """
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$v$tag1/gencode.v${v}${tag2}.annotation.gtf.gz
  gunzip -c *gtf.gz
  """

}

process get_ppi {

  output:
    file 'BIOGRID-MV-Physical-*mitab.txt' into mitab

  """
  wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.mitab.zip
  unzip BIOGRID-MV-Physical-LATEST.mitab.zip
  """

}

process get_fasta {

  output:
    file "*.fa" into fasta_idr, fasta

  """
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$v$tag1/gencode.v${v}${tag2}.pc_translations.fa.gz
  gunzip -c *fa.gz
  sed -E 's/[^>|]+\\|//' gencode.v${v}${tag2}.pc_translations.fa | sed -E 's/\\|.+//' >transcripts.fa
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
  gunzip -c 3did_flat.gz

	grep ID 3did_flat | sed -E 's/^[^(]+//' | sed 's/[()]//g' | sed 's/@Pfam//g' | sed -E 's/\\.[0-9]+//g' >>ddi.tmp

	# DOMINE
	sed -E 's/\\|/\t/g' $domine | cut -f1,2 >>ddi.tmp

	echo -e "domain1\tdomain2" >ddi.tsv
	sort ddi.tmp | uniq >>ddi.tsv
	"""

}

// GET ISOFORM FEATURES
////////////////////////////////////////
feature_dbs = ['pfam','scanprosite']

process download_features {

  input:
    each db from feature_dbs
    val url from ensemblUrls["$v_ens"]

  output:
    file "${db}_features.tsv" into features

  """
  #!/usr/bin/env python

  from biomart import BiomartServer

  server = BiomartServer( "$url" )
  tags = { 'pfam': 'Pfam', 'scanprosite': 'Prosite' }

  ensembl = server.datasets['hsapiens_gene_ensembl']

  response = ensembl.search({
    'attributes': [ 'ensembl_transcript_id', 'transcript_version', '$db', '${db}_start', '${db}_end' ]
  })

  with open('${db}_features.tsv', 'w') as OUT:
    for line in response.iter_lines():
      line = line.decode('utf-8')
      t,v,f,s,e = line.split("\\t")
      if f:
        OUT.write('{}.{}\\t{}\\t{}\\t{}\\t{}\\n'. format(t, tags['$db'], v, f, s, e))
  """

}

process runIupred {

  errorStrategy 'retry'
  maxRetries 3

  input:
    file "protein.fa" from fasta_idr.splitFasta(file: true)

  output:
    file "idr.tsv" into iupred

  """
  iupred protein.fa long | grep -v '#' | sed -E 's/^ +//g' | sed -E 's/ +/ /g' | awk '\$3 > 0.5' >iupred.out
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

process collectIupred {

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
    file '*tsv' from features.collect()
    file idr

  output:
    file 'features' into all_features

  """
  echo -e "Transcript\tFeature_type\tFeature_id\tStart\tEnd" >features
  cat *tsv >>features
  """

}

// RUN SPADA
////////////////////////////////////////
process create_spada_annotation {

  publishDir "$params.out", overwrite: true, mode: "copy"

	input:
		file gtf
		file mitab
		file ddi
		file fasta
		file all_features

	output:
		file "gencode_v${$v}.pklz"

	"""
  spada init --name gencode_v${$v} --new --gtf $gtf --annotation gencode \
--ppi $ppi --ddi $ddi --seq $fasta --features $features
	mv annotation.pklz gencode_v${$v}.pklz
	"""

}

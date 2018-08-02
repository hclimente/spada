#!/usr/bin/env nextflow

// Help message
helpMessage = """
Generate annotation files for spada based on GENCODE annotation.
Usage:
  ./get_gencode_annotation.nf --v 28 [ --genome = GRCh38 ]

PARAMETERS
- v                    Version of ENSEMBL (from 74 to 92).
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
params.genome = 'GRCh38'
tag = (params.genome == 'GRCh38')? params.genome : "${params.genome}.${params.v}"

urls = [ '92': 'http://apr2018.archive.ensembl.org/biomart', '91': 'http://dec2017.archive.ensembl.org/biomart',
         '90': 'http://aug2017.archive.ensembl.org/biomart', '89': 'http://may2017.archive.ensembl.org/biomart',
         '88': 'http://mar2017.archive.ensembl.org/biomart', '87': 'http://dec2016.archive.ensembl.org/biomart',
         '86': 'http://oct2016.archive.ensembl.org/biomart', '85': 'http://jul2016.archive.ensembl.org/biomart',
         '84': 'http://mar2016.archive.ensembl.org/biomart', '83': 'http://dec2015.archive.ensembl.org/biomart',
         '82': 'http://sep2015.archive.ensembl.org/biomart', '81': 'http://jul2015.archive.ensembl.org/biomart',
         '80': 'http://may2015.archive.ensembl.org/biomart', '79': 'http://mar2015.archive.ensembl.org/biomart',
         '78': 'http://dec2014.archive.ensembl.org/biomart', '77': 'http://oct2014.archive.ensembl.org/biomart',
         '76': 'http://aug2014.archive.ensembl.org/biomart', '75': 'http://feb2014.archive.ensembl.org/biomart',
         '74': 'http://dec2013.archive.ensembl.org/biomart' ]

if (params.complete) {
// DOWNLOAD GTF
////////////////////////////////////////
process get_gtf {

  publishDir "$params.out", overwrite: true, mode: "copy"

  output:
    file 'gtf'

  """
  wget ftp://ftp.ensembl.org/pub/release-$params.v/gtf/homo_sapiens/Homo_sapiens.${params.genome}.${params.v}.gtf.gz
  gunzip -c *gtf.gz >gtf
  """

}

// DOWNLOAD FASTA
////////////////////////////////////////
process get_fasta {

  publishDir "$params.out", overwrite: true, mode: "copy"

  input:
    val url from urls["$params.v"]
    
  output:
    file 'fasta' into fasta

  """
  #!/usr/bin/env python

  from biomart import BiomartServer

  server = BiomartServer( "$url" )

  ensembl = server.datasets['hsapiens_gene_ensembl']

  response = ensembl.search({
    'attributes': [ 'ensembl_transcript_id', 'peptide' ]
  })

  with open('fasta', 'w') as OUT:
    for line in response.iter_lines():
      line = line.decode('utf-8')
      p,t = line.split("\\t")
      if p and p != 'Sequence unavailable':
        OUT.write('>{}\\n{}\\n'. format(t, p.replace('*', '')))
  """

}
}

// GET ISOFORM FEATURES
////////////////////////////////////////
if (params.v.toInteger() > 78) {
  
  feature_dbs = (params.v.toInteger() > 88)? ['pfam','scanprosite'] : ['pfam','prosite']

  process download_features {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
      each db from feature_dbs
      val url from urls["$params.v"]

    output:
      file "${db}_features.tsv"

    """
    #!/usr/bin/env python

    from biomart import BiomartServer

    server = BiomartServer( "$url" )
    tags = { 'pfam': 'Pfam', 'scanprosite': 'Prosite', 'prosite': 'Prosite' }

    ensembl = server.datasets['hsapiens_gene_ensembl']

    response = ensembl.search({
      'attributes': [ 'ensembl_transcript_id', 'transcript_version', '$db', '${db}_start', '${db}_end' ]
    })

    with open('${db}_features.tsv', 'w') as OUT:
      for line in response.iter_lines():
        line = line.decode('utf-8')
        t,f,v,s,e = line.split("\\t")
        if f:
          OUT.write('{}.{}\\t{}\\t{}\\t{}\\t{}\\n'. format(t, tags['$db'], v, f, s, e))
    """

  }

} else {

  process compute_features {

    cpus 4
    memory '15 GB'

    input:
      file 'proteins.fa' from fasta.splitFasta( file: true, by: 20 )

    output:
      file 'features' into features

    """
    interproscan.sh -i proteins.fa -appl Pfam,ProSitePatterns,ProSiteProfiles  -f TSV -o features -dp
    """

  }

  process join_features {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
      file 'features*' from features.collect()

    output:
      file 'pfam_features.tsv'
      file 'prosite_features.tsv'

    """
    cat features* >tmp
    cut -f1,4,5,7,8 tmp | sed 's/ProSiteProfiles/Prosite/' | sed 's/ProSitePatterns/Prosite/' >features.tsv
    grep Pfam features.tsv >pfam_features.tsv
    grep Prosite features.tsv >prosite_features.tsv
    """

  }
}


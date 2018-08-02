#!/usr/bin/env python
'''
Input variables:
    - DB    Database. Either Pfam or Prosite.
    - ENSEMBL_VERSION   Version of Ensembl. 

Output file:
    - ${DB}_features.tsv
'''

from biomart import BiomartServer

urls = { '92': 'http://apr2018.archive.ensembl.org/biomart', '91': 'http://dec2017.archive.ensembl.org/biomart',
         '90': 'http://aug2017.archive.ensembl.org/biomart', '89': 'http://may2017.archive.ensembl.org/biomart',
         '88': 'http://mar2017.archive.ensembl.org/biomart', '87': 'http://dec2016.archive.ensembl.org/biomart',
         '86': 'http://oct2016.archive.ensembl.org/biomart', '85': 'http://jul2016.archive.ensembl.org/biomart',
         '84': 'http://mar2016.archive.ensembl.org/biomart', '83': 'http://dec2015.archive.ensembl.org/biomart',
         '82': 'http://sep2015.archive.ensembl.org/biomart', '81': 'http://jul2015.archive.ensembl.org/biomart',
         '80': 'http://may2015.archive.ensembl.org/biomart', '79': 'http://mar2015.archive.ensembl.org/biomart',
         '78': 'http://dec2014.archive.ensembl.org/biomart', '77': 'http://oct2014.archive.ensembl.org/biomart',
         '76': 'http://aug2014.archive.ensembl.org/biomart', '75': 'http://feb2014.archive.ensembl.org/biomart',
         '74': 'http://dec2013.archive.ensembl.org/biomart' }

server = BiomartServer( urls['${ENSEMBL_VERSION}'] )
tags = { 'pfam': 'Pfam', 'scanprosite': 'Prosite', 'prosite': 'Prosite' }

ensembl = server.datasets['hsapiens_gene_ensembl']

response = ensembl.search({
    'attributes': [ 'ensembl_transcript_id', 'transcript_version', '${DB}', '${DB}_start', '${DB}_end' ]
})

with open('${DB}_features.tsv', 'w') as OUT:
    for line in response.iter_lines():
        line = line.decode('utf-8')
        t,f,v,s,e = line.split("\\t")
        if f:
            OUT.write('{}.{}\\t{}\\t{}\\t{}\\t{}\\n'. format(t, tags['${DB}'], v, f, s, e))
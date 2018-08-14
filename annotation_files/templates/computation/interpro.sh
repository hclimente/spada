# Input variables:
#    - FASTA    Path to the fasta file.
# Output file:
#    - features

interproscan.sh -i ${FASTA} -appl Pfam,ProSitePatterns,ProSiteProfiles  -f TSV -o interpro_out -dp
cut -f1,4,5,7,8 interpro_out | sed 's/ProSiteProfiles/Prosite/' | sed 's/ProSitePatterns/Prosite/' >interpro_features.tsv

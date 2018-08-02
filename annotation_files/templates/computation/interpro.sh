# Input variables:
#    - FASTA    Path to the fasta file.
# Output file:
#    - features

interproscan.sh -i ${FASTA} -appl Pfam,ProSitePatterns,ProSiteProfiles  -f TSV -o features -dp
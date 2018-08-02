# Input variables:
#    - GENCODE_VERSION
#    - TAG1
#    - TAG2
# Output file:
#    - fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}${TAG1}/gencode.v${GENCODE_VERSION}${TAG2}.pc_translations.fa.gz
gunzip -c *fa.gz | sed -E 's/[^>|]+\\|//' | sed -E 's/\\|.+//' >fasta
# Input variables:
#    - ENSEMBL_VERSION    Version of Ensembl. 
#    - GENOME_RELEASE     
# Output file:
#     - gtf

wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/gtf/homo_sapiens/Homo_sapiens.${GENOME_RELEASE}.${ENSEMBL_VERSION}.gtf.gz
gunzip -c *gtf.gz >gtf
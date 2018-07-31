from spada.io import io
from spada.network.ensembl_gene_network import ENSEMBLGeneNetwork

class GENCODEGeneNetwork(ENSEMBLGeneNetwork):
	def __init__(self, name):

		self._accepted_status = ['KNOWN','NOVEL','PUTATIVE','KNOWN_BY_PROJECTION']
		ENSEMBLGeneNetwork.__init__(self, name)

	def accept(self, line):

		try:
			return line['gene_status'] in self._accepted_status
		except KeyError:
			return True

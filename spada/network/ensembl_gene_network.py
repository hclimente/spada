from spada.io import io
from spada.network.gene_network import GeneNetwork

class ENSEMBLGeneNetwork(GeneNetwork):
	def __init__(self, name):

		self._accepted_status = ['KNOWN','NOVEL','PUTATIVE','KNOWN_BY_PROJECTION']
		GeneNetwork.__init__(self, name)

	def accept(self, line):
		return True

	def nameFilter(self, full_name="", gene_id="", gene_symbol=""):

		geneSymbol 	= None
		geneID 		= None

		if full_name:
			geneID 	= full_name
		if not gene_id and gene_symbol:
			self.update_idx()
			geneID = self.symbol2id.get(gene_symbol)
		else:
			if gene_id:
				geneID 		= gene_id
			if gene_symbol:
				geneSymbol 	= gene_symbol

		return (geneID, geneSymbol)

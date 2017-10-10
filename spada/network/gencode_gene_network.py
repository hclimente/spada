from spada import utils
from spada.network import gene_network

import abc

class GENCODEGeneNetwork(gene_network.GeneNetwork):
	def __init__(self):
		gene_network.GeneNetwork.__init__(self, __name__)

	def nameFilter(self, full_name="", gene_id="", gene_symbol=""):
		geneSymbol 	= None
		geneID 		= None

		if full_name:
			geneID 	= full_name

		return (geneID, geneSymbol)

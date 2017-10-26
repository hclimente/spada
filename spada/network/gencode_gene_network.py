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
		if not gene_id and gene_symbol:
			assumedGeneId = [ x for x,y in self.nodes(data=True) if y["symbol"] == gene_symbol ]

			if assumedGeneId:
				geneID = assumedGeneId[0]
		else:
			if gene_id:
				geneID 		= gene_id
			if gene_symbol:
				geneSymbol 	= gene_symbol


		return (geneID, geneSymbol)

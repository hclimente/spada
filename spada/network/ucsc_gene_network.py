from spada.io import io
from spada.network import gene_network

import abc

class UCSCGeneNetwork(gene_network.GeneNetwork):
	def __init__(self, name):
		gene_network.GeneNetwork.__init__(self, name)

	def nameFilter(self, full_name="", gene_id="", gene_symbol=""):
		geneSymbol 	= None
		geneID 		= None

		if full_name:
			nameComponents 	= full_name.split("|")

			if len(nameComponents) > 1:
				geneSymbol 	= nameComponents[0]
				geneID 		= nameComponents[1]

		elif "locus" not in gene_id and "locus" not in gene_symbol:
			if gene_id:
				geneID 		= gene_id
			if gene_symbol:
				geneSymbol 	= gene_symbol

		if not geneID:
			self.update_idx()
			geneID = self.symbol2id.get(gene_symbol)

		return (geneID, geneSymbol)

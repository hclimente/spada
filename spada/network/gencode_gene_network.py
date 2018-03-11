from spada.io import io
from spada.network.gene_network import GeneNetwork

class GENCODEGeneNetwork(GeneNetwork):
	def __init__(self):

		self._accepted_status = ['KNOWN','NOVEL','PUTATIVE','KNOWN_BY_PROJECTION']
		GeneNetwork.__init__(self, __name__)

	def accept(self, line):
		return line['gene_status'] in self._accepted_status

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

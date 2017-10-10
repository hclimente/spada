from spada.network import isoform_network

class GENCODEIsoformNetwork(isoform_network.IsoformNetwork):
	def __init__(self):
		isoform_network.IsoformNetwork.__init__(self, __name__)

	def genenameFilter(self, full_name="", gene_id="", gene_symbol=""):
		geneSymbol 	= None
		geneID 		= None

		if full_name:
			geneID 		= full_name

		return (geneID, geneSymbol)

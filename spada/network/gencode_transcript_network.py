from spada.network import transcript_network

class GENCODETranscriptNetwork(transcript_network.TranscriptNetwork):
	def __init__(self):
		transcript_network.TranscriptNetwork.__init__(self, __name__)

	def genenameFilter(self, full_name="", gene_id="", gene_symbol=""):
		geneSymbol 	= None
		geneID 		= None

		if full_name:
			geneID 		= full_name

		return (geneID, geneSymbol)

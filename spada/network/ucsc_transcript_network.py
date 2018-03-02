from spada.network import transcript_network

class UCSCTranscriptNetwork(transcript_network.TranscriptNetwork):
	def __init__(self):
		transcript_network.TranscriptNetwork.__init__(self, __name__)

	def genenameFilter(self, full_name="", gene_id="", gene_symbol=""):
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

		return (geneID, geneSymbol)

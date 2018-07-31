from spada.network.transcript_network import TranscriptNetwork

class ENSEMBLTranscriptNetwork(TranscriptNetwork):
	def __init__(self, name):

		TranscriptNetwork.__init__(self, name)

	def accept(self, line):
		return True

	def acceptCDS(self, line):
		return True

	def isMain(self, line):
		return True

	def genenameFilter(self, full_name="", gene_id="", gene_symbol=""):
		geneSymbol 	= None
		geneID 		= None

		if full_name:
			geneID 	= full_name

		return (geneID, geneSymbol)

	def txFilter(self, name):

		ids = dict([ (y[0], '.'.join(y)) for y in map(lambda x: x.split('.'), self.nodes()) ])
		transcript = name if name in self.nodes() else ids.get(name)

		return transcript

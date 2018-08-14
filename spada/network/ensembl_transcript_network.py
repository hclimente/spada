from spada.network.transcript_network import TranscriptNetwork

class ENSEMBLTranscriptNetwork(TranscriptNetwork):
	def __init__(self, name):

		TranscriptNetwork.__init__(self, name)

		self.usedIds = {}

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

		if '.' in name:
			transcript, version = name.split('.')
			if transcript not in self.usedIds.keys():
				self.usedIds[transcript] = name
			return name
		else:
			return self.usedIds.get(name)
from spada.network.transcript_network import TranscriptNetwork

class ENSEMBLTranscriptNetwork(TranscriptNetwork):
	def __init__(self, name):

		TranscriptNetwork.__init__(self, name)

		self.usedIds = {}
		self.skip_filter = True

		self._rejected = ['cds_end_NF', 'cds_start_NF']

	def accept(self, line):
		
		reject = bool([ t for t in line['tags'] if t in self._rejected ])

		return not reject

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

		if self.skip_filter or name in self.nodes():
			return name
		if '.' in name:
			transcript, version = name.split('.')
			if transcript not in self.usedIds.keys():
				self.usedIds[transcript] = name
			return name
		else:
			return self.usedIds.get(name)
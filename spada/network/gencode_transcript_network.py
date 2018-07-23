from spada.network.transcript_network import TranscriptNetwork

class GENCODETranscriptNetwork(TranscriptNetwork):
	def __init__(self, name):

		TranscriptNetwork.__init__(self, name)

		self._main_tx = ['appris_principal_1', 'appris_principal_2']
		self._accepted_status = ['KNOWN','NOVEL','PUTATIVE','KNOWN_BY_PROJECTION']
		self._accepted_support_level = ['1','2','3']
		self._accepted = ['basic', 'CCDS']
		self._rejected = ['cds_end_NF', 'cds_start_NF']

	def accept(self, line):
		accept = bool([ t for t in line['tags'] if t in self._accepted ])
		reject = bool([ t for t in line['tags'] if t in self._rejected ])
		
		try:
			status = line['transcript_status'] in self._accepted_status
		except KeyError:
			try:
				status = line['transcript_support_level'] in self._accepted_support_level
			except KeyError:
				status = True

		return accept and status and not reject

	def acceptCDS(self, line):
		return line['transcript_type'] != "nonsense_mediated_decay"

	def isMain(self, line):
		return bool([ x for x in line['tags'] if x in self._main_tx ])

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

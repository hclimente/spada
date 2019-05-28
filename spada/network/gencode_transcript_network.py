from spada.network.ensembl_transcript_network import ENSEMBLTranscriptNetwork

class GENCODETranscriptNetwork(ENSEMBLTranscriptNetwork):
	def __init__(self, name):

		ENSEMBLTranscriptNetwork.__init__(self, name)

		self._main_tx = ['appris_principal_1', 'appris_principal_2']
		self._accepted_status = ['KNOWN','NOVEL','PUTATIVE','KNOWN_BY_PROJECTION']
		self._accepted_support_level = ['1','2','3']
		self._accepted = ['basic', 'CCDS']
		self._rejected = ['cds_end_NF', 'cds_start_NF']

		self.skip_filter = False

	def accept(self, line):
		#accept = bool([ t for t in line['tags'] if t in self._accepted ])
		reject = bool([ t for t in line['tags'] if t in self._rejected ])
		
		try:
			valid_status = line['transcript_status'] in self._accepted_status
		except KeyError:
			try:
				valid_status = line['transcript_support_level'] in self._accepted_support_level
			except KeyError:
				valid_status = True

		return valid_status and not reject

	def acceptCDS(self, line):
		return line['transcript_type'] != "nonsense_mediated_decay"

	def isMain(self, line):
		return bool([ x for x in line['tags'] if x in self._main_tx ])

import protein
import transcript

class IsoformSwitch:
	def __init__(self, nTx, tTx, score, patients):
		self._normal_transcript_name 	= nTx
		self._tumor_transcript_name 	= tTx
		self._percent 					= score
		self._patients 					= patients

		self._normal_transcript 		= None
		self._tumor_transcript 			= None
		self._normal_protein 			= None
		self._tumor_protein 			= None
		
	@property
	def nTx(self): return self._normal_transcript
	@property
	def tTx(self): return self._tumor_transcript
	@property
	def score(self): return self._percent

	@property 
	def nIsoform(self): return self._normal_protein
	@property 
	def tIsoform(self): return self._tumor_protein
	@property
	def nTranscript(self): return self._normal_transcript
	@property
	def tTranscript(self): return self._tumor_transcript

	@property 
	def cds_overlap(self):
		if self._normal_transcript.cds_exclusive and self._normal_transcript.cds_exclusive < 1:
			return True
		elif self._tumor_transcript.cds_exclusive and self._tumor_transcript.cds_exclusive < 1:
			return True
		else:
			return False 

	@property 
	def cds_diff(self):
		cdsDiff = [ x for x in self._normal_transcript.cds if self._normal_transcript.cds[x] ]
		cdsDiff.extend( [ x for x in self._tumor_transcript.cds if self._tumor_transcript.cds[x] ])
		return cdsDiff

	@property 
	def utr_diff(self):
		utrDiff = [ x for x in self._normal_transcript.utr if self._normal_transcript.utr[x] ]
		utrDiff.extend( [ x for x in self._tumor_transcript.utr if self._tumor_transcript.utr[x] ])
		return utrDiff

	def addTxs(self, nInfo, tInfo):
		self._normal_transcript = transcript.Transcript( self._normal_transcript_name, nInfo )
		self._tumor_transcript 	= transcript.Transcript( self._tumor_transcript_name, tInfo )

		self.get_cdsDiff()
		self.get_utrDiff()

	def addIsos(self, nInfo, tInfo):
		self._normal_protein 	= protein.Protein( self._normal_transcript_name, nInfo)
		self._tumor_protein 	= protein.Protein( self._tumor_transcript_name, tInfo)

	def get_cdsDiff(self):
		"""Makes a list with the genomic position not shared between the CDS of two transcripts."""
		for gPos in self._normal_transcript.cds:
			if gPos not in self._tumor_transcript.cds: 	
				self._normal_transcript.cds[gPos] = True
			else: 										
				self._normal_transcript.cds[gPos] = False

		for gPos in self._tumor_transcript.cds:
			if gPos not in self._normal_transcript.cds:
				self._normal_transcript.cds[gPos] = True
			else:
				self._normal_transcript.cds[gPos] = False

	def get_utrDiff(self):
		"""Makes a list with the genomic position not shared between the UTR of two transcripts."""
		for gPos in self._normal_transcript.utr:
			if gPos not in self._tumor_transcript.utr:
				self._normal_transcript.utr[gPos] = True
			else:
				self._normal_transcript.utr[gPos] = False

		for gPos in self._tumor_transcript.utr:
			if gPos not in self._normal_transcript.utr:
				self._tumor_transcript.utr[gPos] = True
			else:
				self._tumor_transcript.utr[gPos] = False
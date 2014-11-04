import protein
import transcript

class IsoformSwitch:
	def __init__(self, nTx, tTx, score, patients, pval):
		self._normal_transcript_name= nTx
		self._tumor_transcript_name = tTx
		self._percent 				= score
		self._patients 				= patients
		self._p 					= pval

		self._normal_transcript 	= None
		self._tumor_transcript 		= None
		self._normal_protein 		= None
		self._tumor_protein 		= None

		#Structural analysis
		self._functional_change 	= None
		self._disorder_change 		= None
		self._broken_surfaces 		= None
		#Network analysis
		self._guild_top1 			= None
		self._guild_top5 			= None
		#Neighborhood analysis
		self._neighborhood_change 	= None
		#iLoops analysis
		self._iloops_change 		= None
		
	@property
	def nTx(self): return self._normal_transcript_name
	@property
	def tTx(self): return self._tumor_transcript_name
	@property
	def score(self): return self._percent
	@property
	def patients(self): return self._patients
	@property
	def p(self): return self._p
	
	@property 
	def nIsoform(self): return self._normal_protein
	@property 
	def tIsoform(self): return self._tumor_protein
	@property
	def nTranscript(self): return self._normal_transcript
	@property
	def tTranscript(self): return self._tumor_transcript

	@property
	def functionalChange(self): return self._functional_change
	@property
	def disorderChange(self): return self._disorder_change
	@property
	def brokenSurfaces(self): return self._broken_surfaces
	@property
	def guildTop1(self): return self._guild_top1
	@property
	def guildTop5(self): return self._guild_top5
	@property
	def neighborhoodChange(self): return self._neighborhood_change
	@property
	def iloopsChange(self): return self._iloops_change

	@property
	def is_relevant(self):
		if self.cds_overlap and (self.iloopsChange or self.brokenSurfaces or self.functionalChange):
			return True
		else:
			return False

	@property 
	def cds_overlap(self):
		"""Returns True if there is an overlap between the transcripts coding sequence."""
		#Check that is not None and thad the exclusive region is not the whole CDS.
		if self._normal_transcript.cds_exclusive and self._normal_transcript.cds_exclusive < 1:
			return True
		elif self._tumor_transcript.cds_exclusive and self._tumor_transcript.cds_exclusive < 1:
			return True
		else:
			return False 

	@property 
	def cds_diff(self):
		cdsDiff = [ x for x in self._normal_transcript.cds if x not in self._tumor_transcript.cds ]
		cdsDiff.extend( [ x for x in self._tumor_transcript.cds if x not in self._normal_transcript.cds ])
		return cdsDiff

	@property 
	def utr_diff(self):
		utrDiff = [ x for x in self._normal_transcript.utr if x not in self._tumor_transcript.utr ]
		utrDiff.extend( [ x for x in self._tumor_transcript.utr if x not in self._normal_transcript.utr ])
		return utrDiff

	def addTxs(self, nInfo, tInfo):
		self._normal_transcript = transcript.Transcript( self._normal_transcript_name, nInfo )
		self._tumor_transcript 	= transcript.Transcript( self._tumor_transcript_name, tInfo )

		self.get_cdsDiff()
		self.get_utrDiff()

	def addIsos(self, nInfo, tInfo):
		if nInfo["Uniprot"]:
			self._normal_protein = protein.Protein( self._normal_transcript_name, nInfo)
			self._normal_protein.checkInteractome3DStructures()
		if tInfo["Uniprot"]:
			self._tumor_protein  = protein.Protein( self._tumor_transcript_name, tInfo)
			self._tumor_protein.checkInteractome3DStructures()

		if self._normal_protein and self._tumor_protein:
			self.getAlteredRegions()

	def get_cdsDiff(self):
		"""Changes the values of the CDS dictionary of the transcripts to
			a bool, indicating if they are transcript specific or not."""
		for gPos in self._normal_transcript.cds:
			if gPos not in self._tumor_transcript.cds: 	
				self._normal_transcript._cds[gPos] = True
			else: 										
				self._normal_transcript._cds[gPos] = False

		for gPos in self._tumor_transcript.cds:
			if gPos not in self._normal_transcript.cds:
				self._normal_transcript._cds[gPos] = True
			else:
				self._normal_transcript._cds[gPos] = False

	def get_utrDiff(self):
		"""Changes the values of the UTR dictionary of the transcripts to
			a bool, indicating if they are transcript specific or not."""
		for gPos in self._normal_transcript.utr:
			if gPos not in self._tumor_transcript.utr:
				self._normal_transcript._utr[gPos] = True
			else:
				self._normal_transcript._utr[gPos] = False

		for gPos in self._tumor_transcript.utr:
			if gPos not in self._normal_transcript.utr:
				self._tumor_transcript._utr[gPos] = True
			else:
				self._tumor_transcript._utr[gPos] = False

	def getAlteredRegions(self):
		for res in self._normal_protein._structure:
			if res.genomicPosition not in [ y.genomicPosition for y in self._tumor_protein._structure]:
				res.setIsoformSpecific(True)

		for res in self._tumor_protein._structure:
			if res.genomicPosition not in [ y.genomicPosition for y in self._normal_protein._structure]:
				res.setIsoformSpecific(True)
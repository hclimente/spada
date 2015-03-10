from libs import options
from libs import utils
import protein
import transcript

import os

class IsoformSwitch:
	def __init__(self, nTx, tTx, score, patients,precision,sensitivity):
		self._normal_transcript_name= nTx
		self._tumor_transcript_name = tTx
		self._percent 				= score
		self._patients 				= patients
		self._precision 			= precision
		self._sensitivity 			= sensitivity

		self._normal_transcript 	= None
		self._tumor_transcript 		= None
		self._normal_protein 		= None
		self._tumor_protein 		= None

		#Relevance measure
		self._iloops_change		 	= None
		self._functional_change 	= None
		self._disorder_change 		= None
		self._anchor_change 		= None
		self._broken_surfaces 		= None
		self._ptm_change 			= None
		self._gps_change			= None
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
	def precision(self): return self._precision
	@property
	def sensitivity(self): return self._sensitivity
	
	@property 
	def nIsoform(self): return self._normal_protein
	@property 
	def tIsoform(self): return self._tumor_protein
	@property
	def nTranscript(self): return self._normal_transcript
	@property
	def tTranscript(self): return self._tumor_transcript

	@property
	def functionalChange(self): 
		if self._functional_change is None:
			self.readRelevanceAnalysis()
		return self._functional_change
	@property
	def disorderChange(self): 
		if self._disorder_change is None:
			self.readRelevanceAnalysis()
		return self._disorder_change
	@property
	def anchorChange(self): 
		if self._anchor_change is None:
			self.readRelevanceAnalysis()
		return self._anchor_change
	@property
	def brokenSurfaces(self): 
		if self._broken_surfaces is None:
			self.readRelevanceAnalysis()
		return self._broken_surfaces
	@property
	def gpsChange(self): 
		if self._gps_change is None:
			self.readRelevanceAnalysis()
		return self._gps_change
	@property
	def ptmChange(self): 
		if self._ptm_change is None:
			self.readRelevanceAnalysis()
		return self._ptm_change

	@property
	def guildTop1(self): return self._guild_top1
	@property
	def guildTop5(self): return self._guild_top5
	@property
	def neighborhoodChange(self): return self._neighborhood_change
	@property
	def iloopsChange(self): 
		if self._iloops_change is None:
			self.readRelevanceAnalysis()
		return self._iloops_change

	@property
	def is_relevant(self):
		"""We define as relevant a switch that:
			* Is considered significant in the statistical analysis.
			* There is an overlap in the CDS regions, so the switch
		 	  can be attributed to splicing.
			* Involves a change in the CDS or in the UTR.
			* Involves one change in:
				* Disordered regions.
				* Loops mapped with iLoops.
				* I3D Interaction surfaces altered.
				* Differences in mapped structural features.
		"""
		# cds_diff is required if there is any feature
		# utr_diff only is impossible if a feature change is required
		if self.cds_overlap and (self.disorderChange or self.iloopsChange or self.anchorChange or self.functionalChange or self.gpsChange or self.ptmChange):
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
		"""Returns a list with the differencen between the transcripts coding sequences."""
		cdsDiff = [ x for x in self._normal_transcript.cds if x not in self._tumor_transcript.cds ]
		cdsDiff.extend( [ x for x in self._tumor_transcript.cds if x not in self._normal_transcript.cds ])
		return cdsDiff

	@property 
	def utr_diff(self):
		"""Returns True if there is a difference between the transcripts utr sequences."""
		utrDiff = [ x for x in self._normal_transcript.utr if x not in self._tumor_transcript.utr ]
		utrDiff.extend( [ x for x in self._tumor_transcript.utr if x not in self._normal_transcript.utr ])
		return utrDiff

	def addTxs(self,nInfo,tInfo):
		"""Creates the transcript objects for the transcripts involved
		in the switch and calculates UTR and CDS differences."""
		self._normal_transcript = transcript.Transcript( self._normal_transcript_name, nInfo )
		self._tumor_transcript 	= transcript.Transcript( self._tumor_transcript_name, tInfo )

		self.get_cdsDiff()
		self.get_utrDiff()

	def addIsos(self,nInfo,tInfo,partialCreation=False):
		"""Creates the isoform objects for the transcripts involved
		in the switch if they have an UniProt identifier. If both do,
		it calculates the shared and specific regions."""
		if nInfo["Uniprot"]:
			self._normal_protein = protein.Protein( self._normal_transcript_name, nInfo)
			if not partialCreation:
				self._normal_protein.checkInteractome3DStructures()
		if tInfo["Uniprot"]:
			self._tumor_protein  = protein.Protein( self._tumor_transcript_name, tInfo)
			if not partialCreation:
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
				self._tumor_transcript._cds[gPos] = True
			else:
				self._tumor_transcript._cds[gPos] = False

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
		"""Calculates the specific and non-specific residues of the isoforms 
		involved in an isoform switch."""
		for res in self._normal_protein._structure:
			if res.genomicPosition not in [ y.genomicPosition for y in self._tumor_protein._structure]:
				res.setIsoformSpecific(True)

		for res in self._tumor_protein._structure:
			if res.genomicPosition not in [ y.genomicPosition for y in self._normal_protein._structure]:
				res.setIsoformSpecific(True)

	def readRelevanceAnalysis(self):
		if not os.path.exists("{0}structural_analysis/structural_summary{1}.tsv".format(options.Options().qout,options.Options().filetag)):
			raise Exception("Relevance information not generated.")
			return False

		for elements in utils.readTable("{0}structural_analysis/structural_summary{1}.tsv".format(options.Options().qout,options.Options().filetag)):
			if elements[1] == self.nTx and elements[2] == self.tTx:
				# QUITAR CUANDO SE ACTUALIZEN LOS RELEVANTES
				if elements[3] == "True": self._iloops_change = True
				elif elements[3] == "False": self._iloops_change = False

				if elements[4] == "True": self._broken_surfaces = True
				elif elements[4] == "False": self._broken_surfaces = False

				if elements[5] == "True": self._functional_change =True 
				elif elements[5] == "False": self._functional_change =False

				if elements[6] == "True": self._disorder_change = True
				elif elements[6] == "False": self._disorder_change = False

				# if elements[3] == "True": self._iloops_change = True
				# elif elements[3] == "False": self._iloops_change = False

				# if elements[4] == "True": self._functional_change =True 
				# elif elements[4] == "False": self._functional_change =False

				# if elements[5] == "True": self._disorder_change = True
				# elif elements[5] == "False": self._disorder_change = False

				# if elements[6] == "True": self._anchor_change = True
				# elif elements[6] == "False": self._anchor_change = False

				# if elements[7] == "True": self._ptm_change = True
				# elif elements[7] == "False": self._ptm_change = False

				# if elements[8] == "True": self._gps_change = True
				# elif elements[8] == "False": self._gps_change = False
from libs import options
from libs import utils
from biological_entities import protein
from biological_entities import transcript

import os
import operator

class IsoformSwitch:
	def __init__(self, nTx, tTx, patients):
		self._normal_transcript_name= nTx
		self._tumor_transcript_name = tTx
		self._patients 				= patients

		self._normal_transcript 	= None
		self._tumor_transcript 		= None
		self._normal_protein 		= None
		self._tumor_protein 		= None

		#Relevance measure
		self._iloops_change		 	= None
		self._domain_change 		= None
		self._disorder_change 		= None
		self._anchor_change 		= None
		self._broken_surfaces 		= None
		self._ptm_change 			= None

		#Relevance measure
		self._deep_domain_change 	= { }
		self._deep_disorder_change 	= { }
		self._deep_anchor_change 	= { }
		self._deep_ptm_change 		= { }

		#Neighborhood analysis
		self._neighborhood_change 	= None
		#iLoops analysis
		self._iloops_change 		= None
		
	@property
	def nTx(self): return self._normal_transcript_name
	@property
	def tTx(self): return self._tumor_transcript_name
	@property
	def patients(self): return self._patients
	
	@property 
	def nIsoform(self): return self._normal_protein
	@property 
	def tIsoform(self): return self._tumor_protein
	@property
	def nTranscript(self): return self._normal_transcript
	@property
	def tTranscript(self): return self._tumor_transcript

	@property
	def iloopsChange(self):
		if self._iloops_change is None:
			self.readRelevanceAnalysis()
		return self._iloops_change
	@property
	def domainChange(self): 
		if self._domain_change is None:
			self.readRelevanceAnalysis()
		return self._domain_change
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
	def ptmChange(self): 
		if self._ptm_change is None:
			self.readRelevanceAnalysis()
		return self._ptm_change

	@property
	def neighborhoodChange(self): return self._neighborhood_change

	@property
	def is_functional(self):
		"""We define as functional a switch that:
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
		if self.cds_overlap and (self.disorderChange or self.iloopsChange or self.anchorChange or self.domainChange or self.ptmChange):
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
		if nInfo["proteinSequence"]:
			self._normal_protein = protein.Protein( self._normal_transcript_name, nInfo)
			if not partialCreation:
				self._normal_protein.checkInteractome3DStructures()
		if tInfo["proteinSequence"]:
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
		randomTag = "_random" if self.patients==0.0 else ""
		if not os.path.exists("{0}structural_analysis/structural_summary{1}.tsv".format(options.Options().qout,randomTag)):
			raise Exception("Relevance information not generated.")
			return False

		for elements in utils.readTable("{0}structural_analysis/structural_summary{1}.tsv".format(options.Options().qout,randomTag)):
			if elements[2] == self.nTx and elements[3] == self.tTx:
				if elements[4] == "True": self._iloops_change = True
				elif elements[4] == "False": self._iloops_change = False

				if elements[5] == "True": self._domain_change = True 
				elif elements[5] == "False": self._domain_change = False

				if elements[6] == "True": self._disorder_change = True
				elif elements[6] == "False": self._disorder_change = False

				if elements[7] == "True": self._anchor_change = True
				elif elements[7] == "False": self._anchor_change = False

				if elements[8] == "True": self._ptm_change = True
				elif elements[8] == "False": self._ptm_change = False

	def readDeepRelevanceAnalysis(self,skipDomain=False,skipIupred=False,skipAnchor=False,skipPtm=False,filetag=""):
		if not self.is_functional:
			return

		if self.domainChange and not skipDomain:
			for line in utils.readTable("{}structural_analysis/interpro_analysis{}.tsv".format(options.Options().qout,filetag)):
				if line[2] == self.nTx and line[3] == self.tTx:
					self._deep_domain_change.setdefault(line[5].replace(" ","_"),[])
					if "Gained" in line[4]:
						self._deep_domain_change[line[5].replace(" ","_")].append("Gained_in_tumor")
					elif "Lost" in line[4]:
						self._deep_domain_change[line[5].replace(" ","_")].append("Lost_in_tumor")
					elif "Nothing" in line[4]:
						self._deep_domain_change[line[5].replace(" ","_")].append("Nothing")

		if self.disorderChange and not skipIupred:
			for line in utils.readTable("{}structural_analysis/iupred_analysis{}.tsv".format(options.Options().qout,filetag)):
				if line[2] == self.nTx and line[3] == self.tTx and float(line[-1]):
					self._deep_disorder_change.setdefault(line[5],[])
					if "Gained" in line[4]:
						self._deep_disorder_change[line[5]].append("Gained_in_tumor")
					elif "Lost" in line[4]:
						self._deep_disorder_change[line[5]].append("Lost_in_tumor")
					elif "Nothing" in line[4]:
						self._deep_disorder_change[line[5]].append("Nothing")
		
		if self.anchorChange and not skipAnchor:
			for line in utils.readTable("{}structural_analysis/anchor_analysis{}.tsv".format(options.Options().qout,filetag)):
				if line[2] == self.nTx and line[3] == self.tTx and float(line[-1]):
					self._deep_anchor_change.setdefault(line[5],[])
					if "Gained" in line[4]:
						self._deep_anchor_change[line[5]].append("Gained_in_tumor")
					elif "Lost" in line[4]:
						self._deep_anchor_change[line[5]].append("Lost_in_tumor")
					elif "Nothing" in line[4]:
						self._deep_anchor_change[line[5]].append("Nothing")
			
		if self.ptmChange and not skipPtm:
			for line in utils.readTable("{}structural_analysis/prosite_analysis{}.tsv".format(options.Options().qout,filetag)):
				if line[2] == self.nTx and line[3] == self.tTx:
					self._deep_ptm_change.setdefault(line[5],[])
					if "Gained" in line[4]:
						self._deep_ptm_change[line[5]].append("Gained_in_tumor")
					elif "Lost" in line[4]:
						self._deep_ptm_change[line[5]].append("Lost_in_tumor")
					elif "Nothing" in line[4]:
						self._deep_ptm_change[line[5]].append("Nothing")

	def analyzeSplicing(self):

		# if one of the transcripts doesn't have an isoform described, there is no use
		if not self.nIsoform or not self.tIsoform:
			return []

		allN = sorted(self.nTranscript.getSegments("isoform-specific") + self.nTranscript.getSegments("non-isoform-specific"),key=operator.itemgetter(0))
		allT = sorted(self.tTranscript.getSegments("isoform-specific") + self.tTranscript.getSegments("non-isoform-specific"),key=operator.itemgetter(0))

		correspondence = []
		c = []
		for x1,x2 in enumerate(allN):
			correspondenceFound = False
			for y1,y2 in enumerate(allT):
				intersect = list(set(x2) & set(y2))
				if intersect:
					correspondence.append([intersect,intersect])
					c.append([x1,y1])

		for thisIso,otherIso,thisIdx,otherIdx in zip([allN,allT],[allT,allN],[0,1],[1,0]):
			for pos in set(range(len(thisIso))) - set(x[thisIdx] for x in c):
				thisCorrespondence = [None,None]
				thisCorrespondence[thisIdx] = thisIso[pos]
				if pos == 0:
					if not [ x for x in c if x[otherIdx] == 0 ]:
						thisCorrespondence[otherIdx] = otherIso[pos]
				elif pos == len(thisIso) - 1:
					if not [ x for x in c if x[otherIdx] == len(otherIso) - 1 ]:
						thisCorrespondence[otherIdx] = otherIso[len(otherIso) - 1]
				else:

					prevItem = sorted([ x for x in c if x[thisIdx]==pos-1 ],key=operator.itemgetter(otherIdx),reverse=True)[0]
					nextItem = sorted([ x for x in c if x[thisIdx]==pos+1 ],key=operator.itemgetter(otherIdx))[0]
		
					if prevItem[otherIdx]+2 == nextItem[otherIdx]:
						thisCorrespondence[otherIdx] = otherIso[prevItem[otherIdx]+2]

				if thisCorrespondence not in correspondence:
					correspondence.append(thisCorrespondence)

		correspondence = sorted(correspondence,key=lambda x: x[0][0] if x[0] else x[1][0],reverse=True if self.nTranscript._strand=='-' else False )

		return correspondence

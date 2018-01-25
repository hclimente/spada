from spada import utils
from spada.biological_entities import protein, transcript

from itertools import zip_longest
import os
import operator

class IsoformSwitch:
	def __init__(self, nTx, tTx, samples):
		self._normal_transcript_name= nTx
		self._tumor_transcript_name = tTx
		self._samples 				= samples

		self._normal_transcript 	= None
		self._tumor_transcript 		= None
		self._normal_protein 		= None
		self._tumor_protein 		= None

		self._candidate 			= None
		self._noise 				= None

		#Relevance measure
		self._domain_change 		= None
		self._disorder_change 		= None
		self._ptm_change 			= None

		#Candidate measures
		self._recurrent				= None
		self._mutatedFeature		= None
		self._mutuallyExclusive		= None
		self._coocurrent			= None
		self._ppi					= None

		#Relevance measure
		self._deep_domain_change 	= { }
		self._deep_disorder_change 	= { }
		self._deep_ptm_change 		= { }

		#Neighborhood analysis
		self._neighborhood_change 	= None

	@property
	def nTx(self): return self._normal_transcript_name
	@property
	def tTx(self): return self._tumor_transcript_name
	@property
	def samples(self): return self._samples

	@property
	def isCandidate(self): return self._candidate
	def setCandidate(self, candidate): self._candidate = candidate
	@property
	def isNoise(self): return self._noise
	def setNoise(self, noise): self._noise = noise

	@property
	def nIsoform(self): return self._normal_protein
	@property
	def tIsoform(self): return self._tumor_protein
	@property
	def nTranscript(self): return self._normal_transcript
	@property
	def tTranscript(self): return self._tumor_transcript

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
		if self.disorderChange or self.domainChange or self.ptmChange:
			if self.tx_overlap:
				return True
			elif None in [self.nIsoform,self.tIsoform]:
				return True

		return False

	@property
	def recurrent(self):
		if self._recurrent is None:
			self.readCandidateAnalysis()
		return self._recurrent

	@property
	def coocurrent(self):
		if self._coocurrent is None:
			self.readCandidateAnalysis()
		return self._coocurrent

	@property
	def mutatedFeature(self):
		if self._mutatedFeature is None:
			self.readCandidateAnalysis()
		return self._mutatedFeature

	@property
	def mutuallyExclusive(self):
		if self._mutuallyExclusive is None:
			self.readCandidateAnalysis()
		return self._mutuallyExclusive

	@property
	def ppi(self):
		if self._ppi is None:
			self.readCandidateAnalysis()
		return self._ppi

	@property
	def is_candidate(self):
		if self.recurrent or self.coocurrent or self.mutatedFeature or self.mutuallyExclusive or self.ppi:
			return True

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
	def tx_overlap(self):
		nTx = set(self._normal_transcript._cds) | set(self._normal_transcript._utr)
		tTx = set(self._tumor_transcript._cds) | set(self._tumor_transcript._utr)
		return bool(nTx & tTx)

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

	def addTxInfo(self,nInfo,tInfo):
		"""Creates the transcript objects for the transcripts involved
		in the switch and calculates UTR and CDS differences."""
		self._normal_transcript = transcript.Transcript( self._normal_transcript_name, nInfo )
		self._tumor_transcript 	= transcript.Transcript( self._tumor_transcript_name, tInfo )

		self.computeCdsDiff()
		self.computeUtrDiff()

		if nInfo["proteinSequence"] and nInfo["cdsCoords"]:
			self._normal_protein = protein.Protein( self._normal_transcript_name, nInfo)

		if tInfo["proteinSequence"] and tInfo["cdsCoords"]:
			self._tumor_protein  = protein.Protein( self._tumor_transcript_name, tInfo)

		if self._normal_protein and self._tumor_protein:
			self.getAlteredRegions()

	def computeCdsDiff(self):
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

	def computeUtrDiff(self):
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
			if res.genomicPosition not in [ y.genomicPosition for y in self._tumor_protein._structure ]:
				res.setIsoformSpecific(True)

		for res in self._tumor_protein._structure:
			if res.genomicPosition not in [ y.genomicPosition for y in self._normal_protein._structure ]:
				res.setIsoformSpecific(True)

	def readRelevanceAnalysis(self):

		randomTag = "_random" if self.samples==[] else ""

		if not os.path.exists("structural_analysis/structural_summary{}.tsv".format(randomTag)):
			raise Exception("Relevance information not generated.")
			return False

		for elements in utils.readTable("structural_analysis/structural_summary{}.tsv".format(randomTag)):
			if elements[2] == self.nTx and elements[3] == self.tTx:

				if elements[5] == "True": self._domain_change = True
				elif elements[5] == "False": self._domain_change = False

				if elements[6] == "True": self._disorder_change = True
				elif elements[6] == "False": self._disorder_change = False

				if elements[8] == "True": self._ptm_change = True
				elif elements[8] == "False": self._ptm_change = False

				break

	def readDeepRelevanceAnalysis(self,skipDomain=False,skipIupred=False,skipPtm=False,filetag=""):
		if not self.is_functional:
			return

		if self.domainChange and not skipDomain:
			for line in utils.readTable("structural_analysis/interpro_analysis{}.tsv".format(filetag)):
				if line[2] == self.nTx and line[3] == self.tTx:
					self._deep_domain_change.setdefault(line[5].replace(" ","_"),[])
					if "Gained" in line[4]:
						self._deep_domain_change[line[5].replace(" ","_")].append("Gained_in_tumor")
					elif "Lost" in line[4]:
						self._deep_domain_change[line[5].replace(" ","_")].append("Lost_in_tumor")
					elif "Nothing" in line[4]:
						self._deep_domain_change[line[5].replace(" ","_")].append("Nothing")

		if self.disorderChange and not skipIupred:
			for line in utils.readTable("structural_analysis/iupred_analysis{}.tsv".format(filetag)):
				if line[2] == self.nTx and line[3] == self.tTx and float(line[-1]):
					self._deep_disorder_change.setdefault(line[5],[])
					if "Gained" in line[4]:
						self._deep_disorder_change[line[5]].append("Gained_in_tumor")
					elif "Lost" in line[4]:
						self._deep_disorder_change[line[5]].append("Lost_in_tumor")
					elif "Nothing" in line[4]:
						self._deep_disorder_change[line[5]].append("Nothing")

		if self.ptmChange and not skipPtm:
			for line in utils.readTable("structural_analysis/prosite_analysis{}.tsv".format(filetag)):
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

	def readCandidateAnalysis(self):

		if not os.path.exists("candidateList_driverEvidence.tsv"):
			raise Exception("Candidate information not generated.")
			return False

		for elements in utils.readTable("candidateList_driverEvidence.tsv"):
			if elements[3] == self.nTx and elements[4] == self.tTx:
				if elements[5] == "1": self._recurrent = True
				elif elements[5] == "0": self._recurrent = False

				if elements[6] == "1": self._mutatedFeature = True
				elif elements[6] == "0": self._mutatedFeature = False

				if elements[7] == "1": self._mutuallyExclusive = True
				elif elements[7] == "0": self._mutuallyExclusive = False

				if elements[8] == "1": self._coocurrent = True
				elif elements[8] == "0": self._coocurrent = False

				if elements[9] == "1": self._ppi = True
				elif elements[9] == "0": self._ppi = False

				break

	def analyzeDomains(self, featureType):

		featureInfo = []
		features = set()

		if not featureType in ["Pfam", "Prosite"]:
			raise Exception(featureType + ' not recognized. Use Pfam or Prosite.')

		for isoform in [self.nIsoform,self.tIsoform]:
			if isoform:
				f = isoform._prosite if featureType == "Prosite" else isoform._pfam
				[ features.add(x) for x in f ]

		for feature in features:

			featInfo = { self.nTx: [], self.tTx: []}

			for isoform in [self.nIsoform, self.tIsoform]:
				if not isoform:
					continue

				featureRegions = isoform.getFeature(featureType, feature)
				specificRegions = isoform.getSegments("isoform-specific")

				for region in featureRegions:

					thisIsosp = []
					if None not in [self.nIsoform,self.tIsoform]:
						[ thisIsosp.extend(x) for x in specificRegions if set(x) & set(region) ]
					else:
						thisIsosp = set(region)
					intersection = float(len(set(region) & set(thisIsosp)))
					domainLenght = float(len(set(region)))
					specificLength = float(len(set(thisIsosp)))

					macroScore = intersection/domainLenght
					microScore =  float("nan") if specificLength == 0 else intersection/specificLength
					jaccard = intersection/len(set(region) | set(thisIsosp))

					featInfo[isoform.tx].append({"macro": macroScore,"micro": microScore, "jaccard": jaccard})

				featInfo[isoform.tx] = sorted(featInfo[isoform.tx], key=operator.itemgetter("macro"))

			i = 1
			emptyDict = {"macro": float("nan"), "micro": float("nan"), "jaccard": float("nan")}
			for nDict,tDict in zip_longest(featInfo[self.nTx], featInfo[self.tTx], fillvalue = emptyDict):

				what = "Nothing"
				if nDict["macro"] > 0 and tDict == emptyDict:
					what = "Lost_in_tumor"
				elif tDict["macro"] > 0 and nDict == emptyDict:
					what = "Gained_in_tumor"

				f = { "feature": feature, "index": i, "what": what, \
					  "nM": nDict["macro"], "nm": nDict["micro"], "nJ": nDict["jaccard"], \
					  "tM": tDict["macro"], "tm": tDict["micro"], "tJ": tDict["jaccard"] }
				featureInfo. append(f)
				i += 1

		return featureInfo

	def analyzeIDR(self, idr_threshold):

		idrInfo = []
		features = set()

		for isoform in [self.nIsoform,self.tIsoform]:
			if isoform:
				[ features.add(x) for x in isoform._idr ]

		for feature in features:
			for protein,what in zip([self.nIsoform,self.tIsoform], ["Lost_in_tumor","Gained_in_tumor"]):
				if not protein:
					continue

				idrs = protein.getFeature("IDR", feature)
				isoform = protein.getSegments("isoform-specific")

				for idr in idrs:

					whatsHappening = what
					idrSet = set(idr)

					overlappingIsoSpecific = []
					if None not in [self.nIsoform,self.tIsoform]:
						[ overlappingIsoSpecific.extend(x) for x in isoform if set(x) & idrSet]
					else:
						overlappingIsoSpecific = idrSet

					if not overlappingIsoSpecific:
						jaccard = float("nan")
						macroScore = float("nan")
						microScore = float("nan")
						whatsHappening = "Nothing"
						significant = 0

					else:
						overlappingIsoSpecificSet = set(overlappingIsoSpecific)

						intersection = float(len(overlappingIsoSpecificSet & idrSet))
						union = float(len(overlappingIsoSpecificSet | idrSet))

						jaccard = intersection/union
						microScore = intersection/len(overlappingIsoSpecificSet)
						macroScore = intersection/len(idrSet)
						significant = int(max(microScore,macroScore) > idr_threshold)

					motifSequence = ""
					start = float("inf")
					end = float("-inf")
					for thisRes in idr:
						res = thisRes.res.upper() if thisRes.isoformSpecific else thisRes.res.lower()
						motifSequence += res
						if thisRes.num > end:
							end = thisRes.num
						if thisRes.num < start:
							start = thisRes.num

					if significant:
						f = { "feature": motifSequence, "start": start, "end": end, \
							  "what": what, "M": macroScore, "m": microScore, "J": jaccard }
						idrInfo.append(f)

		return idrInfo

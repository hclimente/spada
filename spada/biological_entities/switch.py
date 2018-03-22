from spada.io import io
from spada.biological_entities import protein, transcript

from itertools import zip_longest
import os
import operator

class IsoformSwitch:
	def __init__(self, nTx, tTx, samples):
		self._normal_transcript_name= nTx
		self._tumor_transcript_name = tTx
		self._samples 				= samples

		self._normal_transcript = None
		self._tumor_transcript 	= None
		self._normal_protein 	= None
		self._tumor_protein 	= None

		self._functional		= None
		self._main 				= None
		self._noise 			= None

		self._pfamChange		= None
		self._prositeChange 	= None
		self._idrChange			= None

	@property
	def nTx(self): return self._normal_transcript_name
	@property
	def tTx(self): return self._tumor_transcript_name
	@property
	def samples(self): return self._samples

	@property
	def isMain(self): return self._main
	def setMain(self, candidate): self._main = candidate
	@property
	def isNoise(self): return self._noise
	def setNoise(self, noise): self._noise = noise
	@property
	def isFunctional(self):

		if self._pfamChange is None:
			self.analyzeDomains('Pfam')
		if self._prositeChange is None:
			self.analyzeDomains('Prosite')
		if self._idrChange is None:
			self.analyzeIDR(0.2)

		self._functional = self._pfamChange or self._prositeChange or self._idrChange
		return self._functional

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

		if nInfo["proteinSequence"]:
			self._normal_protein = protein.Protein( self._normal_transcript_name, nInfo)

		if tInfo["proteinSequence"]:
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

					start = min([ x.num for x in region ])
					end = max([ x.num for x in region ])

					featInfo[isoform.tx].append({'macro': macroScore,'micro': microScore,
												 'jaccard': jaccard, 'start': start, 'end': end})

				featInfo[isoform.tx] = sorted(featInfo[isoform.tx], key=operator.itemgetter("macro"))

			i = 1
			emptyDict = {'macro': float('nan'), 'micro': float('nan'),
						 'jaccard': float('nan'), 'start': float('nan'),
						 'end': float('nan')}
			for nDict,tDict in zip_longest(featInfo[self.nTx], featInfo[self.tTx], fillvalue = emptyDict):

				what = "Nothing"
				if nDict["macro"] > 0 and tDict == emptyDict:
					what = "Lost_in_tumor"
				elif tDict["macro"] > 0 and nDict == emptyDict:
					what = "Gained_in_tumor"

				f = { "feature": feature, "index": i, "what": what, \
					  "nStart": nDict["start"], "nEnd": nDict["end"], \
					  "tStart": tDict["start"], "tEnd": tDict["end"], \
					  "nM": nDict["macro"], "nm": nDict["micro"], "nJ": nDict["jaccard"], \
					  "tM": tDict["macro"], "tm": tDict["micro"], "tJ": tDict["jaccard"] }
				featureInfo. append(f)
				i += 1

		if featureType == "Prosite":
			self._prositeChange = bool([ x for x in featureInfo if f["what"] != 'Nothing' ])
		else:
			self._pfamChange = bool([ x for x in featureInfo if f["what"] != 'Nothing' ])

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

		self._idrChange = bool([ x for x in idrInfo if f["what"] != 'Nothing' ])

		return idrInfo

from spada.bio import polypeptide, protein, transcript
from spada.io import io
from spada.io.error import SpadaError

from collections import namedtuple
from itertools import zip_longest
import os
import operator

LiteSwitch = namedtuple('LiteSwitch', ['ctrl', 'case', 'samples'])

class IsoformSwitch:
	def __init__(self, ctrl, case, samples):
		self._ctrl_transcript_name = ctrl
		self._case_transcript_name = case
		self._samples              = samples

		self._ctrl_transcript      = None
		self._case_transcript      = None
		self._ctrl_protein 	       = None
		self._case_protein 	       = None

		self._functional		   = None

		self._pfamChange		   = None
		self._prositeChange 	   = None
		self._idrChange			   = None

	def __repr__(self):
		functional =  '' if self.isFunctional else 'non-'
		return '{} - {} {}functional switch'.format(self.ctrl, self.case, functional)

	@property
	def ctrl(self): return self._ctrl_transcript_name
	@property
	def case(self): return self._case_transcript_name
	@property
	def samples(self): return self._samples

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
	def ctrlIsoform(self): return self._ctrl_protein
	@property
	def caseIsoform(self): return self._case_protein
	@property
	def ctrlTranscript(self): return self._ctrl_transcript
	@property
	def caseTranscript(self): return self._case_transcript

	@property
	def cds_diff(self):
		"""Returns a list with the differencen between the transcripts coding sequences."""
		cdsDiff = [ x for x in self._ctrl_transcript.cds if x not in self._case_transcript.cds ]
		cdsDiff.extend( [ x for x in self._case_transcript.cds if x not in self._ctrl_transcript.cds ])
		return cdsDiff

	def utr_diff(self, which):
		"""Returns True if there is a difference between the transcripts utr sequences."""

		if which == "5'":
			ctrl_utr = self._ctrl_transcript._5utr
			case_utr = self._case_transcript._5utr
		if which == "3'":
			ctrl_utr = self._ctrl_transcript._3utr
			case_utr = self._case_transcript._3utr

		utrDiff = [ x for x in ctrl_utr if x not in case_utr ]
		utrDiff.extend( [ x for x in case_utr if x not in ctrl_utr ])
		return utrDiff

	def addTxInfo(self,nInfo,tInfo):
		"""Creates the transcript objects for the transcripts involved
		in the switch and calculates UTR and CDS differences."""
		self._ctrl_transcript = transcript.Transcript( self._ctrl_transcript_name, nInfo )
		self._case_transcript 	= transcript.Transcript( self._case_transcript_name, tInfo )

		self.computeCdsDiff()
		self.computeUtrDiff("3'")
		self.computeUtrDiff("5'")

		if nInfo["proteinSequence"] and ( nInfo["CDS"] or not nInfo['canonical'] ):
			self._ctrl_protein = protein.Protein( self._ctrl_transcript_name, nInfo)

		if tInfo["proteinSequence"] and ( tInfo["CDS"] or not tInfo['canonical'] ):
			self._case_protein  = protein.Protein( self._case_transcript_name, tInfo)

		if self._ctrl_protein and self._case_protein:
			self.getAlteredRegions()

	def computeCdsDiff(self):
		"""Changes the values of the CDS dictionary of the transcripts to
			a bool, indicating if they are transcript specific or not."""
		for gPos in self._ctrl_transcript.cds:
			if gPos not in self._case_transcript.cds:
				self._ctrl_transcript._cds[gPos] = True
			else:
				self._ctrl_transcript._cds[gPos] = False

		for gPos in self._case_transcript.cds:
			if gPos not in self._ctrl_transcript.cds:
				self._case_transcript._cds[gPos] = True
			else:
				self._case_transcript._cds[gPos] = False

	def computeUtrDiff(self, which):
		"""Changes the values of the UTR dictionary of the transcripts to
			a bool, indicating if they are transcript specific or not."""

		if which == "5'":
			ctrl_utr = self._ctrl_transcript._5utr
			case_utr = self._case_transcript._5utr
		if which == "3'":
			ctrl_utr = self._ctrl_transcript._3utr
			case_utr = self._case_transcript._3utr

		for gPos in ctrl_utr:
			if gPos not in case_utr:
				ctrl_utr[gPos] = True
			else:
				ctrl_utr[gPos] = False

		for gPos in case_utr:
			if gPos not in ctrl_utr:
				case_utr[gPos] = True
			else:
				case_utr[gPos] = False

	def getAlteredRegions(self):
		"""Calculates the specific and non-specific residues of the isoforms
		involved in an isoform switch."""

		for res in self._ctrl_protein._structure:
			if res.genomicPosition not in [ y.genomicPosition for y in self._case_protein._structure ]:
				res.setIsoformSpecific(True)

		for res in self._case_protein._structure:
			if res.genomicPosition not in [ y.genomicPosition for y in self._ctrl_protein._structure ]:
				res.setIsoformSpecific(True)

	def analyzeSplicing(self):

		# if one of the transcripts doesn't have an isoform described, there is no use
		if not self.ctrlIsoform or not self.caseIsoform:
			return []

		allN = sorted(self.ctrlTranscript.getSegments("isoform-specific") + self.ctrlTranscript.getSegments("non-isoform-specific"),key=operator.itemgetter(0))
		allT = sorted(self.caseTranscript.getSegments("isoform-specific") + self.caseTranscript.getSegments("non-isoform-specific"),key=operator.itemgetter(0))

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

		correspondence = sorted(correspondence,key=lambda x: x[0][0] if x[0] else x[1][0],reverse=True if self.ctrlTranscript._strand=='-' else False )

		return correspondence

	def analyzeDomains(self, featureType):

		featureInfo = []
		features = set()

		if not featureType in ["Pfam", "Prosite"]:
			raise SpadaError(featureType + ' not recognized. Use Pfam or Prosite.')

		for isoform in [self.ctrlIsoform,self.caseIsoform]:
			if isoform:
				f = isoform._prosite if featureType == "Prosite" else isoform._pfam
				[ features.add(x) for x in f ]

		for feature in features:

			featInfo = { self.ctrl: [], self.case: []}

			for isoform in [self.ctrlIsoform, self.caseIsoform]:
				if not isoform:
					continue

				featureRegions = isoform.getFeature(featureType, feature)
				specificRegions = isoform.getSegments("isoform-specific")

				for region in featureRegions:

					thisIsosp = polypeptide.Polypeptide([])
					if None not in [self.ctrlIsoform,self.caseIsoform]:
						for x in specificRegions:
							if len(x & region):
								thisIsosp = thisIsosp | x
					else:
						thisIsosp = region
					intersection = float(len(region & thisIsosp))
					domainLength = float(len(region))
					specificLength = float(len(thisIsosp))

					macroScore = intersection/domainLength
					microScore =  float("nan") if specificLength == 0 else intersection/specificLength
					jaccard = intersection/len(region | thisIsosp)

					start = min([ x.num for x in region ])
					end = max([ x.num for x in region ])

					featInfo[isoform.tx].append({'macro': macroScore,'micro': microScore,
												 'jaccard': jaccard, 'start': start, 'end': end})

				featInfo[isoform.tx] = sorted(featInfo[isoform.tx], key=operator.itemgetter("macro"))

			i = 1
			emptyDict = {'macro': float('nan'), 'micro': float('nan'),
						 'jaccard': float('nan'), 'start': float('nan'),
						 'end': float('nan')}
			for nDict,tDict in zip_longest(featInfo[self.ctrl], featInfo[self.case], fillvalue = emptyDict):

				what = "Nothing"
				if nDict["macro"] > 0 and tDict == emptyDict:
					what = "Lost_in_cases"
				elif tDict["macro"] > 0 and nDict == emptyDict:
					what = "Gained_in_cases"

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

		for isoform in [self.ctrlIsoform,self.caseIsoform]:
			if isoform:
				[ features.add(x) for x in isoform._idr ]

		for feature in features:
			for protein,what in zip([self.ctrlIsoform,self.caseIsoform], ["Lost_in_cases","Gained_in_cases"]):
				if not protein:
					continue

				idrs = protein.getFeature("IDR", feature)
				specificRegions = protein.getSegments("isoform-specific")

				for idr in idrs:

					whatsHappening = what

					overlappingIsoSpecific = polypeptide.Polypeptide([])
					if None not in [self.ctrlIsoform,self.caseIsoform]:
						for x in specificRegions:
							if len(x & idr):
								overlappingIsoSpecific = overlappingIsoSpecific | x
					else:
						overlappingIsoSpecific = idr

					if not len(overlappingIsoSpecific):
						jaccard = float("nan")
						macroScore = float("nan")
						microScore = float("nan")
						whatsHappening = "Nothing"
						significant = 0

					else:

						intersection = float(len(overlappingIsoSpecific & idr))
						union = float(len(overlappingIsoSpecific | idr))

						jaccard = intersection/union
						microScore = intersection/len(overlappingIsoSpecific)
						macroScore = intersection/len(idr)
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

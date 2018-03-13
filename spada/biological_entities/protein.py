from spada.biological_entities.aminoacid import Aminoacid

class Protein:
	def __init__(self, tx, txInfo):

		self._tx				= tx
		self._gene				= txInfo["gene_id"]
		self._sequence			= txInfo["proteinSequence"]
		self._structure			= []
		self._residueCorresp	= {}
		self._pfam				= txInfo["Pfam"]
		self._prosite			= txInfo["Prosite"]
		self._idr				= txInfo["IDR"]

		if self._sequence is not None:
			cds = self.expandCDS(txInfo)
			self.mapCDStoProtein(cds, txInfo)

			self.annotateFeaturesToResidues("Pfam", self._pfam)
			self.annotateFeaturesToResidues("Prosite", self._prosite)
			self.annotateFeaturesToResidues("IDR", self._idr)

	@property
	def tx(self): return self._tx
	@property
	def seq(self): return self._sequence

	@property
	def structure(self):
		if len(self._structure) > 1:
			for res in sorted(self._structure,key=lambda x: x.num):
				yield res

	def expandExons(self, txInfo):

		plus = txInfo["strand"] == "+"
		sign = 1 if plus else -1

		for exon in txInfo["exons"]:
			exonStart = exon[not plus]
			exonEnd = exon[plus]

			for i in range(exonStart, exonEnd + sign, sign):
				yield i

	def expandCDS(self, txInfo):
		"""Assign the position of the first nucleotide of the codon to each Aminoacid."""

		cds = []

		if txInfo["CDS"]:
			cdsStart,cdsEnd = txInfo["CDS"]
			cds = [ x for x in self.expandExons(txInfo) if x >= cdsStart and x <= cdsEnd ][::3]

		return(cds)

	def mapCDStoProtein(self, cds, txInfo):

		for i in range(0, len(self._sequence)):
			aa = Aminoacid(i+1, self._sequence[i])
			self._structure.append(aa)

		if cds:
			if txInfo['start_codon'] and txInfo['stop_codon']:
				if txInfo['stop_codon']:
					if len(self._sequence) != len(cds):
						raise Exception('Transcript {}: lengths of protein sequence and CDS do not match ({} vs. {}).'.format(self._tx, len(self._structure), len(cds)))

					mrna = [ x for x in self.expandExons(txInfo) ]
					if mrna.index(cds[-1]) + 3 != mrna.index(txInfo['stop_codon']):
						raise Exception('Transcript {}: number of nucleotides in the CDS must be multiple of 3.'.format(self._tx))

				for aa,pos in zip(self._structure, cds):
					aa.setGenomicPosition(pos)

			elif txInfo['stop_codon']:
				for aa,pos in zip(reversed(self._structure), reversed(cds)):
					aa.setGenomicPosition(pos)

	def annotateFeaturesToResidues(self, featureType, features):
		for feature,region in features.items():
			for start,end in region:
				for aa in self.structure:
					if aa.num >= start and aa.num <= end:
						aa._features.add(feature)

	def getSegments(self, segmentType, minLength = 1):
		segments = []
		segment = []

		for aa in self.structure:
			flag = False
			if segmentType == "isoform-specific":
				flag = aa.isoformSpecific
			elif segmentType == "non-isoform-specific":
				flag = not aa.isoformSpecific

			if flag:
				segment.append(aa)
			elif segment:
				if len(segment) >= minLength:
					segments.append(segment)
				segment = []

		if len(segment) >= minLength:
			segments.append(segment)

		return segments

	def getFeature(self, featureType, f):

		if featureType == "Pfam": 		regions = self._pfam
		elif featureType == "Prosite": 	regions = self._prosite
		elif featureType == "IDR": 		regions = self._idr
		feature = []

		if f in regions:
			for start, end in regions[f]:
				feature.append([ x for x in self.structure ][(start - 1):end])

		return feature

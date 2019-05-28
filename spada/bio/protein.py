from spada.bio.aminoacid import Aminoacid
from spada.bio.polypeptide import Polypeptide
from spada.io.error import SpadaError, SpadaWarning

class Protein:
	def __init__(self, tx, txInfo):

		self._tx				= tx
		self._gene				= txInfo["gene_id"]
		self._sequence			= txInfo["proteinSequence"]
		self._structure			= []
		self._pfam				= txInfo["Pfam"]
		self._prosite			= txInfo["Prosite"]
		self._idr				= txInfo["IDR"]

		if self._sequence is not None:
			cds = self.expandCDS(txInfo)
			self._structure = self.mapCDStoProtein(cds, txInfo)

			self.annotateFeaturesToResidues("Pfam", self._pfam)
			self.annotateFeaturesToResidues("Prosite", self._prosite)
			self.annotateFeaturesToResidues("IDR", self._idr)

	def __repr__(self):
		return """Protein from {}.
\t- CDS: {} residues.
\t- Pfam domains: {}.
\t- ProSite features: {}.
\t+ {} IDRs.""".format(self.tx, len(self), list(self._pfam.keys()), 
					   list(self._prosite.keys()), len(self._idr.keys()) )

	def __len__(self):
		return len(self._structure)

	@property
	def tx(self): return self._tx
	@property
	def seq(self): return self._sequence

	@property
	def structure(self): return self._structure

	def expandExons(self, txInfo):

		exons = sorted(txInfo["exons"])

		for exonStart,exonEnd  in exons:
			for i in range(exonStart, exonEnd + 1):
				yield i

	def expandCDS(self, txInfo):
		"""Assign the position of the first nucleotide of the codon to each Aminoacid."""

		cds = []

		if txInfo["CDS"]:
			if txInfo['strand'] == '+':
				cdsStart = txInfo["CDS"][0] if not txInfo['start_codon'] else txInfo['start_codon']
				cdsEnd = min(txInfo["CDS"][1], txInfo['stop_codon'] - 1) if txInfo['stop_codon'] else txInfo["CDS"][1]
				cds = [ x for x in self.expandExons(txInfo) if x >= cdsStart and x <= cdsEnd ][::3]
			
			elif txInfo['strand'] == '-':
				cdsStart = txInfo["CDS"][1] if not txInfo['start_codon'] else txInfo['start_codon']
				cdsEnd = min(txInfo["CDS"][0], txInfo['stop_codon'] + 1) if txInfo['stop_codon'] else txInfo["CDS"][0]
				cds = [ x for x in self.expandExons(txInfo) if x <= cdsStart and x >= cdsEnd ][::-3]

		return(cds)

	def mapCDStoProtein(self, cds, txInfo):

		structure = []

		for i in range(0, len(self._sequence)):
			aa = Aminoacid(i+1, self._sequence[i])
			structure.append(aa)

		if cds:
			if txInfo['start_codon']:
				if txInfo['stop_codon']:
					if len(self._sequence) != len(cds):
						SpadaWarning('Transcript {}: lengths of protein sequence ({}) and CDS ({}) do not match.'.format(self._tx, len(structure), len(cds)))
						return []

					mrna = [ x for x in self.expandExons(txInfo) ]
					sign = 1 if txInfo['strand'] == '+' else -1
					if mrna.index(cds[-1]) + (sign * 3) != mrna.index(txInfo['stop_codon']):
						SpadaWarning('Transcript {}: number of nucleotides in the CDS is not a multiple of 3.'.format(self._tx))
						return []

				for aa,pos in zip(structure, cds):
					aa.setGenomicPosition(pos)

			elif txInfo['stop_codon']:
				for aa,pos in zip(reversed(structure), reversed(cds)):
					aa.setGenomicPosition(pos)

			elif len(self._sequence) == len(cds):
				for aa,pos in zip(structure, cds):
					aa.setGenomicPosition(pos)

		return structure

	def annotateFeaturesToResidues(self, featureType, features):
		for feature,region in features.items():
			for start,end in region:
				for i in range(start - 1, end):
					self.structure[i]._features.add(feature)

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
					segments.append(Polypeptide(segment))
				segment = []

		if len(segment) >= minLength:
			segments.append(Polypeptide(segment))

		return segments

	def getFeature(self, featureType, f):

		if featureType == "Pfam": 		regions = self._pfam
		elif featureType == "Prosite": 	regions = self._prosite
		elif featureType == "IDR": 		regions = self._idr
		segments = []

		if f in regions:
			for start, end in regions[f]:
				segments.append(Polypeptide(self.structure[(start - 1):end]))

		return segments

from spada.biological_entities import aminoacid
from spada import utils

import logging
import os

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
			for res in range(0, len(self._sequence) ):
				self._structure.append( aminoacid.AminoAcid(res+1, self._sequence[res]) )

			self.mapResiduesToGenome(txInfo["exonStructure"], txInfo["cdsCoords"], txInfo["strand"])
			self.annotateFeaturesToResidues("Pfam", self._pfam)
			self.annotateFeaturesToResidues("Prosite", self._prosite)
			self.annotateFeaturesToResidues("IDR", self._idr)

	@property
	def tx(self): return self._tx
	@property
	def seq(self): return self._sequence

	@property
	def structure_ordered(self):
		if len(self._structure) > 1:
			for res in sorted(self._structure,key=lambda x: x.num):
				yield res

	def mapResiduesToGenome(self, exons, cds, strand):
		"""Assign the position of the first nucleotide of the codon to each AminoAcid."""

		genomicPositions = []
		gap  = 0

		#Generate a list with the genomic position of all codons of the CDS.
		if strand == "+":
			cdsStart = cds[0]
			cdsEnd	 = cds[1]

			for exonStart,exonEnd in exons:
				if exonEnd < cdsStart or exonStart > cdsEnd: continue

				start = exonStart if exonStart >= cdsStart else cdsStart
				if gap != 0: start += 3 - gap
				end = exonEnd if exonEnd <= cdsEnd else cdsEnd

				for gPos in range(start, end, 3):
					genomicPositions.append(gPos)

				gap = (end - start)%3

		elif strand == "-":
			cdsStart = cds[1]
			cdsEnd	 = cds[0]

			for exonEnd,exonStart in [ (x-1,y-1) for x,y in reversed(exons) ]:
				if exonEnd > cdsStart or exonStart < cdsEnd: continue

				start = exonStart if exonStart <= cdsStart else cdsStart
				if gap != 0: start -= 3 - gap
				end = exonEnd if exonEnd >= cdsEnd else cdsEnd

				for gPos in range(start, end, -3):
					genomicPositions.append(gPos)

				gap = (start - end)%3

		if len(self._structure) != len(genomicPositions) - 1:
			raise Exception('Transcript {}: lengths of protein sequence and CDS do not match ({} vs. {}).'.format(self._tx, len(self._structure), len(genomicPositions) - 1))

		for aa, gPos in zip(self._structure, genomicPositions[:-1]): #Remove the stop codon
			aa.setGenomicPosition(gPos)

	def annotateFeaturesToResidues(self, featureType, features):
		for feature,region in features.items():
			for start,end in region:
				for aa in self.structure_ordered:
					if aa.num >= start and aa.num <= end:
						aa._features.add(feature)

	def getSegments(self, segmentType, minLength = 1):
		segments = []
		segment = []

		for aa in self.structure_ordered:
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
				feature.append([ x for x in self.structure_ordered ][(start - 1):end])

		return feature

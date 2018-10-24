from spada.io import io

from collections import OrderedDict
import logging
import os

class Transcript:
	def __init__(self, tx, txInfo):

		self._name  			= tx
		self._exons 			= txInfo["exons"]
		self._cds_coordinates 	= txInfo["CDS"]
		self._tx_coordinates 	= txInfo["txCoords"]
		self._strand 			= txInfo["strand"]

		# dictionaries clasifying every the nucleotide as either CDS or UTR
		# key: genomic positon; value: is it isoform specific in the switch?
		self._cds = OrderedDict()
		self._5utr = OrderedDict()
		self._3utr = OrderedDict()
		self._exon = OrderedDict()

		if self._strand == "+":
			if self._cds_coordinates is None:
				cdsStart = float("Inf")
				cdsEnd 	 = float("-Inf")
			else:
				cdsStart = self._cds_coordinates[0]
				cdsEnd 	 = self._cds_coordinates[1]

			exon = 1
			for exonStart,exonEnd in self._exons:
				for gPos in range(exonStart, exonEnd + 1):
					if self._cds_coordinates:
						if gPos >= cdsStart and gPos <= cdsEnd:
							self._cds.setdefault(gPos, None)
						elif gPos <= cdsStart:
							self._5utr.setdefault(gPos, None)
						elif gPos >= cdsEnd:
							self._3utr.setdefault(gPos, None)
					else:
						self._5utr.setdefault(gPos, None)

					self._exon[gPos] = exon
				exon += 1

		elif self._strand == "-":
			if self._cds_coordinates is None:
				cdsStart = float("-Inf")
				cdsEnd 	 = float("Inf")
			else:
				cdsStart = self._cds_coordinates[1]
				cdsEnd 	 = self._cds_coordinates[0]

			# iterate in reverse genomic order, still 5'->3' in the - strand
			exon = 1
			for exonEnd,exonStart in reversed(sorted(self._exons)):
				for gPos in range(exonStart, exonEnd - 1, -1):
					if self._cds_coordinates:
						if gPos <= cdsStart and gPos >= cdsEnd:
							self._cds.setdefault(gPos, None)
						elif gPos >= cdsStart:
							self._5utr.setdefault(gPos, None)
						elif gPos <= cdsEnd:
							self._3utr.setdefault(gPos, None)
					else:
						self._5utr.setdefault(gPos, None)
				self._exon[gPos] = exon
				exon += 1

	@property
	def name(self): return self._name
	@property
	def utr5(self): return self._5utr
	@property
	def utr3(self): return self._3utr
	@property
	def cds(self):
		if len(self._cds) == 1:
			return {}
		else:
			return self._cds

	def getSegments(self,thing,minLength=1,gap=0):
		segments = []
		segment = []
		gapped = []

		reverseOrder = False if self._strand=='+' else True

		for nt in sorted(self._cds,reverse=reverseOrder):
			flag = False
			if thing =="isoform-specific":
				flag = self._cds[nt]
			elif thing =="non-isoform-specific":
				flag = not self._cds[nt]

			if flag:
				if gapped:
					segment.extend(gapped)
					gapped = []
				segment.append(nt)
			elif segment:
				if len(gapped) < gap:
					gapped.append(nt)
				else:
					if len(segment) >= minLength:
						segments.append(segment)

					gapped = []
					segment = []

		if len(segment) >= minLength:
			segments.append(segment)

		return segments

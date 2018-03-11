from spada.io import io

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
		self._cds = {}
		self._utr = {}
		self._exon = {}

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
						else:
							self._utr.setdefault(gPos, None)
					else:
						self._utr.setdefault(gPos, None)

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
			for exonEnd,exonStart in reversed(self._exons):
				for gPos in range(exonStart, exonEnd - 1, -1):
					if self._cds_coordinates:
						if gPos <= cdsStart and gPos >= cdsEnd:
							self._cds.setdefault(gPos, None)
						else:
							self._utr.setdefault(gPos, None)
					else:
						self._utr.setdefault(gPos, None)
				self._exon[gPos] = exon
				exon += 1

	@property
	def name(self): return self._name
	@property
	def utr(self): return self._utr
	@property
	def cds(self):
		if len(self._cds) == 1:
			return {}
		else:
			return self._cds

	@property
	def cds_ordered(self):
		if len(self._cds) > 1:
			reverseOrder = False if self._strand=='+' else True

			for nt in sorted(self._cds,reverse=reverseOrder):
				yield nt

	@property
	def cds_exclusive(self):
		exclusive = float(sum([ 1 for x in self._cds if self._cds[x] ]))
		nonExclusive = float(sum([ 1 for x in self._cds if not self._cds[x] ]))
		try:
			return exclusive/(exclusive+nonExclusive)
		except ZeroDivisionError:
			return None

	@property
	def tx_exclusive(self):
		tx = self._cds.copy()
		tx.update(self._utr)

		exclusive = float(sum([ 1 for x in tx if tx[x] ]))
		nonExclusive = float(sum([ 1 for x in tx if not tx[x] ]))
		try:
			return exclusive/(exclusive+nonExclusive)
		except ZeroDivisionError:
			return None

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

import logging

class Transcript:
	def __init__(self, name, properties):
		self._name  			= name
		self._exons 			= properties["exonStructure"]
		self._cds_coordinates 	= properties["cdsCoords"]
		self._tx_coordinates 	= properties["txCoords"]
		self._strand 			= properties["strand"]

		#Create two lists, containing all the nucleotides from the transcript, classified in CDS and UTR.
		self._cds 		= []
		self._utr 		= []
		
		if strand == "+": 
			cdsStart = self._cds_coordinates[0]
			cdsEnd 	 = self._cds_coordinates[1]
			
			for exonStart,exonEnd in self._exons:
				for gPos in range(exonStart, exonEnd):
					if self._cds_coordinates:
						if gPos >= cdsStart and gPos <= cdsEnd:
							self._cds.append(gPos)
						else:
							self._utr.append(gPos)
					else:
						self._utr.append(gPos)

		elif strand == "-": 
			cdsStart = self._cds_coordinates[1]
			cdsEnd 	 = self._cds_coordinates[0]

			#In strand -, iterate in reverse genomic order, still 5'->3'.
			#Subtract 1 from the exon coordinates to convert the UCSC format to the reverse direction.
			for exonEnd,exonStart in [ (x-1,y-1) for x,y in reversed(self._exons) ]:
				for gPos in range(exonStart, exonEnd,-1):
					if self._cds_coordinates:
						if gPos <= cdsStart and gPos >= cdsEnd:
							self._cds.append(gPos)
						else:
							self._utr.append(gPos)
					else:
						self._utr.append(gPos)

	@property
	def name(self): return self._name
	@property
	def utr(self): return self._utr
	@property
	def cds(self): return self._cds

	def get_cdsDiff(self, otherTranscript):
		"""Makes a list with the genomic position not shared between the CDS of two transcripts."""
		cdsDiffs = []
		for gPos in self._cds:
			if gPos not in otherTranscript.cds:
				cdsDiffs.append(gPos)

		for gPos in otherTranscript.cds:
			if gPos not in self._cds:
				cdsDiffs.append(gPos)
		
		return cdsDiffs

	def get_utrDiff(self, otherTranscript):
		"""Makes a list with the genomic position not shared between the UTR of two transcripts."""
		utrDiffs = []
		for gPos in self._utr:
			if gPos not in otherTranscript.utr:
				utrDiffs.append(gPos)

		for gPos in otherTranscript.utr:
			if gPos not in self._utr:
				utrDiffs.append(gPos)

		return utrDiffs

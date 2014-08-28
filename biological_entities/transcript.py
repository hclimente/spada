import logging

class Transcript:
	def __init__(self, name, properties):
		self._name  			= name
		self._exons 			= properties["exonStructure"]
		self._cds_coordinates 	= properties["cdsCoords"]
		self._tx_coordinates 	= properties["txCoords"]
		self._strand 			= properties["strand"]

		#Create two dicts, clasifying all the nucleotides in CDS and UTR
		#Value: isoform specific in this switch
		self._cds = {}
		self._utr = {}
		
		if self._strand == "+": 
			cdsStart = self._cds_coordinates[0]
			cdsEnd 	 = self._cds_coordinates[1]
			
			for exonStart,exonEnd in self._exons:
				for gPos in range(exonStart, exonEnd):
					if self._cds_coordinates:
						if gPos >= cdsStart and gPos <= cdsEnd:
							self._cds.setdefault(gPos, None)
						else:
							self._utr.setdefault(gPos, None)
					else:
						self._utr.setdefault(gPos, None)

		elif self._strand == "-": 
			cdsStart = self._cds_coordinates[1]
			cdsEnd 	 = self._cds_coordinates[0]

			#In strand -, iterate in reverse genomic order, still 5'->3'.
			#Subtract 1 from the exon coordinates to convert the UCSC format to the reverse direction.
			for exonEnd,exonStart in [ (x-1,y-1) for x,y in reversed(self._exons) ]:
				for gPos in range(exonStart, exonEnd,-1):
					if self._cds_coordinates:
						if gPos <= cdsStart and gPos >= cdsEnd:
							self._cds.setdefault(gPos, None)
						else:
							self._utr.setdefault(gPos, None)
					else:
						self._utr.setdefault(gPos, None)

	@property
	def name(self): return self._name
	@property
	def utr(self): return self._utr
	@property
	def cds(self): return self._cds
	@property
	def cds_exclusive(self): 
		exclusive = float(sum([ 1 for x in self._cds if self._cds[x] ]))
		nonExclusive = float(sum([ 1 for x in self._cds if not self._cds[x] ]))
		try:
			return exclusive/(exclusive+nonExclusive)
		except ZeroDivisionError:
			return None
	
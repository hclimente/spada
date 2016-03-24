from libs import utils
from libs import options

import logging
import os

class Transcript:
	def __init__(self, name, properties):

		self._name  			= name
		self._exons 			= properties["exonStructure"]
		self._cds_coordinates 	= properties["cdsCoords"]
		self._tx_coordinates 	= properties["txCoords"]
		self._strand 			= properties["strand"]

		self._pfam 				= {}
		self._ptms 				= {}

		#Create two dicts, clasifying all the nucleotides in CDS and UTR
		#Value: isoform specific in this switch
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
				for gPos in range(exonStart, exonEnd):
					if self._cds_coordinates:
						# notation in ucsc format: first pos included and last excluded
						if gPos >= cdsStart and gPos < cdsEnd:
							self._cds.setdefault(gPos, None)
						else:
							self._utr.setdefault(gPos, None)
					else:
						self._utr.setdefault(gPos, None)

					self._exon[gPos]=exon
				exon +=1

		elif self._strand == "-": 
			if self._cds_coordinates is None:
				cdsStart = float("-Inf")
				cdsEnd 	 = float("Inf")
			else:
				cdsStart = self._cds_coordinates[1]
				cdsEnd 	 = self._cds_coordinates[0]

			#In strand -, iterate in reverse genomic order, still 5'->3'.
			#Subtract 1 from the exon coordinates to convert the UCSC format to the reverse direction.
			exon = 1
			for exonEnd,exonStart in [ (x-1,y-1) for x,y in reversed(self._exons) ]:
				for gPos in range(exonStart, exonEnd,-1):
					if self._cds_coordinates:
						if gPos < cdsStart and gPos >= cdsEnd:
							self._cds.setdefault(gPos, None)
						else:
							self._utr.setdefault(gPos, None)
					else:
						self._utr.setdefault(gPos, None)
				self._exon[gPos]=exon
				exon +=1

	@property
	def name(self): return self._name
	@property
	def utr(self): return self._utr
	@property
	def cds(self): 
		if len(self._cds)==1:
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

	def readPfamDomains(self):
		featFile = "{}data/{}/InterPro/{}.tsv".format(options.Options().wd,options.Options().annotation,self.name)

		if not os.path.exists(featFile):
			return
		elif self._pfam:
			# already read
			return

		for line in utils.readTable(featFile,header=False):

			if line[3] != "Pfam":
				continue

			domainId = "{0}|{1}".format(line[4],line[5]).replace(" ","_")

			self._pfam.setdefault(domainId,[])
			self._pfam[domainId].append((int(line[6]),int(line[7])))

	def readProsite(self):

		featFile = "{}data/{}/ProSite/{}.out".format(options.Options().wd,options.Options().annotation,self.name)

		if not os.path.exists(featFile) or os.stat(featFile).st_size == 0:
			return
		elif self._ptms:
			# already read
			return

		for line in utils.readTable(featFile,header=False):

			prositeId = line[-1].replace(" ","_")

			self._ptms.setdefault(prositeId,[])
			self._ptms[prositeId].append((int(line[-3].replace(" -","")),int(line[-2])))


	def readIupred(self):

		featFile = "{}data/{}/ProSite/{}.out".format(options.Options().wd,options.Options().annotation,self.name)

		if not os.path.exists(featFile) or os.stat(featFile).st_size == 0:
			return
		elif self._ptms:
			# already read
			return

		for line in utils.readTable(featFile,header=False):

			prositeId = line[-1].replace(" ","_")

			self._ptms.setdefault(prositeId,[])
			self._ptms[prositeId].append((int(line[-3].replace(" -","")),int(line[-2])))

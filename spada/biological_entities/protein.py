from spada.biological_entities import aminoacid
from spada.interface import interpro_analysis
from spada import utils

from Bio import pairwise2
import logging
import os

class Protein:
	def __init__(self, tx, txInfo):

		self._tx				= tx
		self._gene				= txInfo["gene_id"]
		self._uniprot			= txInfo["Uniprot"]
		self._sequence			= txInfo["proteinSequence"]
		self._structure			= []
		self._residueCorresp	= {}
		self._has_pdbs			= False
		self._pfam				= []
		self._prosite			= {}
		self._pdbs				= {}

		if self._sequence is not None:
			for res in range(0, len(self._sequence) ):
				self._structure.append( aminoacid.AminoAcid(res+1, self._sequence[res]) )

			self.mapResiduesToGenome(txInfo["exonStructure"], txInfo["cdsCoords"], txInfo["strand"])

	@property
	def tx(self): return self._tx
	@property
	def uniprot(self): return self._uniprot
	@property
	def hasPdbs(self): return self._has_pdbs
	@property
	def seq(self): return self._sequence

	@property
	def structure_ordered(self):
		if len(self._structure) > 1:
			for res in sorted(self._structure,key=lambda x: x.num):
				yield res

	def mapSubsequence(self, querySequence, strict=True):

		#Get alignment
		alignments = pairwise2.align.globalms(self._sequence, querySequence, 2, -1, -.5, -.1)

		oriSeq = alignments[0][0]
		subSeq = alignments[0][1]

		posOri = 0
		posSub = 0
		correspondence = []

		for lPos in range(0, len(subSeq) ):

			if subSeq[lPos] != "-" and posOri < len(self._sequence):
				thisRes = querySequence[posSub]

				if not strict:
					correspondence.append(posOri)
				elif self._sequence[posOri] == thisRes:
					correspondence.append(posOri)
				else:
					logging.debug("Unmatched residue {0} != {1}.".format(self._sequence[posOri], thisRes))

				posSub += 1

			if oriSeq[lPos] != "-":
				posOri += 1

		return correspondence

	def calculateVolumes(self,interactionPdb,chain):

		self._has_pdbs = True

		#Generate the pdb object, extract the chain of interest and calculate volumes

		pdb = PDB(interactionPdb)
		chainObj = pdb.get_chain_by_id(chain)
		try:
			chainObj.calculate_dssp()
		except AttributeError:
			return False

		interacting_surface		= []
		non_interacting_surface = []
		buried					= []

		#Get interface
		ppComplex = Complex(pdb, PNI=False, PHI=False)

		for i in ppComplex.PPInterfaces:
			if i.protein_chain.chain == chain:
				interacting_surface.extend([ x.identifier for x in i.protein_view_interface.keys()])
			elif i.protein_interactor.chain == chain:
				interacting_surface.extend([ x.identifier for x in i.interactor_view_interface.keys()])

		# Get non-interface (surface and core)
		non_interacting_surface.extend( [ x.identifier for x in chainObj.aminoacids if x.identifier not in interacting_surface and x.exposed] )
		buried.extend( [ x.identifier for x in chainObj.aminoacids if x.identifier not in interacting_surface and not x.exposed] )

		#Map the sequence found in the PDB to the full sequence
		a = self.mapSubsequence(chainObj.protein_sequence)

		for x in range(0, len(a)):
			thisRes = chainObj.aminoacids[x]

			if thisRes.identifier in buried:
				self._structure[a[x]].setTag("B")
				self._structure[a[x]]._pdbMapping[interactionPdb] = ( chainObj.chain, thisRes.identifier, "B" )
				logging.debug("{0}: residue {1}-{2} ({3} in sequence) detected as buried.".format(
											interactionPdb, chainObj.chain, thisRes.identifier, x))
			elif thisRes.identifier in non_interacting_surface:
				self._structure[a[x]].setTag("NIS")
				self._structure[a[x]]._pdbMapping[interactionPdb] = ( chainObj.chain, thisRes.identifier, "NIS" )
				logging.debug("{0}: residue {1}-{2} ({3} in sequence) detected as non-interacting surface.".format(
											interactionPdb, chainObj.chain, thisRes.identifier, x))
			elif thisRes.identifier in interacting_surface:
				self._structure[a[x]].setTag("IS")
				self._structure[a[x]]._pdbMapping[interactionPdb] = ( chainObj.chain, thisRes.identifier, "IS" )
				logging.debug("{0}: residue {1}-{2} ({3} in sequence) detected as interacting surface.".format(
											interactionPdb, chainObj.chain, thisRes.identifier, x))

		return True

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

		for aminoAcid, gPos in zip(self._structure, genomicPositions[:-1]): #Remove the stop codon
			aminoAcid.setGenomicPosition(gPos)

	def getFeatures(self,interproOut):
		for featInfo in interpro_analysis.InterproAnalysis().readInterpro(interproOut,self):
			self._pfam.append(featInfo)

	def getSegments(self,thing,minLength=1,gap=0):
		segments = []
		segment = []
		gapped = []

		for res in self.structure_ordered:
			flag = False
			if thing =="isoform-specific":
				flag = res.isoformSpecific
			elif thing =="non-isoform-specific":
				flag = not res.isoformSpecific
			elif thing =="anchor":
				flag = res.isAnchored
			elif thing =="disordered":
				flag = res.isDisordered
			else:
				flag = thing in res._feature

			if flag:
				if gapped:
					segment.extend(gapped)
					gapped = []
				segment.append(res)
			elif segment:
				if len(gapped) < gap:
					gapped.append(res)
				else:
					if len(segment) >= minLength:
						segments.append(segment)

					gapped = []
					segment = []

		if len(segment) >= minLength:
			segments.append(segment)

		return segments

	def readAnchor(self):
		pass
		"""
		outfile = "{0}data/{1}/ANCHOR/{2}.txt".format(options.Options().wd,options.Options().annotation,self.tx)

		if not os.path.isfile(outfile):
			out = []
			fasta = "{0}{1}{2}.fa".format(options.Options().qout,self.tx,options.Options().filetag)
			with open(fasta,"w") as FASTA:
				FASTA.write(">{0}\n{1}\n".format(self.tx,self.seq))

			proc = utils.cmdOut("anchor",fasta)
			out = [ x.strip().split() for x in proc.stdout if "#" not in x ]
			os.remove(fasta)

			with open(outfile,"w") as ANCHORout:
				for line in out:
					resNum  = int(line[0])
					residue	= line[1]
					score	= float(line[2].strip())
					ANCHORout.write("{0}\t{1}\t{2}\n".format(resNum,residue,score))
					thisRes = self._structure[resNum-1]

					if residue != thisRes.res:
						self.logger.error("Not matching residue in ANCHOR analysis, transcript {0}.".format(self.tx))
						continue

					thisRes.set_anchorScore(score)

		else:
			for line in utils.readTable(outfile,header=False):
				resNum  = int(line[0])
				residue	= line[1]
				score	= float(line[2].strip())
				thisRes = self._structure[resNum-1]

				if residue != thisRes.res:
					self.logger.error("Not matching residue in ANCHOR analysis, transcript {0}.".format(self.tx))
					continue

				thisRes.set_anchorScore(score)
		"""

	def readIupred(self,mode):
		pass
		"""
		outfile = "{}data/{}/IUPred/{}.{}.txt".format(options.Options().wd,options.Options().annotation,self.tx,mode)

		if not os.path.isfile(outfile):
			out = []
			fasta = "{}{}{}.fa".format(options.Options().qout,self.tx,options.Options().filetag)
			with open(fasta,"w") as FASTA:
				FASTA.write(">{0}\n{1}\n".format(self.tx,self.seq))

			proc = utils.cmdOut(options.Options().wd+"pipeline/libs/bin/iupred/iupred",fasta,mode)
			out = [ x.strip().split(" ") for x in proc.stdout if "#" not in x ]
			os.remove(fasta)

			with open(outfile,"w") as IUout:
				for line in out:
					resNum  = int(line[0])
					residue	= line[1]
					score	= float(line[-1])
					IUout.write("{0}\t{1}\t{2}\n".format(resNum,residue,score))

					thisRes = self._structure[resNum-1]
					if residue != thisRes.res:
						self.logger.error("Not matching residue in ANCHOR analysis, transcript {0}.".format(self.tx))
						continue
					thisRes.set_iuPredScore(score)

		else:
			for line in utils.readTable(outfile,header=False):
				resNum  = int(line[0])
				residue	= line[1]
				score	= float(line[-1])
				thisRes = self._structure[resNum-1]

				if residue != thisRes.res:
					self.logger.error("Not matching residue in ANCHOR analysis, transcript {0}.".format(self.tx))
					continue

				thisRes.set_iuPredScore(score)
		"""

	def readInterpro(self):
		pass
		"""
		outfile = "{0}data/{1}/InterPro/{2}.tsv".format(options.Options().wd,options.Options().annotation,self.tx)
		acceptedAnalysis = ["Pfam"]

		if not os.path.isfile(outfile):
			return

		else:
			for line in utils.readTable(outfile,header=False):
				# https://code.google.com/p/interproscan/wiki/OutputFormats#Tab-separated_values_format_%28TSV%29
				protein_accession	= line[0] #Protein Accession
				md5_digest			= line[1] #Sequence MD5 digest
				seq_length			= line[2] #Sequence Length
				analysis			= line[3] #Analysis
				signature_accession = line[4] #Signature Accession
				signature_descript  = line[5].replace(" ","_") #Signature Description
				start				= int(line[6]) #Start location
				stop				= int(line[7]) #Stop location
				try:
					 score			= float(line[8]) #Score - is the e-value of the match reported by member database method
				except ValueError:
					 score			= None
				status				= True if line[9] == "T" else False #Status - is the status of the match (T: true)
				date				= line[10] #Date - is the date of the run
				if len(line) > 11:
					interpro_annotation = line[11] #(InterPro annotations - accession)
				else:
					interpro_annotation = ""
				if len(line) > 12:
					interpro_descript	= line[12] #(InterPro annotations - description)
				else:
					interpro_descript	= ""
				if len(line) > 13:
					go_annotation		= line[13] #(GO annotations)
				else:
					go_annotation		= ""
				if len(line) > 14:
					pathway_annotation  = line[14] #(Pathways annotations)
				else:
					pathway_annotation  = ""
				tx = protein_accession.strip().split("#")[0]
				if score and score > 0.01:
					 continue
				elif analysis not in acceptedAnalysis:
					 continue
				elif tx != self.tx:
					 raise Exception("Error reading InterPro features for {0}. Invalid identifier {1} found.".format(self.tx,tx))

				for resNum in range(start,stop):
					thisRes = self._structure[resNum-1]
					thisRes._feature.append("{0}|{1}".format(signature_accession,signature_descript))

				featInfo = {}
				featInfo["region"]		= [start,stop]
				featInfo["accession"]	= signature_accession
				featInfo["description"]	= signature_descript
				featInfo["analysis"]	= analysis
				featInfo["go"]	 		= go_annotation
				self._pfam.append(featInfo)
		"""

	def readProsite(self):
		pass
		"""
		featFile = "{0}data/{1}/ProSite/{2}.out".format(options.Options().wd,options.Options().annotation,self.tx)

		if not os.path.exists(featFile) or os.stat(featFile).st_size == 0:
			return

		for line in utils.readTable(featFile,header=False):

			prositeId = line[-1].replace(" ","_")

			self._prosite.setdefault(prositeId,[])
			self._prosite[prositeId].append((int(line[-3].replace(" -","")),int(line[-2])))

			for ptm in self._prosite:
				for start,end in self._prosite[ptm]:
					for resNum in range(start,end):
						thisRes = self._structure[resNum-1]
						thisRes._feature.append(ptm)
		"""

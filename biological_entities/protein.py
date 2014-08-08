from biological_entities import aminoacid
from libs.SBI.structure import PDB
from libs.SBI.structure.contacts import Complex

from Bio import pairwise2
import logging

class Protein:
	def __init__(self, tx, uniprot, seq, exons, cds, strand):
		self.logger 			= logging.getLogger(__name__)
		self._tx 				= tx
		self._uniprot 			= uniprot
		self._sequence 			= seq
		self._structure			= []
		self._residueCorresp	= {}
		self._has_pdbs			= False

		if self._sequence is not None:
			for res in range(0, len(self._sequence) ):
				self._structure.append( aminoacid.AminoAcid(res+1, self._sequence[res]) )

			self.mapResiduesToGenome(exons, cds, strand)

	@property
	def tx(self): return self._tx
	@property
	def uniprot(self): return self._uniprot
	@property
	def hasPdbs(self): return self._has_pdbs	

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
					self.logger.debug("Unmatched residue {0} != {1}.".format(self._sequence[posOri], thisRes))

				posSub += 1

			if oriSeq[lPos] != "-": 
				posOri += 1

		return correspondence

	def calculateVolumes(self, interaction, chain, exons):

		self._has_pdbs = True

		#Generate the pdb object, extract the chain of interest and calculate volumes

		pdb = PDB(interaction)
		chainObj = pdb.get_chain_by_id(chain)
		try:
			chainObj.calculate_dssp()
		except AttributeError:
			return False
		
		interacting_surface 	= []
		non_interacting_surface = []
		buried 					= []

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
				self._structure[a[x]]._pdbMapping[interaction] = ( chainObj.chain, thisRes.identifier, "B" )
				self.logger.debug("{0}: residue {1}-{2} ({3} in sequence) detected as buried.".format(
											interaction, chainObj.chain, thisRes.identifier, x))
			elif thisRes.identifier in non_interacting_surface:
				self._structure[a[x]].setTag("NIS")
				self._structure[a[x]]._pdbMapping[interaction] = ( chainObj.chain, thisRes.identifier, "NIS" )
				self.logger.debug("{0}: residue {1}-{2} ({3} in sequence) detected as non-interacting surface.".format(
											interaction, chainObj.chain, thisRes.identifier, x))
			elif thisRes.identifier in interacting_surface:
				self._structure[a[x]].setTag("IS")
				self._structure[a[x]]._pdbMapping[interaction] = ( chainObj.chain, thisRes.identifier, "IS" )
				self.logger.debug("{0}: residue {1}-{2} ({3} in sequence) detected as interacting surface.".format(
											interaction, chainObj.chain, thisRes.identifier, x))

		return True

	def mapResiduesToGenome(self, exons, cds, strand):
		"""Assign the position of the first nucleotide of the codon to each AminoAcid."""

		genomicPositions = []
		gap  = 0

		#Generate a list with the genomic position of all codons of the CDS.
		if strand == "+": 
			cdsStart = cds[0]
			cdsEnd 	 = cds[1]
			
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
			cdsEnd 	 = cds[0]

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

	def getAlteredRegions(self, otherIsoform):
		for thisResidue in self._structure:
			if thisResidue.genomicPosition not in [ y.genomicPosition for y in otherIsoform._structure]:
				thisResidue.setIsoformSpecific(True)

	def report(self):
		seq 	= ""
		isoSp 	= ""
		for aResidue in self._structure:
			if not aResidue.tag: 		seq += "*"
			elif aResidue.tag == "NIS":	seq += "S"
			elif aResidue.tag == "IS": 	seq += "I"
			elif aResidue.tag == "B": 	seq += "B"

			if aResidue.isoformSpecific: 	isoSp += "X"
			else: 							isoSp += "-"

		self.logger.debug("{0},{1}".format(self._tx, seq))
		self.logger.debug("{0},{1}".format(self._tx, isoSp))

	def printPDBInfo(self):
		for pdb in set([ y for x in self._structure for y in x._pdbMapping ]):
			
			chainsSet = set([ x._pdbMapping[pdb][0] for x in self._structure if pdb in x._pdbMapping ])
			chain = set(chainsSet).pop()
			if len(chainsSet) > 1:
				self.logger.error("More than one chain for a protein in the PDB {0}, chains {1}.".format(pdb, chainsSet))
			isoSpec = [ x._pdbMapping[pdb][1][:-1] for x in self._structure if pdb in x._pdbMapping and x.isoformSpecific ]
			interact = [ x._pdbMapping[pdb][1][:-1] for x in self._structure if pdb in x._pdbMapping and x._pdbMapping[pdb][2]=="IS" ]

			if not isoSpec or not interact:
				return

			self.logger.debug("{0}, load {1}".format(self._tx,pdb))
			self.logger.debug("{0}, set bg_rgb=[1,1,1]".format(self._tx))
			self.logger.debug("{0}, hide all".format(self._tx))
			self.logger.debug("{0}, show cartoon, all".format(self._tx))
			self.logger.debug("{0}, color black, all".format(self._tx))
			self.logger.debug("{0}, color white, chain {1}".format(self._tx,chain))
			self.logger.debug("{0}, select interact, chain {1} & resi {2}".format(self._tx,chain, "+".join(interact) ))
			self.logger.debug("{0}, select isoSpecific, chain {1} & resi {2}".format(self._tx,chain, "+".join(isoSpec) ))
			self.logger.debug("{0}, color orange, isoSpecific".format(self._tx))
			self.logger.debug("{0}, color red, interact & isoSpecific".format(self._tx))
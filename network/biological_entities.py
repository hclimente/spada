from libs import utils
from libs.SBI.structure import PDB
from libs.SBI.structure.contacts import Complex

from collections import Counter
from Bio import pairwise2
import logging

import pdb

class AminoAcid:
	def __init__(self, resNum, resName):
		self.logger = logging.getLogger(__name__)
		self._num 				= resNum
		self._res 				= resName
		self._tag 				= []
		self._isoformSpecific 	= False
		self._pdbMapping 		= {}
		self._genomicPosition 	= None

	@property
	def res(self): return self._res
	@property
	def genomicPosition(self): return self._genomicPosition
	@property
	def tag(self): return self._tag
	def setTag(self,tag): self._tag.append(tag)

	def setIsoformSpecific(self,isoformSpecific): self._isoformSpecific=isoformSpecific

class Protein:
	def __init__(self, tx, uniprot, seq, exons):
		self.logger = logging.getLogger(__name__)
		self._tx 				= tx
		self._uniprot 			= uniprot
		self._sequence 			= seq
		self._structure			= []
		self._residueCorresp	= {}
		self._has_pdbs			= False

		if self._sequence is not None:
			for res in range(0, len(self._sequence) ):
				self._structure.append( AminoAcid(res+1, self._sequence[res]) )

			self.mapResiduesToGenome(exons)

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
			if subSeq[lPos] != "-": 
				thisRes = querySequence[posSub]

				if not strict:
					correspondence.append(posOri)
				elif self._sequence[posOri] == thisRes:
					correspondence.append(posOri)
				else:
					self.logger.debug("Unmatched residue {0} != {1}".format(self._sequence[posOri], thisRes))

				posSub += 1

			if oriSeq[lPos] != "-": posOri += 1

		return correspondence


	def calculateVolumes(self, interaction, chain, exons):

		self._has_pdbs = True

		#Generate the pdb object, extract the chain of interest and calculate volumes
		pdb = PDB(interaction)
		chainObj = pdb.get_chain_by_id(chain)
		chainObj.calculate_dssp()
		
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
				self._structure[a[x]]._pdbMapping[interaction] = ( chainObj.chain, thisRes.identifier )
				self.logger.debug("{0}: residue {1}-{2} ({3} in sequence) detected as buried.".format(
											interaction, chainObj.chain, thisRes.identifier))
			elif thisRes.identifier in non_interacting_surface:
				self._structure[a[x]].setTag("NIS")
				self._structure[a[x]]._pdbMapping[interaction] = ( chainObj.chain, thisRes.identifier )
				self.logger.debug("{0}: residue {1}-{2} ({3} in sequence) detected as non-interacting surface.".format(
											interaction, chainObj.chain, thisRes.identifier))
			elif thisRes.identifier in interacting_surface:
				self._structure[a[x]].setTag("IS")
				self._structure[a[x]]._pdbMapping[interaction] = ( chainObj.chain, thisRes.identifier )
				self.logger.debug("{0}: residue {1}-{2} ({3} in sequence) detected as interacting surface.".format(
											interaction, chainObj.chain, thisRes.identifier))

	def mapResiduesToGenome(self, exons):
		j = 0
		for exonStart,exonEnd in exons:
			for i in range(exonStart, exonEnd, 3):
				self._structure[j].setGenomicPosition(i)
				j += 1

	def getAlteredRegions(self, otherIsoform):
		for a in self._structure:
			if a.genomicPosition not in [ y.genomicPosition for y in otherIsoform._structure]:
				a.setIsoformSpecific(True)

	def report(self):
		seq 	= ""
		isoSp 	= ""
		for aa in [ [ Counter(x._tag).most_common(1), x._isoformSpecific] for x in self._structure ]:
			if not aa[0]: 	seq += "*"
			else:		
				if aa[0][0][0] == "NIS":	seq += "S"
				elif aa[0][0][0] == "IS": 	seq += "I"
				elif aa[0][0][0] == "B": 	seq += "B"

			if aa[1]: isoSp += "X"
			else: isoSp += "-"

		self.logger.info("{0}".format(seq))
		self.logger.info("{0}".format(isoSp))
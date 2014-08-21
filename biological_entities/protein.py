from biological_entities import aminoacid
from libs import utils
from libs.SBI.structure import PDB
from libs.SBI.structure.contacts import Complex

from Bio import pairwise2
import logging

class Protein:
	def __init__(self, tx, txInfo):

		self._tx 				= tx
		self._gene 				= txInfo["gene_id"]
		self._uniprot 			= txInfo["Uniprot"]
		self._sequence 			= txInfo["proteinSequence"]
		self._structure			= []
		self._residueCorresp	= {}
		self._has_pdbs			= False
		self._features 			= []

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

	def calculateVolumes(self, interaction, chain):

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
				logging.debug("{0}: residue {1}-{2} ({3} in sequence) detected as buried.".format(
											interaction, chainObj.chain, thisRes.identifier, x))
			elif thisRes.identifier in non_interacting_surface:
				self._structure[a[x]].setTag("NIS")
				self._structure[a[x]]._pdbMapping[interaction] = ( chainObj.chain, thisRes.identifier, "NIS" )
				logging.debug("{0}: residue {1}-{2} ({3} in sequence) detected as non-interacting surface.".format(
											interaction, chainObj.chain, thisRes.identifier, x))
			elif thisRes.identifier in interacting_surface:
				self._structure[a[x]].setTag("IS")
				self._structure[a[x]]._pdbMapping[interaction] = ( chainObj.chain, thisRes.identifier, "IS" )
				logging.debug("{0}: residue {1}-{2} ({3} in sequence) detected as interacting surface.".format(
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

	def checkInteractome3DStructures(self):

		noInteractions = True

		for line in utils.readTable("Data/Databases/Interactome3D/2014_01/interactions.dat"):
			pdbFile = "Data/Databases/Interactome3D/2014_01/interactions/" + line[21]

			if self.uniprot == line[0]:
				logging.debug("Relevant interaction for {0} at {1}.".format(self.tx, pdbFile))
				if self.calculateVolumes(pdbFile, "A"):
					noInteractions = False
			elif self.uniprot == line[1]:
				logging.debug("Relevant interaction for {0} at {1}.".format(self.tx, pdbFile))
				if self.calculateVolumes(pdbFile, "B"):
					noInteractions = False

		if noInteractions:
			logging.debug("No relevant structures found for {0}, {1}.".format(self.tx,self.uniprot))
			return False

		return True

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

		return (seq,isoSp)

	def printPDBInfo(self):
		for pdb in set([ y for x in self._structure for y in x._pdbMapping ]):
			
			chainsSet = set([ x._pdbMapping[pdb][0] for x in self._structure if pdb in x._pdbMapping ])
			chain = set(chainsSet).pop()
			if len(chainsSet) > 1:
				logging.error("More than one chain for a protein in the PDB {0}, chains {1}.".format(pdb, chainsSet))
			isoSpec = [ x._pdbMapping[pdb][1][:-1] for x in self._structure if pdb in x._pdbMapping and x.isoformSpecific ]
			interact = [ x._pdbMapping[pdb][1][:-1] for x in self._structure if pdb in x._pdbMapping and x._pdbMapping[pdb][2]=="IS" ]

			if not isoSpec or not interact:
				return

			logging.debug("{0}, load {1}".format(self._tx,pdb))
			logging.debug("{0}, set bg_rgb=[1,1,1]".format(self._tx))
			logging.debug("{0}, hide all".format(self._tx))
			logging.debug("{0}, show cartoon, all".format(self._tx))
			logging.debug("{0}, color black, all".format(self._tx))
			logging.debug("{0}, color white, chain {1}".format(self._tx,chain))
			logging.debug("{0}, select interact, chain {1} & resi {2}".format(self._tx,chain, "+".join(interact) ))
			logging.debug("{0}, select isoSpecific, chain {1} & resi {2}".format(self._tx,chain, "+".join(isoSpec) ))
			logging.debug("{0}, color orange, isoSpecific".format(self._tx))
			logging.debug("{0}, color red, interact & isoSpecific".format(self._tx))

	def readInterpro(self,interproOut):

		for cols in utils.readTable(interproOut, header=False):
			# https://code.google.com/p/interproscan/wiki/OutputFormats#Tab-separated_values_format_%28TSV%29
		    protein_accession 	= cols[0] #Protein Accession
		    md5_digest 			= cols[1] #Sequence MD5 digest
		    seq_length 			= cols[2] #Sequence Length
		    analysis 			= cols[3] #Analysis
		    signature_accession = cols[4] #Signature Accession
		    signature_descript 	= cols[5] #Signature Description
		    start 				= int(cols[6]) #Start location
		    stop 				= int(cols[7]) #Stop location
		    try:
		    	score 			= float(cols[8]) #Score - is the e-value of the match reported by member database method
		    except ValueError:
		    	score 			= None
		    status 				= True if cols[9] == "T" else False #Status - is the status of the match (T: true)
		    date 				= cols[10] #Date - is the date of the run
		    if len(cols) > 11:
		    	interpro_annotation	= cols[11] #(InterPro annotations - accession)
		    if len(cols) > 12:
		    	interpro_descript	= cols[12] #(InterPro annotations - description)
		    if len(cols) > 13:
		    	go_annotation		= cols[13] #(GO annotations)
		    if len(cols) > 14:
		    	pathway_annotation 	= cols[14] #(Pathways annotations)

			if score and score > 0.01: continue

			isoSpecificRes = set([ x._num for x in self._structure if x.isoformSpecific ])
			featureRes = set(range(start,stop+1))
			overlapRes = isoSpecificRes | featureRes

			featInfo = {}
			featInfo["region"] 			= [start,stop]
			featInfo["accession"] 		= signature_accession
			featInfo["description"] 	= signature_descript
			featInfo["analysis"]		= analysis
			featInfo["percentAffected"]	= len(overlapRes)/len(featureRes) * 100

			self._features.append(featInfo)
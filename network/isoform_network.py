import network
from libs import utils
from libs import options

import pandas as pd
import numpy as np
import abc

class IsoformNetwork(network.Network):
	"""docstring for IsoformNetwork
	IsoformNetwork contains a network of isoforms.

	Node information:
		id(str) 					Transcript Id
		gene_id(str) 				Gene Id of the parent gene.
		exonStructure(list,None)	List of lists, each of them containing the limits of an exon.
		txCoords(list,None) 		List with the starting and the ending genome positions
									of the trancript.
		cdsCoords(list,None) 		List with the starting and the ending genome positions of the CDS.
		strand(str,None)			Strand.
		chr(str,None)				Chromosome.
		median_TPM_N(float,None) 	Median TPM of the isoform in the normal patients.
		median_PSI_N(float,None) 	Median PSI of the isoform in the normal patients.
		median_TPM_T(float,None) 	Median TPM of the isoform in the tumor patients.
		median_PSI_T(float,None) 	Median PSI of the isoform in the tumor patients.
		iLoopsFamily(str,None) 		Loop pattern.
		proteinSequence(str,None)	Protein sequence.
		Uniprot(str,None)			Associated UniprotId.

	Edge information:
		Id1(str) 						Transcript id of interactor 1.
		Id2(str) 						Transcript id of interactor 2.
	"""

	__metaclass__ = abc.ABCMeta

	def __init__(self, name):
		network.Network.__init__(self, name)

	@abc.abstractmethod
	def genenameFilter(self, **kwds):
		raise NotImplementedError()

	def add_node(self, tx, gene_full_name):
		
		if tx in self.nodes():
			self.logger.debug("Transcript {0}, gene {1} already in the network.".format(tx, gene_full_name))
			return True

		geneID = self.genenameFilter( full_name=gene_full_name )[0]

		if geneID is None:
			self.logger.debug("No gene could be extracted for transcript {0}, gene {1}.".format(tx, gene_full_name))
			return False
		
		return self._net.add_node( 	tx, 
									gene_id			= geneID,
									exonStructure	= None,
									txCoords		= None,
									cdsCoords		= None,
									strand 			= None,
									chr 			= None,
									median_TPM_N	= None, 
									median_PSI_N	= None, 
									median_TPM_T	= None, 
									median_PSI_T	= None, 
									iLoopsFamily 	= None, 
									proteinSequence	= None,
									Uniprot 		= None)

	def update_node(self, tx, key, value):
		return self._update_node(tx, key, value)

	def add_edge(self, tx1, tx2):
		self._net.add_edge(tx1,tx2)

	def update_edge(self, tx1, tx2, key, value):
		return self._update_edge(tx1, tx2, key, value)

	def importTranscriptome(self):

		for line in utils.readTable(options.Options().qout + "expression_normal.tsv"):
			gene 		= line[0]
			tx 			= line[1]
			median_PSI 	= float(line[2]) if(line[2] != "NA") else None
			median_TPM 	= float(line[5]) if(line[5] != "NA") else None
			
			if not self.add_node(tx, gene): continue
			if median_PSI is not None: self.update_node( tx, "median_PSI_N", median_PSI )
			if median_TPM is not None: self.update_node( tx, "median_TPM_N", median_TPM )

		for line in utils.readTable(options.Options().qout + "expression_tumor.tsv"):
			gene 		= line[0]
			tx 			= line[1]
			median_PSI 	= float(line[2]) if(line[2] != "NA") else None
			median_TPM 	= float(line[5]) if(line[5] != "NA") else None

			if not self.add_node(tx, gene): continue
			if median_PSI is not None: self.update_node( tx, "median_PSI_T", median_PSI )
			if median_TPM is not None: self.update_node( tx, "median_TPM_T", median_TPM )

		# currentLoopFamily 	= ""
		# for line in utils.readTable("Data/TCGA/UnifiedFasta_" + options.Options().iLoopsVersion + "_loopFamilies.txt", header=False):
		# 	txName = ""
		# 	if ">" in line[0]:
		# 		currentLoopFamily = line[0][1:]
		# 		txName = line[1]
		# 	else:
		# 		txName = line.pop()
			
		# 	if txName in self.nodes(): 
		# 		self.update_node(txName, "iLoopsFamily", currentLoopFamily)

	def readTranscriptInfo(self):
		self.logger.debug("Reading transcript info: exon, CDS and UTR structure.")
		for line in utils.readTable("{0}Data/{1}/knownGene.txt".format(options.Options().wd, options.Options().inputType), header=False):
			if line[0] not in self.nodes(): continue

			tx			= line[0]
			chrom		= line[1] 
			strand		= line[2] 
			txStart		= int(line[3])
			txEnd		= int(line[4])
			cdsStart	= int(line[5])
			cdsEnd		= int(line[6])
			exonCount	= int(line[7])
			exonStarts	= map(int, filter(None, line[8].split(",") ) )
			exonEnds	= map(int, filter(None, line[9].split(",") ) )
			proteinID	= line[10]
			alignID		= line[11]

			self.update_node(tx, "exonStructure", [])
			for i in range(0, len(exonStarts)):
				exon = [ exonStarts[i], exonEnds[i] ]
				self.update_node(tx, "exonStructure", exon)

			self.update_node(tx, "txCoords", [txStart, txEnd])
			self.update_node(tx, "cdsCoords", [cdsStart, cdsEnd])
			self.update_node(tx, "strand", strand)
			self.update_node(tx, "chr", chrom)

		self.logger.debug("Reading transcript info: protein sequence, Uniprot and iLoops family.")
		with open("{0}Data/{1}/UnifiedFasta_{2}.fa".format(options.Options().wd, options.Options().inputType, options.Options().iLoopsVersion)) as FASTA:
			txName 			= ""
			geneFullName 	= ""
			sequence 		= ""
			iLoopsFamily 	= ""
			Uniprot 		= ""

			for line in FASTA:
				if ">" in line:
					if txName:
						if sequence:		self.update_node(txName, "proteinSequence", sequence)	
						if Uniprot: 		self.update_node(txName, "Uniprot", Uniprot[:-1])
						if iLoopsFamily: 	self.update_node(txName, "iLoopsFamily", iLoopsFamily)
						
					elements = line[1:].strip().split("#")
					txName 			= elements[0]
					geneFullName 	= elements[1]
					Uniprot 		= elements[2]
					iLoopsFamily 	= elements[3]
					sequence 		= ""

				else:
					sequence += line.strip()


	def iterate_transcripts(self):
		'''
		Iterate transcripts.
		'''

		txsByGene = {}

		for tx,info in self.nodes(data=True):
			txsByGene.setdefault(info["gene_id"],[])

			txsByGene[info["gene_id"]].append(tx)

		for gene in txsByGene:
			yield gene,txsByGene[gene]
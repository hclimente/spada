from libs import utils
from libs import options
import network

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
		Id1(str) 					Transcript id of interactor 1.
		Id2(str) 					Transcript id of interactor 2.
	"""

	__metaclass__ = abc.ABCMeta

	def __init__(self, name):
		network.Network.__init__(self, name)

	@abc.abstractmethod
	def genenameFilter(self, **kwds):
		raise NotImplementedError()

	def add_node(self,tx,genesFullName):
		
		if tx in self.nodes():
			self.logger.debug("Transcript {0}, gene {1} already in the network.".format(tx, genesFullName))
			return True

		genes = set()
		for g in genesFullName:
			gene = self.genenameFilter(full_name=g)[0]
			if gene != None:
				genes.add(gene)

		if not genes:
			self.logger.debug("No gene could be extracted for transcript {0}, gene {1}.".format(tx, genesFullName))
			return False
		
		self._net.add_node( tx, 
							gene_id			= genes,
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

		return True

	def update_node(self, tx, key, value):
		return self._update_node(tx, key, value)

	def add_edge(self, tx1, tx2):
		self._net.add_edge(tx1,tx2)

	def update_edge(self, tx1, tx2, key, value):
		return self._update_edge(tx1, tx2, key, value)

	def importTranscriptome(self):
		# create transcripts from expression info
		for line in utils.readTable(options.Options().qout + "transcript_expression.tsv",header=False):
			gene = line[0]
			tx = line[1]
			median_TPM_n = float(line[2]) if(line[2] != "NA") else None
			median_TPM_t = float(line[3]) if(line[3] != "NA") else None
			median_PSI_n = float(line[4]) if(line[4] != "NA") else None
			median_PSI_t = float(line[5]) if(line[5] != "NA") else None

			if not self.add_node(tx, gene): 
				continue

			if median_TPM_n is not None: self.update_node( tx, "median_TPM_N", median_TPM_n )
			if median_TPM_t is not None: self.update_node( tx, "median_TPM_T", median_TPM_t )
			if median_PSI_n is not None: self.update_node( tx, "median_PSI_N", median_PSI_n )
			if median_PSI_t is not None: self.update_node( tx, "median_PSI_T", median_PSI_t )

		# exon and CDS info
		for line in utils.readTable("{}data/{}/annotation.gaf".format(options.Options().wd,options.Options().annotation),header=False):
			if line[2]!="transcript" or line[1] not in self.nodes() or line[7]!="GRCh37":
				continue

			entryNumber 		 = line[0]
			featureID 			 = line[1]
			featureType 		 = line[2]
			featureDBSource 	 = line[3]
			featureDBVersion 	 = line[4]
			featureDBDate 		 = line[5]
			featureSeqFileName 	 = line[6]
			composite 			 = line[7]
			compositeType 		 = line[8]
			compositeDBSource 	 = line[9]
			compositeDBVersion 	 = line[10]
			compositeDBDate 	 = line[11]
			alignmentType 		 = line[12]
			featureCoordinates 	 = line[13]
			compositeCoordinates = line[14]
			gene 				 = None
			geneLocus 			 = None
			featureAliases	 	 = None
			featureInfo  		 = None
			if len(line) > 15:
				gene 				 = line[15]
				geneLocus 			 = line[16]
				if len(line) > 17:
					featureAliases	 = line[17]
					if len(line) > 18:
						featureInfo  = line[18]

			chromosome,exons,strand=compositeCoordinates.split(":")
			exons = [ [ int(y) for y in x.split("-") ] for x in exons.split(",")  ]
			txStart = min([ x[0] for x in exons ])
			txEnd = max([ x[1] for x in exons ])

			self.update_node(tx, "strand", strand)
			self.update_node(tx, "chr", chromosome)
			self.update_node(tx,"exonStructure",exons)
			self.update_node(tx, "txCoords", [txStart, txEnd])

			if featureInfo and "CDSstart" in featureInfo:

				cds = [ x for x in featureInfo.split(";") if "CDSstart" in x  or "CDSstop" in x ]

				cdsStart = int([ x.split("=")[1] for x in cds if "CDSstart" in x ][0]) - 1
				cdsEnd = int([ x.split("=")[1] for x in cds if "CDSstop" in x ][0]) - 1
				self.update_node(tx, "cdsCoords", [txStart+cdsStart, txStart+cdsEnd])

		self.logger.debug("Reading transcript info: protein sequence, Uniprot and iLoops family.")
		with open("{}data/{}/sequences.uniprot.loops.fa".format(options.Options().wd, options.Options().annotation)) as FASTA:
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

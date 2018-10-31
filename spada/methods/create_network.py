from spada.bio.transcript import Transcript
from spada.bio.protein import Protein
from spada.io import io
from spada.methods import method
from spada.network import ucsc_gene_network, ucsc_transcript_network
from spada.network import gencode_gene_network, gencode_transcript_network
from spada.network import ensembl_gene_network, ensembl_transcript_network
from spada.io.error import SpadaError

from itertools import product
from networkx import get_node_attributes
import pickle
import numpy as np

class CreateNetwork(method.Method):
	def __init__(self, name, annotation = 'annotation.pklz', new = True):

		if new:
			method.Method.__init__(self, __name__, None)

			if annotation == "ucsc":
				self._genes = ucsc_gene_network.UCSCGeneNetwork(name)
				self._txs = ucsc_transcript_network.UCSCTranscriptNetwork(name)
			elif annotation == "gencode":
				self._genes = gencode_gene_network.GENCODEGeneNetwork(name)
				self._txs = gencode_transcript_network.GENCODETranscriptNetwork(name)
			elif annotation == "ensembl":
				self._genes = ensembl_gene_network.ENSEMBLGeneNetwork(name)
				self._txs = ensembl_transcript_network.ENSEMBLTranscriptNetwork(name)
			else:
				raise SpadaError("Unrecognized annotation: {}.".format(annotation))
		else:
			method.Method.__init__(self, __name__, annotation)
			self._genes._name = name
			self._txs._name = name

		self._new = new

	def run(self, gtf, sequences, ppi, ddi, features, aberrant):

		# read annotation
		self.createNetworks(gtf)
		self.addAberrant(aberrant)

		# isoform data
		self.getIsoformSequences(sequences)
		self.getIsoformFeatures(features)
		self.getInteractions(ppi)
		self.getDomainInteractions(ddi)

		# QC and save
		self.check()

	def createNetworks(self, gtf):

		if not self._new:
			if gtf:
				raise SpadaError("gtf provided when previous network is to be used.")
			else:
				return

		self.logger.info("Importing genes and transcripts from GTF.")

		txLines = ['transcript','exon','CDS','start_codon','stop_codon']

		for line in io.readGTF(gtf):

			if line["feature"] == "gene" and self._genes.accept(line):
				self._genes.add_node(gene_id = line["gene_id"], gene_symbol = line["gene_name"])
			elif line["feature"] in txLines and self._txs.accept(line):
				if line["feature"] == "transcript":
					self._txs.add_node(line["transcript_id"], line["gene_id"])
					self._txs.update_node(line["transcript_id"], "txCoords", [int(line["start"]), int(line["end"]) ])
					self._txs.update_node(line["transcript_id"], "strand", line["strand"] )
					self._txs.update_node(line["transcript_id"], "chr", line["chromosome"] )
					self._txs.update_node(line["transcript_id"], "main", self._txs.isMain(line) )
				elif line["feature"] == "exon":
					self._txs.update_node(line["transcript_id"], "exons", [int(line["start"]), int(line["end"]) ])
				elif line["feature"] == "start_codon":
					pos = line["start"] if line['strand'] == '+' else line['end']
					self._txs.update_node(line["transcript_id"], "start_codon", pos)
				elif line["feature"] == "stop_codon":
					pos = line["start"] if line['strand'] == '+' else line['end']
					self._txs.update_node(line["transcript_id"], "stop_codon", pos)
				elif line["feature"] == "CDS" and self._txs.acceptCDS(line):
					self._txs.update_node(line["transcript_id"], "CDS", [int(line["start"]), int(line["end"]) ])

	def getInteractions(self, ppi):

		if self._new and not ppi:
			raise SpadaError("A file containing the protein-protein interactions must be provided.")
		elif not ppi:
			self.logger.info("Protein-protein interactions from the provided network will be used.")
			return

		self.logger.info("Building protein-protein interaction network.")
		symbols = [ y["symbol"] for x,y in self._genes.nodes(data=True) ]

		for line in io.readPSIMITAB(ppi):

			if line["organismA"][0]["id"] != "9606" or line["organismB"][0]["id"] != "9606":
				continue

			symbolA = [ x["id"] for x in line["symbolA"] if x["type"] == 'entrez gene/locuslink' ]
			symbolB = [ x["id"] for x in line["symbolB"] if x["type"] == 'entrez gene/locuslink' ]

			added = False
			for _ in range(2):
				for A,B in product(symbolA, symbolB):
					if self._genes.add_edge(symbol1 = A, symbol2 = B):
						added = True
						break

				if added:
					break
				
				symbolA.extend([ x["id"] for x in line["aliasA"] if x.get("extra", None) in ["gene name", "gene name synonym"] and x["id"] in symbols ])
				symbolB.extend([ x["id"] for x in line["aliasB"] if x.get("extra", None) in ["gene name", "gene name synonym"] and x["id"] in symbols ])

	def getDomainInteractions(self, ddi):

		if self._new and not ddi:
			raise SpadaError("A file containing the domain-domain interactions must be provided.")
		elif not ddi:
			self.logger.info("Domain-domain interactions from the provided network will be used.")
			return

		self.logger.info("Building isoform-isoform interaction network.")
		allDDIs = { frozenset([x['Pfam1'],x['Pfam2']]) for x in io.readTable(ddi, keys = ['Pfam1','Pfam2']) }
		gene2tx = io.getGene2Tx(self._txs)

		for gene1, gene2 in self._genes._net.edges():
			for tx1,tx2 in product(gene2tx.get(gene1, set()), gene2tx.get(gene2, set())):

				possibleDDIs = { frozenset(x) for x in product(self._txs[tx1]["Pfam"], self._txs[tx2]["Pfam"]) }
				matches = possibleDDIs & allDDIs

				if matches:
					self._txs.add_edge(tx1, tx2, ddi = matches)

	def getIsoformSequences(self, proteins):

		if self._new and not proteins:
			raise SpadaError("A FASTA file with protein sequences file must be provided.")
		elif not proteins:
			self.logger.info("Protein sequences from the provided network will be used.")
			return

		self.logger.info("Reading protein sequences.")

		for tx, sequence in io.readFasta(proteins):
			self._txs.update_node(tx, "proteinSequence", sequence)

	def getIsoformFeatures(self, features):

		if self._new and not features:
			raise SpadaError("A file with the protein features must be provided.")
		elif not features:
			self.logger.info("Protein features from the provided network will be used.")
			return

		self.logger.info("Reading isoform features.")

		featureFields = ['Transcript','Feature_type', 'Feature', 'Start', 'End']
		for line in io.readTable(features, keys = featureFields):

			tx = line['Transcript']
			featureType = line['Feature_type']
			feature = line['Feature']
			start = int(line['Start'])
			end = int(line['End'])

			self._txs.update_node(tx, featureType, (start,end), feature)

	def addAberrant(self, aberrant):

		if aberrant:

			self.logger.info("Import aberrant isoforms absent in GTF.")

			prev = self._txs.skip_filter
			self._txs.skip_filter = True

			for line in io.readTable(aberrant, keys = ['GeneId','Transcript']):
				self._txs.add_node(line['Transcript'], line['GeneId'])
				self._txs.update_node(line["Transcript"], "canonical", False )

			self._txs.skip_filter = prev

	def check(self):

		self.logger.info("Performing consistency checks.")
		for tx, info in self._txs.transcripts():
			Transcript(tx, info)

			if info['CDS']:
				Protein(tx, info)

		self.logger.info("Save annotation into annotation.pklz.")
		self.saveNetworks(filename = 'annotation.pklz')

if __name__ == '__main__':
	pass

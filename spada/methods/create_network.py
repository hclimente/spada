from spada.biological_entities.transcript import Transcript
from spada.biological_entities.protein import Protein
from spada.io import io
from spada.methods import method
from spada.network import ucsc_gene_network, ucsc_transcript_network
from spada.network import gencode_gene_network, gencode_transcript_network
from spada.io.io import SpadaError

from itertools import product
from networkx import get_node_attributes
from numpy import exp2
import pickle
import pandas as pd

import pkg_resources

class CreateNetwork(method.Method):
	def __init__(self, tumor, annotation, new = True):

		if new:
			method.Method.__init__(self, __name__, None, None)

			if annotation == "ucsc":
				self._genes = ucsc_gene_network.UCSCGeneNetwork()
				self._txs = ucsc_transcript_network.UCSCTranscriptNetwork()
			elif annotation == "gencode":
				self._genes = gencode_gene_network.GENCODEGeneNetwork()
				self._txs = gencode_transcript_network.GENCODETranscriptNetwork()
			else:
				raise SpadaError("Unrecognized annotation: {}.".format(annotation), self.logger)
		else:
			if annotation in ["ucsc", "gencode"]:
				genes = pkg_resources.resource_filename('spada', 'data/{}_genes.pkl'.format(annotation))
				txs = pkg_resources.resource_filename('spada', 'data/{}_transcripts.pkl'.format(annotation))
				method.Method.__init__(self, __name__, genes, txs)
			else:
				method.Method.__init__(self, __name__, True, True)

		self._new = new
		self._genes.tumor = tumor
		self._txs.tumor = tumor

	def run(self, gtf, normalExpression, tumorExpression, minExpression, sequences, ppi, ddi, drivers, features, aberrant):

		# read annotation
		self.createNetworks(gtf)
		self.addAberrant(aberrant)

		# gene data
		self.getInteractions(ppi)

		drivers, specificDrivers = self.readDrivers(drivers)
		if drivers:
			self._genes.update_nodes("driver", self.symbol2ids(drivers))
		if specificDrivers:
			self._genes.update_nodes("specificDriver", self.symbol2ids(specificDrivers))

		# isoform data
		for expression,origin in zip([normalExpression, tumorExpression], ["N", "T"]):
			self.measureExpression(expression, minExpression, origin)

		self.getIsoformSequences(sequences)
		self.getIsoformFeatures(features)
		self.getDomainInteractions(ddi)

		# QC and save
		self.check()

	def createNetworks(self, gtf):

		if not self._new:
			if gtf:
				raise SpadaError("gtf provided when previous network is to be used.", self.logger)
			else:
				return()

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
					self._txs.update_node(line["transcript_id"], "start_codon", int(pos))
				elif line["feature"] == "stop_codon":
					pos = line["start"] if line['strand'] == '+' else line['end']
					self._txs.update_node(line["transcript_id"], "stop_codon", int(pos))
				elif line["feature"] == "CDS" and self._txs.acceptCDS(line):
					self._txs.update_node(line["transcript_id"], "CDS", [int(line["start"]), int(line["end"]) ])

	def measureExpression(self, expression, minExpression, origin):

		if self._new and not expression:
			raise SpadaError("An expression file must be provided.", self.logger)
		elif not expression:
			self.logger.info("Expression from the provided network will be used.")
			return()

		self.logger.info("Reading {} expression and calculating PSI.".format(origin))
		expression = pd.read_csv(expression, sep = '\t', index_col = 0)

		# filter out readings on excluded transcripts
		txs = [ tx for tx, txInfo in self._txs.transcripts() if tx in expression.index ]
		expression = expression.loc[txs]

		medianExpression = self.readExpression(expression)
		medianPSI = self.calculatePSI(expression)
		expressed = self.isExpressed(medianExpression, minExpression)

		self._txs.update_nodes("median_TPM_" + origin, medianExpression)
		self._txs.update_nodes("median_PSI_" + origin, medianPSI)
		self._genes.update_nodes("expressedTxs" + origin, expressed)

	def readExpression(self, expression):

		medians = expression.median(axis = 1)
		medians = medians.to_dict()

		return(medians)

	def calculatePSI(self, expression):

		tx2gene = get_node_attributes(self._txs._net, 'gene_id')

		psi = lambda x: exp2(x) / exp2(x).sum()

		psis = expression.groupby(tx2gene).transform(psi)
		medians = psis.median(axis = 1)
		medians = medians.to_dict()

		return(medians)

	def isExpressed(self, expression, threshold):

		expressed = {}

		for tx, median in expression.items():
			if median > threshold:
				gene = self._txs._net.node[tx]["gene_id"]

				expressed.setdefault(gene, set())
				expressed[gene].add(tx)

		return(expressed)

	def getInteractions(self, ppi):

		if self._new and not ppi:
			raise SpadaError("A file containing the protein-protein interactions must be provided.", self.logger)
		elif not ppi:
			self.logger.info("Protein-protein interactions from the provided network will be used.")
			return()

		self.logger.info("Building protein-protein interaction network.")
		symbols = [ y["symbol"] for x,y in self._genes.nodes(data=True) ]

		for line in io.readPSIMITAB(ppi):

			if line["organismA"][0]["id"] != "9606" or line["organismB"][0]["id"] != "9606":
				next

			symbolA = [ x["id"] for x in line["symbolA"] if x["type"] == 'entrez gene/locuslink' ]
			symbolA.extend([ x["id"] for x in line["aliasA"] if x.get("extra", None) in ["gene name", "gene name synonym"] and x["id"] in symbols ])
			symbolB = [ x["id"] for x in line["symbolB"] if x["type"] == 'entrez gene/locuslink' ]
			symbolB.extend([ x["id"] for x in line["aliasB"] if x.get("extra", None) in ["gene name", "gene name synonym"] and x["id"] in symbols ])

			for A in symbolA:
				for B in symbolB:
					self._genes.add_edge(symbol1 = A, symbol2 = B)

	def getDomainInteractions(self, ddi):

		if self._new and not ddi:
			raise SpadaError("A file containing the domain-domain interactions must be provided.", self.logger)
		elif not ddi:
			self.logger.info("Domain-domain interactions from the provided network will be used.")
			return()

		self.logger.info("Building isoform-isoform interaction network.")
		allDDIs = { frozenset([x['Pfam1'],x['Pfam2']]) for x in io.readTable(ddi, keys = ['Pfam1','Pfam2']) }
		
		for gene1, gene2 in self._genes._net.edges():
			txs1 = [ (t,i["Pfam"]) for t,i in self._txs.nodes(data=True) if i["gene_id"] == gene1 and i["Pfam"] ]
			txs2 = [ (t,i["Pfam"]) for t,i in self._txs.nodes(data=True) if i["gene_id"] == gene2 and i["Pfam"] ]

			for tx1,d1 in txs1:
				for tx2,d2 in txs2:
					possibleDDIs = { frozenset(x) for x in product(d1, d2) }
					matches = possibleDDIs & allDDIs

					if matches:
						self._txs.add_edge(tx1, tx2, ddi = matches)

	def readDrivers(self, drivers):

		if self._new and not drivers:
			raise SpadaError("A file containing the tumor drivers must be provided.", self.logger)
		elif not drivers:
			self.logger.info("Drivers from the provided network will be used.")
			return(None,None)

		self.logger.info("Reading known driver genes.")

		allDrivers = {}
		specificDrivers = {}

		for line in io.readTable(drivers, keys = ['gene_id','tumor']):

			allDrivers[line['gene_id']] = True

			if line['tumor'] == self._genes.tumor:
				specificDrivers[line['gene_id']] = True

		return allDrivers, specificDrivers

	def getIsoformSequences(self, proteins):

		if self._new and not proteins:
			raise SpadaError("A FASTA file with protein sequences file must be provided.", self.logger)
		elif not proteins:
			self.logger.info("Protein sequences from the provided network will be used.")
			return()

		self.logger.info("Reading protein sequences.")

		for tx, sequence in io.readFasta(proteins):
			self._txs.update_node(tx, "proteinSequence", sequence)

	def getIsoformFeatures(self, features):

		if self._new and not features:
			raise SpadaError("A file with the protein features must be provided.", self.logger)
		elif not features:
			self.logger.info("Protein features from the provided network will be used.")
			return()

		self.logger.info("Reading isoform features.")

		featureFields = ['tx','featureType', 'feature', 'start', 'end']
		for line in io.readTable(features, keys = featureFields):

			tx = line['tx']
			featureType = line['featureType']
			feature = line['feature']
			start = int(line['start'])
			end = int(line['end'])

			self._txs.update_node(tx, featureType, (start,end), feature)

	def addAberrant(self, aberrant):

		if aberrant:

			self.logger.info("Import aberrant isoforms absent in GTF.")

			for line in io.readTable(aberrant, keys = ['gene','tx']):
				self._txs.add_node(line['tx'], line['gene'])

	def check(self):

		self.logger.info("Quality checks and save networks.")
		for tx, info in self._txs.transcripts():
			Transcript(tx, info)

			if info['CDS']:
				Protein(tx, info)

		self._genes.saveNetwork("genes.pkl")
		self._txs.saveNetwork("transcripts.pkl")

	def symbol2ids(self, symbols):

		ids = {}

		for symbol, value in symbols.items():
			gene_id = [ g for g,i in self._genes.nodes(data=True) if i["symbol"] == symbol ]

			if gene_id:
				ids[gene_id[0]] = value

		return(ids)

if __name__ == '__main__':
	pass

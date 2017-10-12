from spada.interface import out_network
from spada import utils
from spada.methods import method
from spada.network import ucsc_gene_network, ucsc_isoform_network
from spada.network import gencode_gene_network, gencode_isoform_network

from networkx import get_node_attributes
import pickle
import pandas as pd

class CreateNetwork(method.Method):
	def __init__(self, tumor, annotation):

		method.Method.__init__(self, __name__, None, None, None)

		if annotation == "ucsc":
			self._genes = ucsc_gene_network.UCSCGeneNetwork()
			self._txs = ucsc_isoform_network.UCSCIsoformNetwork()
		elif annotation == "gencode":
			self._genes = gencode_gene_network.GENCODEGeneNetwork()
			self._txs = gencode_isoform_network.GENCODEIsoformNetwork()
		else:
			self.logger.error("Unrecognized annotation: {}.".format(annotation))
			exit()

		self._genes.tumor = tumor
		self._txs.tumor = tumor

	def run(self, gtf, normalExpression, tumorExpression, minExpression, proteins, ppi, drivers):

		self.logger.info("Importing genes and transcripts from GTF.")
		self.createNetworks(gtf)

		self.logger.info("Reading expression and calculating PSI.")
		for exprFile,origin in zip([normalExpression, tumorExpression], ["N", "T"]):
			expression = pd.DataFrame.from_csv(exprFile, sep='\t')
			medianExpression = self.readExpression(expression)
			medianPSI = self.calculatePSI(expression)
			expressed = self.isExpressed(medianExpression, minExpression)

			self._txs.update_nodes("median_TPM_" + origin, medianExpression)
			self._txs.update_nodes("median_PSI_" + origin, medianPSI)
			self._genes.update_nodes("expressedTxs" + origin, expressed)

		self.logger.info("Building protein-protein interaction network.")
		self.getInteractions(ppi)

		self.logger.info("Reading known driver genes.")
		drivers, specificDrivers = self.readDrivers(drivers)
		self._genes.update_nodes("driver", drivers)
		self._genes.update_nodes("specificDriver", specificDrivers)

		self.logger.info("Saving networks.")
		self._genes.saveNetwork("genes.pkl")
		self._txs.saveNetwork("transcripts.pkl")

	def createNetworks(self, gtf):

		for line in utils.readGTF(gtf):

			if line["feature"] == "gene":
				self._genes.add_node(gene_id = line["gene_id"], gene_symbol = line["gene_name"])
			elif line["feature"] == "transcript":
				self._txs.add_node(line["transcript_id"], line["gene_id"])
				self._txs.update_node(line["transcript_id"], "txCoords", [line["start"], line["end"]] )
				self._txs.update_node(line["transcript_id"], "strand", line["strand"] )
				self._txs.update_node(line["transcript_id"], "chr", line["chromosome"] )
			elif line["feature"] == "exon":
				self._txs.update_node(line["transcript_id"], "exonStructure", [line["start"], line["end"]] )
			elif line["feature"] == "CDS":
				self._txs.update_node(line["transcript_id"], "cdsCoords", [line["start"], line["end"]] )
			else:
				pass

	def readExpression(self, expression):

		medians = expression.median(axis = 1)
		medians = medians.to_dict()

		return(medians)

	def calculatePSI(self, expression):

		tx2gene = get_node_attributes(self._txs._net, 'gene_id')

		psi = lambda x: x / x.sum()

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

		edges = set()

		for line in utils.readPSIMITAB(ppi):

			if line["organismA"][0]["id"] != "9606" or line["organismB"][0]["id"] != "9606":
				next

			symbolA = [ x["id"] for x in line["symbolA"] if "extra" in x and x["extra"] == "gene name" ]
			symbolB = [ x["id"] for x in line["symbolB"] if "extra" in x and x["extra"] == "gene name" ]

			if symbolA and symbolB:
				self._genes.add_edge(symbol1 = symbolA[0], symbolB = symbolB[0])

		return(edges)

	def readDrivers(self, drivers):

		allDrivers = {}
		specificDrivers = {}

		for line in utils.readTable(drivers):
			gene_id = line[0]
			tumor   = line[1]

			allDrivers[gene_id] = True

			if tumor == self._genes.tumor:
				specificDrivers[gene_id] = True

		return allDrivers, specificDrivers

if __name__ == '__main__':
	pass

from spada.interface import out_network
from spada import utils
from spada.methods import method
from spada.network import ucsc_gene_network, ucsc_isoform_network
from spada.network import gencode_gene_network, gencode_isoform_network

from networkx import get_node_attributes
import pickle
import pandas as pd

class CreateNetwork(method.Method):
	def __init__(self, annotation):

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

	def run(self, gtf, normalExpression, tumorExpression, ppi, drivers):

		self.readGTF(gtf)

		for exprFile,origin in zip([normalExpression, tumorExpression], ["N", "T"]):
			expression = pd.DataFrame.from_csv(exprFile, sep='\t')
			self.updateNodes(self._txs, "median_TPM_" + origin, self.readExpression(expression))
			self.updateNodes(self._txs, "median_PSI_" + origin, self.calculatePSI(expression))

		self.annotateGenes(ppi, drivers)

		self._genes.saveNetwork("genes.pkl")
		self._txs.saveNetwork("transcripts.pkl")

	def readGTF(self, gtf):

		self.logger.info("Importing genes and transcripts from GTF.")
		for line in self.readGTFline(gtf):

			if line["feature"] == "gene":
				self._genes.add_node(gene_id = line["gene_id"], gene_symbol = line["gene_name"])
			elif line["feature"] == "transcript":
				self._txs.add_node(line["transcript_id"], line["gene_id"])
				self._txs.update_node(line["transcript_id"], "txCoords", [line["start"], line["end"]] )
				self._txs.update_node(line["transcript_id"], "strand", line["strand"] )
				self._txs.update_node(line["transcript_id"], "chr", line["chromosome"] )
			elif line["feature"] == "exon":
				self._txs.update_node(line["transcript_id"], "exonStructure", [line["start"], line["end"]] )
			else:
				pass

	def readGTFline(self, gtf):

		for line in utils.readTable(gtf, header=False):

			# read common fields
			parsed = {}
			parsed["chromosome"] = line[0]
			parsed["source"] 	 = line[1]
			parsed["feature"] 	 = line[2]
			parsed["start"] 	 = line[3]
			parsed["end"] 		 = line[4]
			parsed["score"] 	 = line[5]
			parsed["strand"] 	 = line[6]
			parsed["phase"] 	 = line[7]

			# read additional fields
			for field in line[8].strip().split(";"):
				if len(field):
					field = field.strip().split(" ")
					parsed[field[0]] = field[1].strip('"')

			yield parsed

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

	def annotateGenes(self, ppi, drivers):

		self._genes.getInteractome(ppi)
		self._genes.getDriverInfo(drivers)

	def updateNodes(self, net, field, vals):
		for node, value in vals.items():
			net.update_node(node, field, value )

if __name__ == '__main__':
	pass

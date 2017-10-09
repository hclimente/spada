from interface import out_network
from libs import options
from libs import utils
from methods import method
from network import ucsc_gene_network, ucsc_isoform_network

import pickle

class CreateNetwork(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):

		self.createGeneNetwork()
		self.createTranscriptNetwork()

	def createGeneNetwork(self):

		self.logger.info("Creating gene network.")

		if options.Options().annotation == "ucsc":
			self._gene_network = ucsc_gene_network.UCSCGeneNetwork()
		else:
			self.logger.error("Unrecognized input type {0}.".format(options.Options().annotation))
			exit()

		self.logger.debug("Reading gene info.")

		utils.cmd('/soft/R/R-3.2.3/bin/Rscript',
			'pipeline/methods/get_expression.r',
			"{}data/{}/rawdata/{}_iso_tpm_paired-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_iso_tpm_tumor-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_iso_psi_paired-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_iso_psi_tumor-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}transcript_expression.tsv".format(options.Options().qout))

		self._gene_network.readGeneInfo()
		self._gene_network.importSpecificDrivers()

		self._gene_network.importKnownInteractions()

		self._gene_network.saveNetwork("geneNetwork.pkl")

	def createTranscriptNetwork(self,recover=False):

		if options.Options().annotation == "ucsc":
			self._transcript_network = ucsc_isoform_network.UCSCIsoformNetwork()
		else:
			self.logger.error("Unrecognized input type {0}.".format(options.Options().annotation))
			exit()

		self.logger.info("Creating transcript network.")
		self._transcript_network.importTranscriptome()

		self._transcript_network.saveNetwork("txNetwork.pkl")

if __name__ == '__main__':
	pass

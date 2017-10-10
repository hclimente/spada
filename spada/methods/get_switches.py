from interface import out_network
from spada import options
from spada import utils
from methods import method
from network import ucsc_gene_network, ucsc_isoform_network

import pickle

class GetSwitches(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):

		self.createGeneNetwork()
		self.createTranscriptNetwork()

		if not options.Options().externalSwitchesFile:
			switchesFile = self.calculateSwitches()

			self._gene_network.importCandidates(switchesFile)
			self._gene_network.calculateCompatibilityTable()
		
		else:
			self._gene_network.importExternalCandidates()

			# copy random structural_analysis, as it is computationally expensive
			utils.cmd("cp",
					  "{}analyses/{}/structural_analysis/*random*".format(options.Options().wd,options.Options().parentTag),
					  "{}structural_analysis".format(options.Options().qout))

			# copy random switches
			utils.cmd("cp",
					  "{}analyses/{}/randomgenes.pkl".format(options.Options().wd,options.Options().parentTag),
					  options.Options().qout)

		self._gene_network.saveNetwork("genes.pkl")

		#out_network.outputGTF(self._gene_network,self._transcript_network)
		out_network.outCandidateList(self._gene_network,self._transcript_network)

	def calculateSwitches(self):

		switchesFile = "{}candidateList.tsv".format(options.Options().qout)

		self.logger.info("Calculating switches.")
		utils.cmd('/soft/R/R-3.2.3/bin/Rscript', 
			'pipeline/methods/calculate_switches.r', 
			"{}data/{}/rawdata/{}_iso_tpm_paired-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_iso_tpm_tumor-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_gene_read_paired-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_gene_read_tumor-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_iso_psi_paired-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_iso_psi_tumor-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			switchesFile)

		return switchesFile

	def createGeneNetwork(self):

		self.logger.info("Creating gene network.")

		if options.Options().externalSwitchesFile and options.Options().parentTag:
			self._gene_network = pickle.load(open("{}analyses/{}/genes.pkl".format(options.Options().wd,options.Options().parentTag),"rb"))
			self._gene_network.createLogger()
			self._gene_network.cleanNetwork()
		else:
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

		self._gene_network.saveNetwork("genes.pkl")

	def createTranscriptNetwork(self,recover=False):
		if options.Options().externalSwitchesFile and options.Options().parentTag:
			self._transcript_network = pickle.load(open("{}analyses/{}/transcripts.pkl".format(options.Options().wd,options.Options().parentTag),"r"))
			self._transcript_network.createLogger()
		else:
			if options.Options().annotation == "ucsc":
				self._transcript_network = ucsc_isoform_network.UCSCIsoformNetwork()
			else:
				self.logger.error("Unrecognized input type {0}.".format(options.Options().annotation))
				exit()

			self.logger.info("Creating transcript network.")
			self._transcript_network.importTranscriptome()
		
		self._transcript_network.saveNetwork("transcripts.pkl")

if __name__ == '__main__':
	pass

from interface import export_to_MSAnalysis
from interface import out_network
from libs import options
from libs import utils
from methods import method
from network import ucsc_gene_network, ucsc_isoform_network

import cPickle

class GetSwitches(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):

		if not options.Options().externalSwitchesFile:
			switchesFile = self.calculateSwitches()
			externalSwitches = False
		else:
			switchesFile = options.Options().externalSwitchesFile
			externalSwitches = True

			# copy random structural_analysis, as it is computationally expensive
			utils.cmd("cp",
					  "{0}testResults/{1}/{2}/structural_analysis/*random*".format(options.Options().wd,options.Options().annotation,options.Options().parentTag),
					  "{0}structural_analysis".format(options.Options().qout))

			# copy random switches
			utils.cmd("cp",
					  "{0}testResults/{1}/{2}/randomGeneNetwork.pkl".format(options.Options().wd,options.Options().annotation,options.Options().parentTag),
					  options.Options().qout)

		self.createGeneNetwork(switchesFile,externalSwitches=externalSwitches)
		self.createTranscriptNetwork(externalSwitches=externalSwitches)

		out_network.outputGTF(self._gene_network,self._transcript_network)
		out_network.outCandidateList(self._gene_network,self._transcript_network)

		options.Options().printToFile(initialStep="get-relevant-switches")

	def calculateSwitches(self):

		switchesFile = "{}candidateList.tsv".format(options.Options().qout)

		self.logger.info("Calculating switches.")
		utils.cmd('/soft/R/R-3.2.3/bin/Rscript', 
			'pipeline/methods/calculate_switches.r', 
			"{}data/{}/rawdata/{}_iso_tpm_paired-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_iso_tpm_tumor-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_gene_tpm_paired-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_gene_tpm_tumor-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_iso_psi_paired-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			"{}data/{}/rawdata/{}_iso_psi_tumor-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
			switchesFile)

		return switchesFile

	def createGeneNetwork(self,switchesFile,externalSwitches):

		self.logger.info("Creating gene network.")

		if externalSwitches and options.Options().parentTag:
			self._gene_network = cPickle.load(open("{}results/{}/geneNetwork.pkl".format(options.Options().wd,options.Options().parentTag),"r"))
			self._gene_network.createLogger()
			self._gene_network.cleanNetwork()
		else:
			if options.Options().annotation == "TCGA": 
				self._gene_network = ucsc_gene_network.UCSCGeneNetwork()
			else:
				self.logger.error("Unrecognized input type {0}.".format(options.Options().annotation))
				exit()
		
			self.logger.debug("Reading gene info.")

			utils.cmd('/soft/R/R-3.2.3/bin/Rscript', 
				'pipeline/methods/get_expression.r', 
				"{}data/{}/rawdata/{}_iso_tpm_paired-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
				"{}data/{}/rawdata/{}_iso_tpm_tumor-filtered.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag),
				"{}transcript_expression.tsv".format(options.Options().qout))

			self._gene_network.readGeneInfo()
			if options.Options().specificDrivers:
				self._gene_network.importSpecificDrivers()
			
			self._gene_network.importKnownInteractions()

		if not externalSwitches:
			self._gene_network.importCandidates(switchesFile)
			self._gene_network.calculateCompatibilityTable()
		else:
			self._gene_network.importExternalCandidates(switchesFile)

		self._gene_network.saveNetwork("geneNetwork.pkl")

	def createTranscriptNetwork(self,externalSwitches,recover=False):
		if externalSwitches and options.Options().parentTag:
			self._transcript_network = cPickle.load(open("{0}testResults/{1}/{2}/txNetwork.pkl".format(options.Options().wd,options.Options().annotation,options.Options().parentTag),"r"))
			self._transcript_network.createLogger()
		else:
			if options.Options().annotation == "TCGA":
				self._transcript_network = ucsc_isoform_network.UCSCIsoformNetwork()
			else:
				self.logger.error("Unrecognized input type {0}.".format(options.Options().annotation))
				exit()

			self.logger.info("Creating transcript network.")
			self._transcript_network.importTranscriptome()
			self._transcript_network.readTranscriptInfo()
		
		self._transcript_network.saveNetwork("txNetwork.pkl")

if __name__ == '__main__':
	pass
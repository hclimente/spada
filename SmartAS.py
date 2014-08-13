#!/soft/devel/python-2.7/bin/python

from interface import standarize_input
from interface import out_network
from libs import options
from libs import utils
from methods import analyze_interactions
from methods import network_analysis
from methods import structural_analysis
from network import ucsc_gene_network, ucsc_isoform_network

import cPickle
import logging

class SmartAS:
	def __init__(self):
		self.logger = logging.getLogger()

		self.logger.info("SmartAS - Finding significant AS events")
		self.logger.info("Hector Climente - GRIB 2014")

		self._gene_network = None
		self._transcript_network = None
		self._gene_subnetwork = None

	def importData(self):
		self.logger.info("Importing data to a compatible format.")
		standarize_input.standarizeInput()
		self.logger.info("Done. Relaunch please.")

	def exploreData(self):
		self.logger.info("Reading and summarizing input files: computing PSI values and intereplicate agreement.")
		utils.cmd("Pipeline/methods/ExploreData.r", options.Options().qout, 
				  "Data/Input/{0}/{1}/".format(options.Options().inputType, options.Options().tag) )

	def getCandidates(self):
		self.logger.info("Extracting transcripts with high variance and high expression.")
		utils.cmd( "Pipeline/methods/GetCandidates.r", options.Options().minExpression, 
				   options.Options().qout, options.Options().unpairedReplicates )

		self.createTranscriptNetwork(False)
		self.createGeneNetwork(False)

		sortedNodes = sorted( self._gene_network.nodes(data=True), key=lambda (a,dct): dct['score'], reverse=True)
		out_network.outputGTF(sortedNodes, self._transcript_network )
		out_network.outCandidateList(self._gene_network, self._transcript_network)

	def networkAnalysis(self, onlyExperimental):
		
		n = network_analysis.NetworkAnalysis( self._gene_network, self._transcript_network, self._gene_subnetwork )
		n.run(onlyExperimental=onlyExperimental)
		self._gene_subnetwork = n.getGeneSubnetwork(1)

		self._gene_subnetwork.saveNetwork("geneSubnetwork.pkl")
		self._gene_network.saveNetwork("geneNetwork.pkl")
		self._transcript_network.saveNetwork("txNetwork.pkl")

	def launchiLoops(self):

		self.logger.info("Selecting isoforms suitable for {0}.".format( options.Options().iLoopsVersion) )
		utils.pickUniqPatterns(self._transcript_network, self._gene_network)

		self.logger.info("Sending list to Gaudi and performing the iLoops analysis.")
		gaudiThread = utils.cmdOut(
									"ssh", "hectorc@gaudi", \
									"'{0}Pipeline/methods/calculate_interactions.py {1} {2} {3} {4}'".format(
											options.Options().gwd,
											options.Options().gwd,
											options.Options().inputType,
											options.Options().gout,
											options.Options().iLoopsVersion 
										 ), 
									">{0}calculateInteractions.log".format(options.Options().qout)
								  )

	def analyzeInteractions(self):

		a = analyze_interactions.AnalyzeInteractions( self._gene_network, self._transcript_network, self._gene_subnetwork )
		a.run()
		
		self._gene_network.saveNetwork("geneNetwork.pkl")
		self._transcript_network.saveNetwork("txNetwork.pkl")

	def structuralAnalysis(self):

		s = structural_analysis.StructuralAnalysis( self._gene_network, self._transcript_network, self._gene_subnetwork )
		s.run()
		
		self._gene_network.saveNetwork("geneNetwork.pkl")
		self._transcript_network.saveNetwork("txNetwork.pkl")

	def createGeneNetwork(self, recover=False):

		if recover:
			self.logger.info("Recovering gene network from file.")
			self._gene_network = cPickle.load(open(options.Options().qout + "geneNetwork.pkl", "r"))
			self._gene_network.createLogger()
		else:
			self.logger.info("Creating gene network.")

			if options.Options().inputType == "TCGA": 
				self._gene_network = ucsc_gene_network.UCSCGeneNetwork()
			else:
				self.logger.error("Unrecognized input type {0}.".format(options.Options().inputType))
				exit()
			
			self.logger.debug("Reading gene info.")
			self._gene_network.readGeneInfo()
			self._gene_network.importDiffExpression()
			if options.Options().specificDrivers:
				self._gene_network.importSpecificDrivers()
			
			self._gene_network.importKnownInteractions()
			self._gene_network.importCandidates()

			for switch in [ p["isoformSwitches"] for x,p in self._gene_network.nodes(data=True) ]:
				nInfo = self._transcript_network.net.node[switch.nTx]
				tInfo = self._transcript_network.net.node[switch.tTx]
				switch.addTxs(nInfo,tInfo)
				switch.addIsos(nInfo,tInfo)

			self._gene_network.saveNetwork("geneNetwork.pkl")

	def createTranscriptNetwork(self, recover=False):
		if recover:
			self.logger.info("Recovering transcript network from file.")
			self._transcript_network = cPickle.load(open(options.Options().qout + "txNetwork.pkl", "r"))
			self._transcript_network.createLogger()
		else:
			if options.Options().inputType == "TCGA":
				self._transcript_network = ucsc_isoform_network.UCSCIsoformNetwork()
			else:
				self.logger.error("Unrecognized input type {0}.".format(options.Options().inputType))
				exit()

			self.logger.info("Creating transcript network.")
			self._transcript_network.importTranscriptome()
			self._transcript_network.readTranscriptInfo()
			self._transcript_network.saveNetwork("txNetwork.pkl")

	def createGeneSubnetwork(self):
		self.logger.info("Recovering GUILD gene subnetwork from file.")
		self._gene_subnetwork = cPickle.load(open(options.Options().qout + "geneSubnetwork.pkl", "r"))
		self._gene_subnetwork.createLogger()

if __name__ == '__main__':

	logging.basicConfig(
							level=logging.DEBUG,
							format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
							datefmt='%m-%d %H:%M',
	                    	filename=options.Options().qout + 'smartAS.log',
	                    	filemode='w'
					   )

	console = logging.StreamHandler()
	console.setLevel(logging.INFO)

	formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
	console.setFormatter(formatter)
	logging.getLogger().addHandler(console)

	S = SmartAS()

	if options.Options().initialStep <= 0:
		S.importData()
		exit()
	else:
		utils.setEnvironment()

	if options.Options().initialStep <= 1:
		S.exploreData()
	if options.Options().initialStep <= 2:
	 	S.getCandidates()
		if not options.Options().external:
			options.Options().printToFile(initialStep=3)
			exit()
	else:
		S.createTranscriptNetwork(True)
		S.createGeneNetwork(True)

	if options.Options().initialStep <= 3:
		S.networkAnalysis(True)
	else:
		S.createGeneSubnetwork()
	if options.Options().initialStep <= 4:
		S.launchiLoops()
	if options.Options().initialStep <= 5:
		S.structuralAnalysis()
	if options.Options().initialStep <= 6:
		S.analyzeInteractions()
		#S.networkAnalysis(False)
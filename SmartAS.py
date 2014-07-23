#!/soft/devel/python-2.7/bin/python

from interface import out_network
from libs import utils
from libs import options
from methods import analyze_interactions
from methods import network_analysis
from methods import structural_analysis
from network import ucsc_gene_network, ucsc_isoform_network

import cPickle
import logging

import pdb

class SmartAS:
	def __init__(self):
		self.logger = logging.getLogger()

		self.logger.info("SmartAS - Finding significant AS events")
		self.logger.info("Hector Climente - GRIB 2014")

	def exploreData(self):
		self.logger.info("Reading and summarizing input files: computing PSI values and intereplicate agreement.")
		utils.cmd(	
					"Pipeline/methods/ExploreData.r", 
					options.Options().qout, 
					"Data/Input/{0}/{1}/".format(options.Options().inputType, options.Options().tag)
				 )

	def getCandidates(self):
		self.logger.info("Extracting transcripts with high variance and high expression.")
		utils.cmd(
					"Pipeline/methods/GetCandidates.r", 
					options.Options().minExpression, 
					options.Options().qout, 
					options.Options().unpairedReplicates
				 )

		self.createGeneNetwork(False)
		self.createTranscriptNetwork(False)

		out_network.outputGTF(
				sorted(self._gene_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True),
				self._transcript_network
				)
		out_network.outCandidateList(self._gene_network, self._transcript_network)

	def launchiLoops(self):
		self.logger.info("Selecting isoforms suitable for {0}".format( options.Options().iLoopsVersion) )
		utils.pickUniqPatterns(self._transcript_network)

		self.logger.info("Sending list to Gaudi and performing the iLoops analysis.")
		gaudiThread = utils.cmdOut(
									"ssh", "hectorc@gaudi", \
									"'{0}Pipeline/methods/CalculateInteractions.py {1} {2} {3}'".format(
											options.Options().gwd,
											options.Options().gwd,
											options.Options().inputType,
											options.Options().iLoopsVersion 
										 ), 
									">{0}calculateInteractions.log".format(options.Options().qout)
								  )


	def analyzeInteractions(self):

		a = analyzeInteractions.AnalyzeInteractions( self._gene_network, self._transcript_network )
		a.run()

	def networkAnalysis(self):
		
		n = networkAnalysis.NetworkAnalysis( self._gene_network, self._transcript_network )
		n.run()

	def structuralAnalysis(self):

		s = structural_analysis.StructuralAnalysis( self._gene_network, self._transcript_network )
		s.run()

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
			if options.Options().specificDrivers:
				self._gene_network.importSpecificDrivers()
			
			self._gene_network.importKnownInteractions()
			self._gene_network.importCandidates()

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

if __name__ == '__main__':

	utils.setEnvironment()

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

	if options.Options().initialStep <= 1:
		S.exploreData()
	if options.Options().initialStep <= 2:
	 	S.getCandidates()
		if not options.Options().external:
			exit()
	else:
		S.createGeneNetwork(True)
		S.createTranscriptNetwork(True)

	# if options.Options().initialStep <= 3:
	# 	S.launchiLoops()
	if options.Options().initialStep <= 4:
		S.structuralAnalysis()
	if options.Options().initialStep <= 5:
		S.analyzeInteractions()
	if options.Options().initialStep <= 6:
		S.networkAnalysis()
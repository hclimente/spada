from interface import out_network
from libs import utils
from libs import options
from methods import analyzeInteractions
from methods import networkAnalysis
from methods import structural_analysis
from network import ucsc_gene_network, ucsc_isoform_network

import cPickle
import logging
import pdb

class SmartAS:
	def __init__(self, loadGenes=False, loadTranscripts=False ):

		logging.basicConfig(
							level=logging.DEBUG,
        		            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                		    datefmt='%m-%d %H:%M',
                    		filename='smartas.log',
                    		filemode='w'
                    	   )

		console = logging.StreamHandler()
		console.setLevel(logging.INFO)

		formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
		console.setFormatter(formatter)
		logging.getLogger('').addHandler(console)

		if loadGenes:
			logging.info("Recovering already created gene network.")
			self._gene_network = cPickle.load(open(options.Options().qout + "geneNetwork.pkl", "r"))

		else:
			logging.info("Creating gene network.")

			if options.Options().inputType == "TCGA":
				self._gene_network = ucsc_gene_network.UCSCGeneNetwork()
			else:
				logging.error("Unrecognized input type {0}.".format(options.Options().inputType))
				exit()
			
			logging.debug("Reading gene info.")
			self._gene_network.readGeneInfo()
			if options.Options().specificDrivers:
				self._gene_network.importSpecificDrivers(drivers_file=options.Options().specificDrivers, otherDrivers = False)
			
			self._gene_network.importKnownInteractions()
			if options.Options().initialStep >= 2:
				self._gene_network.importCandidates()
			self._gene_network.saveNetwork("geneNetwork.pkl")

		if loadTranscripts:
			logging.info("Recovering already created transcript network.")
			self._transcript_network = cPickle.load(open(options.Options().qout + "txNetwork.pkl", "r"))
		else:
			if options.Options().inputType == "TCGA":
				self._transcript_network = ucsc_isoform_network.UCSCIsoformNetwork()
			else:
				logging.error("Unrecognized input type {0}.".format(options.Options().inputType))
				exit()

			logging.info("Creating transcript network.")
			self._transcript_network.importTranscriptome()
			self._transcript_network.saveNetwork("txNetwork.pkl")
			
	def exploreData(self):
		logging.info("Reading and summarizing input files: computing PSI values and intereplicate agreement.")
		utils.cmd(	
					"Pipeline/methods/ExploreData.r", 
					options.Options().qout, 
					"Data/Input/{0}/{1}/".format(options.Options().inputType, options.Options().tag)
				 )

	def getCandidates(self):
		logging.info("Extracting transcripts with high variance and high expression.")
		utils.cmd(
						"Pipeline/methods/GetCandidates.r", 
						options.Options().minExpression, 
						options.Options().qout, 
						options.Options().unpairedReplicates
				   )

		self._gene_network.importCandidates()

	def candidatePrioritization(self):
		logging.info("Prioritizing candidates.")
		utils.cmd(
					"Pipeline/methods/CandidatePrioritization.py", 
					options.Options().qout, 
					options.Options().inputType
				 )

	def launchiLoops(self):
		logging.info("Selecting isoforms suitable for {0}".format( options.Options().iLoopsVersion) )
		utils.pickUniqPatterns( 
								options.Options().gout,
								options.Options().qout,
								options.Options().inputType,
								options.Options().iLoopsVersion,
								options.Options().replicates * 0.1
							  )

		logging.info("Sending list to Gaudi and performing the iLoops analysis.")
		gaudiThread = utils.cmdOut(
									"ssh", "hectorc@gaudi", \
									"'{0}/Pipeline/methods/CalculateInteractions.py {1} {2} {3} {4}'".format(
											options.Options().gwd,
											options.Options().gaudiWd,
											options.Options().out,
											options.Options().inputType,
											options.Options().iLoopsVersion 
										 ), 
									">Results/{0}/calculateInteractions.log".format(options.Options().qout)
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

if __name__ == '__main__':
	S = SmartAS(True, False)

	# if options.Options().initialStep <= 1:
	# 	S.exploreData()
	# if options.Options().initialStep <= 2:
	#  	S.getCandidates()
	# if options.Options().initialStep <= 3:
	# 	S.candidatePrioritization()
	# 	if not options.Options().external:
	# 		exit()
	# if options.Options().initialStep <= 4:
	# 	S.launchiLoops()
	if options.Options().initialStep <= 5:
		S.analyzeInteractions()
	if options.Options().initialStep <= 6:
		S.structuralAnalysis()
	if options.Options().initialStep <= 7:
		S.networkAnalysis()
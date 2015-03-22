#!/soft/devel/python-2.7/bin/python

from interface import standarize_input
from interface import out_network
from libs import options
from libs import utils
from methods import analyze_interactions
from methods import get_random_switches
from methods import get_switches
from methods import neighborhood_analysis
from methods import network_analysis
from methods import structural_analysis
from methods import result_summary

from methods import testMethod

import logging


class SmartAS:
	def __init__(self):

		self.logger = logging.getLogger()

		self.logger.info("SmartAS - Finding significant AS events")
		self.logger.info("Hector Climente - GRIB 2014")

	def importData(self):

		self.logger.info("Importing data to a compatible format.")
		standarize_input.standarizeInput()
		self.logger.info("Done. Relaunch please.")

	def getSwitches(self):

		g = get_switches.GetSwitches(None,None,None)
		g.run()
			
	def networkAnalysis(self,onlyExperimental):
		
		n = network_analysis.NetworkAnalysis(True,True)
		n.run(onlyExperimental=onlyExperimental)
		
		geneSubnetwork = n.getGeneSubnetwork(1)
		geneSubnetwork.saveNetwork("geneSubnetwork.pkl")

	def launchiLoops(self):

		self.logger.info("Selecting isoforms suitable for {0}.".format( options.Options().iLoopsVersion) )
		a = analyze_interactions.AnalyzeInteractions(True,True)
		a.selectIloopsSwitches("Driver")

		self.logger.info("Sending list to Gaudi and performing the iLoops analysis.")
		utils.cmd("ssh", "hectorc@gaudi",
				  "'{0}Pipeline/methods/calculate_interactions.py {1} {2} {3} {4}'".format(
		 		  options.Options().gwd,options.Options().gwd,options.Options().inputType,
		 		  options.Options().gout,options.Options().iLoopsVersion ))

	def analyzeInteractions(self):

		a = analyze_interactions.AnalyzeInteractions(True,True)
		a.run()

	def structuralAnalysis(self):

		s = structural_analysis.StructuralAnalysis(True,True)
		if not options.Options().parallelRange:
			import glob
			files=glob.glob("{0}structural_analysis/interpro_analysis_[0-9]*.tsv".format(options.Options().qout))
			if files:
				s.joinFiles()
			else:
				utils.launchJobs(s._gene_network,'structural_analysis')
		else:
			s.run()

	def neighborhoodAnalysis(self):

		n = neighborhood_analysis.NeighborhoodAnalysis(True,True)
		n.run()

	def summarizeResults(self):

		s = result_summary.ResultSummary(True,True)
		s.run()

	def createRandomSwitches(self):
		
		r = get_random_switches.GetRandomSwitches(True,True)
		r.run()

	def testing(self):

		test = testMethod.Test(True,True)
		test.run()

if __name__ == '__main__':

	logging.basicConfig(level=logging.DEBUG,
						format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
						datefmt='%m-%d %H:%M',
					   	filename='{0}logs/{1}_smartAS_{2}{3}.log'.format(
					   					options.Options().qout,
					   					options.Options().tag,
					   					options.Options().initialStep,
					   					options.Options().filetag),
					   	filemode='w')

	console = logging.StreamHandler()
	console.setLevel(logging.INFO)

	formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
	console.setFormatter(formatter)
	logging.getLogger().addHandler(console)

	S = SmartAS()
	utils.setEnvironment()

	if options.Options().initialStep == "import-data":
		S.importData()
	
	# Get and characterize switches
	elif options.Options().initialStep == "get-switches":
	 	S.getSwitches()

	# analyze switches
	elif options.Options().initialStep == "get-relevant-switches":
		S.structuralAnalysis()
	elif options.Options().initialStep == "random-switches":
		S.createRandomSwitches()
	
	# analyze model switches
	elif options.Options().initialStep == "launch-iloops":
		S.launchiLoops()
	elif options.Options().initialStep == "get-interaction-changes":
		S.analyzeInteractions()
	elif options.Options().initialStep == "experimental-network-analysis":
		S.networkAnalysis(False)
	elif options.Options().initialStep == "predicted-network-analysis":
		S.networkAnalysis(True)
	elif options.Options().initialStep == "neighborhood-analysis":
		S.neighborhoodAnalysis()

	# summarize results
	elif options.Options().initialStep == "summary":
		S.summarizeResults()

	# test commands
	elif options.Options().initialStep == "test":
		S.testing()
	
	S.logger.info("SmartAS will close.")
#!/usr/bin/env python

from interface import out_network
from libs import options
from libs import utils
from methods import get_i3d_broken_interactions
from methods import get_random_switches
from methods import get_switches
from methods import mutation_comparison
from methods import mutation_feature_overlap
from methods import neighborhood_analysis
from methods import network_analysis
from methods import structural_analysis
from methods import result_summary
from methods import wgs_mutations

from methods import testMethod

import logging


class SmartAS:
	def __init__(self):

		self.logger = logging.getLogger()

		self.logger.info("SmartAS - Finding significant AS events")
		self.logger.info("Hector Climente - GRIB 2014-2015")

	def getSwitches(self):

		g = get_switches.GetSwitches(None,None,None)
		g.run()
			
	def networkAnalysis(self,onlyExperimental):
		
		n = network_analysis.NetworkAnalysis(True,True)
		n.run(onlyExperimental=onlyExperimental)
		
		geneSubnetwork = n.getGeneSubnetwork(1)
		geneSubnetwork.saveNetwork("geneSubnetwork.pkl")

	def structuralAnalysis(self):

		s = structural_analysis.StructuralAnalysis(True,True)

		if not options.Options().parallelRange:
			# non-parallelized operation
			
			import glob
			files=glob.glob("{}structural_analysis/interpro_analysis_[0-9]*.tsv".format(options.Options().qout))
			
			if files:
				# all analyses ran already
				s.joinFiles()
				out_network.outCandidateList(s._gene_network,s._transcript_network)
			else:
				# launch analyses
				utils.launchJobs(s._gene_network,'structural_analysis')
		else:
			# analyze a chunk of the switches
			s.run()

	def I3DBrokenInteractions(self):

		i = get_i3d_broken_interactions.GetI3DBrokenInteractions(True,True)
		i.clean()
		i.run()

	def neighborhoodAnalysis(self):

		n = neighborhood_analysis.NeighborhoodAnalysis(True,True)
		n.run()

	def compareSwitchesAndMutations(self):

		m = mutation_comparison.MutationComparison(True,True)
		m.clean()
		m.run()

	def summarizeResults(self):

		s = result_summary.ResultSummary(True,True)
		s.run()

	def createRandomSwitches(self):
		
		r = get_random_switches.GetRandomSwitches(True,True)
		r.run()

	def searchMutationFeatureOverlap(self):
		
		m = mutation_feature_overlap.MutationFeatureOverlap(True,True)
		m.clean()
		m.run()

	def studyWGSMutations(self):
		w = wgs_mutations.WGSMutations(True,True)
		w.clean()
		w.run()

	def testing(self):

		test = testMethod.Test(True,True)
		test.run()

if __name__ == '__main__':

	logging.basicConfig(level=logging.DEBUG,
						format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
						datefmt='%m-%d %H:%M',
					   	filename='{}logs/{}_{}{}.log'.format(options.Options().qout,
					   		options.Options().tag, options.Options().initialStep,
					   		options.Options().filetag), filemode='w')

	console = logging.StreamHandler()
	console.setLevel(logging.INFO)

	formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
	console.setFormatter(formatter)
	logging.getLogger().addHandler(console)

	S = SmartAS()

	# Get and characterize switches
	if options.Options().initialStep == "get-switches":
	 	S.getSwitches()

	# analyze switches
	elif options.Options().initialStep == "get-functional-switches":
		S.structuralAnalysis()
	elif options.Options().initialStep == "random-switches":
		S.createRandomSwitches()
	
	# analyze model switches
	elif options.Options().initialStep == "experimental-network-analysis":
		S.networkAnalysis(False)
	elif options.Options().initialStep == "predicted-network-analysis":
		S.networkAnalysis(True)
	elif options.Options().initialStep == "neighborhood-analysis":
		S.neighborhoodAnalysis()
	elif options.Options().initialStep == "get-i3d-broken-interactions":
		S.I3DBrokenInteractions()
	elif options.Options().initialStep == "mutation-comparison":
		S.compareSwitchesAndMutations()
	elif options.Options().initialStep == "mutation-feature-overlap":
		S.searchMutationFeatureOverlap()
	elif options.Options().initialStep == "wgs-mutations":
		S.studyWGSMutations()

	# summarize results
	elif options.Options().initialStep == "summary":
		S.summarizeResults()

	# test commands
	elif options.Options().initialStep == "test":
		S.testing()
	
	S.logger.info("SmartAS will close.")
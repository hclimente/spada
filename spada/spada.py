#!/usr/bin/env python

from interface import out_network
from libs import options
from libs import utils

# calculation
from methods import create_network
from methods import get_random_switches
from methods import get_switches
from methods import structural_analysis
# candidates
from methods import me_analysis
from methods import me_geneset_analysis
from methods import wes_mutations_feature_overlap
from methods import annotate_switches
from methods import me_ppi
# validation
from methods import explore_pannegative
from methods import candidates_pathways
# data
from methods import result_summary
from methods import neighborhood_analysis
# other
from methods import get_i3d_broken_interactions
from methods import test

import logging
import os

class spada:
	def __init__(self):

		self.logger = logging.getLogger()
		self.logger.info("SPADA - Finding significant splicing changes.")

	def createNetwork(self):

		c.create_network.CreateNetwork(None, None, None)
		c.run()

	# initial
	def getSwitches(self):

		g = get_switches.GetSwitches(None,None,None)
		g.run()

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

	# candidates calculation
	def recurrenceAnalysis(self):

		utils.cmd('Rscript',
			"pipeline/meta/recurrence_analysis.R",
			"{}candidateList_info.tsv".format(options.Options().qout),
			"{}candidateList_recurrence.tsv".format(options.Options().qout))

	def I3DBrokenInteractions(self):

		i = get_i3d_broken_interactions.GetI3DBrokenInteractions(True,True)
		i.clean()
		i.run()

	def neighborhoodAnalysis(self):

		n = neighborhood_analysis.NeighborhoodAnalysis(True,True)
		n.run()

	def studyMutualExclusion(self):

		m = me_analysis.MEAnalysis(True,True)
		m.clean()
		m.run()

	def studyGenesetMutualExclusion(self):
		m = me_geneset_analysis.MEGenesetAnalysis(True,True)
		m.clean()
		m.run()

	def annotateSwitches(self):
		a = annotate_switches.AnnotateSwitches(True,True)
		a.run()

	# validation
	def explorePannegative(self):

		p = explore_pannegative.ExplorePannegative(True,True)
		p.clean()
		p.run()

	def candidatesPathways(self):

		c = candidates_pathways.CandidatesPathways(True,True)
		c.run()

	def summarizeResults(self):

		s = result_summary.ResultSummary(True,True)
		s.clean()
		s.run()

	def createRandomSwitches(self):

		r = get_random_switches.GetRandomSwitches(True,True)
		r.run()

	def studyWESMutationsFeatureOverlap(self):

		m = wes_mutations_feature_overlap.WESMutationsFeatureOverlap(True,True)
		m.clean()
		m.run()

	def testing(self):

		t = test.Test(True,True)
		t.run()

if __name__ == '__main__':

	if not os.path.exists("{}".format(options.Options().qout)):
		utils.cmd("mkdir","-p","{}logs".format(options.Options().qout))

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

	S = spada()

	# Get and characterize switches
	if options.Options().initialStep == "get-switches":
	 	S.getSwitches()

	# analyze switches
	elif options.Options().initialStep == "get-functional-switches":
		S.structuralAnalysis()
	elif options.Options().initialStep == "random-switches":
		S.createRandomSwitches()

	# get candidates
	elif options.Options().initialStep == "recurrence-analysis":
		S.recurrenceAnalysis()
	elif options.Options().initialStep == "wes-mutations-feature-overlap":
		S.studyWESMutationsFeatureOverlap()
	elif options.Options().initialStep == "me-analysis":
		S.studyMutualExclusion()
	elif options.Options().initialStep == "me-geneset-analysis":
		S.studyGenesetMutualExclusion()
	elif options.Options().initialStep == "annotate-candidate-switches":
		S.annotateSwitches()

	# validation
	elif options.Options().initialStep == "explore-pannegative":
		S.explorePannegative()
	elif options.Options().initialStep == "candidates-pathways":
		S.candidatesPathways()

	# analyze model switches
	elif options.Options().initialStep == "neighborhood-analysis":
		S.neighborhoodAnalysis()

	# summarize results
	elif options.Options().initialStep == "summary":
		S.summarizeResults()

	# test commands
	elif options.Options().initialStep == "test":
		S.testing()

	# deprecated
	elif options.Options().initialStep == "get-i3d-broken-interactions":
		S.I3DBrokenInteractions()

	S.logger.info("SPADA will close.")

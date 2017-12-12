#!/usr/bin/env python

from spada.methods import create_network
from spada import utils
from spada.interface import out_network

import argparse
import logging
from math import log2
import os

def structuralAnalysis():

	s = structural_analysis.StructuralAnalysis(True,True)

	if not options.parallelRange:
		# non-parallelized operation

		import glob
		files=glob.glob("{}structural_analysis/interpro_analysis_[0-9]*.tsv".format(options.qout))

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
def recurrenceAnalysis():

	utils.cmd('Rscript',
		"pipeline/meta/recurrence_analysis.R",
		"{}candidateList_info.tsv".format(options.qout),
		"{}candidateList_recurrence.tsv".format(options.qout))

def neighborhoodAnalysis():

	n = neighborhood_analysis.NeighborhoodAnalysis(True,True)
	n.run()

def studyMutualExclusion():

	m = me_analysis.MEAnalysis(True,True)
	m.clean()
	m.run()

def studyGenesetMutualExclusion():
	m = me_geneset_analysis.MEGenesetAnalysis(True,True)
	m.clean()
	m.run()

def annotateSwitches():
	a = annotate_switches.AnnotateSwitches(True,True)
	a.run()

# validation
def explorePannegative():

	p = explore_pannegative.ExplorePannegative(True,True)
	p.clean()
	p.run()

def candidatesPathways():

	c = candidates_pathways.CandidatesPathways(True,True)
	c.run()

def summarizeResults():

	s = result_summary.ResultSummary(True,True)
	s.clean()
	s.run()

def createRandomSwitches():

	r = get_random_switches.GetRandomSwitches(True,True)
	r.run()

def studyWESMutationsFeatureOverlap():

	m = wes_mutations_feature_overlap.WESMutationsFeatureOverlap(True,True)
	m.clean()
	m.run()

def testing():

	t = methods.test.Test(True,True)
	t.run()

parser = argparse.ArgumentParser(prog = "spada.py",
				description = "Find significant alternative splicing switches. Analyze their functional impact.",
				epilog= "Héctor Climente-González, 2014-2017")

parser.add_argument('-wd', '--working-directory', dest='wd', action='store',
					default='.',
					help='Root file of SmartAS folder in the current machine.')
parser.add_argument('-a', '--all-switches', dest='onlyModels', action='store_false',
					help='Only use the model switches.')

subparsers = parser.add_subparsers(help='sub-command help')

################################################
###   INITIALIZE THE NETWORKS               ####
################################################

def createNetwork(o):
	c = create_network.CreateNetwork(o.tumor, o.annotation)
	c.run(o.gtf, o.normalExpression, o.tumorExpression, log2(o.minExpression), o.seq, o.ppi, o.ddi, o.drivers, o.features)

subparser_init = subparsers.add_parser('init', help='Initialize help')
subparser_init.add_argument('-T', '--tumor', dest='tumor' ,action='store',
							help='Identifier of the analysis e.g. the tumor type.')
subparser_init.add_argument('-a', '--annotation', dest='annotation', action='store', default='gencode',
							choices=['ucsc', 'gencode'],
							help='Used annotation.')
subparser_init.add_argument('-g', '--gtf', dest='gtf', action='store',
							help='GTF with the gene, transcript, exon and CDS annotation.')
subparser_init.add_argument('-n', '--expression-normal', dest='normalExpression', action='store',
							help='Tab-separated table with transcript-level expression data of the normal samples \
							in log2(TPM); a row per annotated transcript, indicated in the first column.')
subparser_init.add_argument('-t', '--expression-tumor', dest='tumorExpression', action='store',
							help='Equivalent to --expression-normal for the tumor samples.')
subparser_init.add_argument('-m', '--minimum-expression', dest='minExpression', action='store',
					 		default='0.1', type=float,
					 		help='Minimum expression value, in TPM, to consider a transcript expressed.')
subparser_init.add_argument('-p', '--ppi', dest='ppi', action='store',
							help='File with protein-protein interactions, in PSI-MI TAB format >= 2.5.')
subparser_init.add_argument('-d', '--ddi', dest='ddi', action='store',
							help='Pairs of interacting domains, in TSV format.')
subparser_init.add_argument('-s', '--seq', dest='seq', action='store',
							help='Fasta file with the protein sequences of each transcript.')
subparser_init.add_argument('-f', '--features', dest='features', action='store',
							help='Tab-separated table with transcript-level features.')
subparser_init.add_argument('-D', '--drivers', dest='drivers', action='store',
							help='Tab-separated table containing the genes to be considered tumor drivers. The first column \
							contains the gene symbol, and the second the tumor type where it was detected (which must match \
							--tumor when required).')

subparser_init.set_defaults(task="createNetwork")
subparser_init.set_defaults(func=createNetwork)

################################################
###   IMPORT SWITCHES                       ####
################################################

def getSwitches():
	g = get_switches.GetSwitches()
	g.run(o.switchesFile)

subparser_switches = subparsers.add_parser('switches', help='Calculate switches help')
subparser_switches.add_argument('-s', '--switches', dest='switchesFile', action='store',
								default=None,type=str,
								help='File containing switches as TSV.')

################################################
###   MAIN                                  ####
################################################

options = parser.parse_args()

logging.basicConfig(level=logging.DEBUG,
					format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
					datefmt='%m-%d %H:%M',
				   	filename='.{}.log'.format(options.task), filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.INFO)

formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger().addHandler(console)

logger = logging.getLogger()
logger.info("SPADA - Finding significant splicing changes.")
options.func(options)
logger.info("SPADA will close.")

#!/usr/bin/env python

from spada.methods import create_network
from spada import utils
from spada.interface import out_network

import argparse
import logging
from math import log2

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
				description = "Find significant alternative splicing switches. Analyze their functional impact.")

parser.add_argument('-wd', '--working-directory', dest='wd', action='store',
					default='.',
					help='Root file of spada folder in the current machine.')
parser.add_argument('-a', '--all-switches', dest='onlyModels', action='store_false',
					help='Only use the model switches.')

subparsers = parser.add_subparsers(help='sub-command help')

################################################
###   INITIALIZE THE NETWORKS               ####
################################################

def createNetwork(o):
	if not o.network:
		c = create_network.CreateNetwork(o.tumor, o.annotation)
		c.run(o.gtf, o.normalExpression, o.tumorExpression, log2(o.minExpression), o.seq, o.ppi, o.ddi, o.drivers, o.features)
	else:
		c = create_network.CreateNetwork(o.tumor, o.annotation)
		c.run(o.gtf, o.normalExpression, o.tumorExpression, log2(o.minExpression), o.seq, o.ppi, o.ddi, o.drivers, o.features)

subparser_init = subparsers.add_parser('init', help='Initialize help')

subparser_init.add_argument('-T', '--tumor', dest='tumor', action='store',
							help='Identifier of the analysis e.g. the tumor type.')
subparser_init.add_argument('-N', '--recycle-net', dest='network', action='store_true',
							help='Use previous networks.')
subparser_init.add_argument('-a', '--annotation', dest='annotation', action='store',
							choices=['ucsc', 'gencode'], default=None,
							help='Used annotation. Required for newly generated networks.')
subparser_init.add_argument('-g', '--gtf', dest='gtf', action='store', default=None,
							help='GTF with the gene, transcript, exon and CDS annotation.')
subparser_init.add_argument('-n', '--expression-normal', dest='normalExpression',
							action='store', default=None,
							help='Tab-separated table with transcript-level expression data of the normal samples \
							in log2(TPM); a row per annotated transcript, indicated in the first column.  Required \
							for newly generated networks.')
subparser_init.add_argument('-t', '--expression-tumor', dest='tumorExpression',
							action='store', default=None,
							help='Equivalent to --expression-normal for the tumor samples. Required for newly \
							generated networks.')
subparser_init.add_argument('-m', '--minimum-expression', dest='minExpression',
							action='store', default=None, type=float,
					 		help='Minimum expression value, in TPM, to consider a transcript expressed. Required \
							for newly generated networks.')
subparser_init.add_argument('-p', '--ppi', dest='ppi', action='store', default=None,
							help='File with protein-protein interactions, in PSI-MI TAB format >= 2.5. Required \
							for newly generated networks.')
subparser_init.add_argument('-d', '--ddi', dest='ddi', action='store', default=None,
							help='Pairs of interacting domains, in TSV format. Required for newly generated networks.')
subparser_init.add_argument('-s', '--seq', dest='seq', action='store', default=None,
							help='Fasta file with the protein sequences of each transcript. Required for newly \
							generated networks.')
subparser_init.add_argument('-f', '--features', dest='features', action='store', default=None,
							help='Tab-separated table with transcript-level features. Required for newly \
							generated networks.')
subparser_init.add_argument('-D', '--drivers', dest='drivers', action='store', default=None,
							help='Tab-separated table containing the genes to be considered tumor drivers. The first column \
							contains the gene symbol, and the second the tumor type where it was detected (which must match \
							--tumor when required). Required for newly generated networks.')
subparser_init.add_argument('-A', '--aberrant', dest='aberrant', action='store', default=None,
							help='Tab-separated table containing gene-aberrant transcript pairs. The genes must be in the \
							annotation, but not the transcripts. Hence, they do not have known genomic coordinates, \
							exon structure or CDS. They can, however, present other features (protein sequence, \
							expression, etc.).')

subparser_init.set_defaults(task="createNetwork")
subparser_init.set_defaults(func=createNetwork)

################################################
###   IMPORT SWITCHES                       ####
################################################

def readSwitches():
	g = get_switches.GetSwitches()
	g.run(o.switchesFile)

subparser_switches = subparsers.add_parser('switches', help='Read switches help')
subparser_switches.add_argument('-s', '--switches', dest='switchesFile', action='store', required=True,
								default=None, type=str, help='File containing switches as TSV.')
subparser_switches.set_defaults(task="readSwitches")
subparser_switches.set_defaults(func=readSwitches)

def functionalAnalysis():
	s = structural_analysis.StructuralAnalysis(True,True)
	s.run()

subparser_functional = subparsers.add_parser('function', help='Functional analysis help')
subparser_functional.set_defaults(task="functionalAnalysis")
subparser_functional.set_defaults(func=functionalAnalysis)

################################################
###   MAIN                                  ####
################################################

options = parser.parse_args()

if "task" not in options:
	print("No sub-command selected.")
	parser.print_help()
	exit()

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

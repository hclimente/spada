#!/usr/bin/env python

from methods import me_ppi

import logging
import os

class metaspada:
	def __init__(self):

		self.logger = logging.getLogger()

		self.logger.info("METASPADA - Finding significant AS events")
		self.logger.info("Héctor Climente-González - GRIB 2014-2017")

	def pancancerGetSwitches(self):
		utils.cmd('Rscript',"pipeline/meta/pancancer_calculate_switches.R","candidateList_info")
		utils.cmd('Rscript',"pipeline/meta/pancancer_calculate_switches.R","random.candidateList_info")

	def pancancerRecurrenceAnalysis(self):
		utils.cmd('Rscript',
				  'pipeline/meta/get_switch_origin.R')
		utils.cmd('Rscript',
			"pipeline/meta/recurrence_analysis.R",
			"{}candidateList_info.agg.tsv".format(options.Options().qout),
			"{}candidateList.origin.tsv".format(options.Options().qout),
			"{}candidateList_recurrence.tsv".format(options.Options().qout))

	def pancancerStudyWESMutationsFeatureOverlap(self):
		utils.cmd("mkdir","-p","{}structural_analysis".format(options.Options().qout))
		utils.cmd("mkdir","-p","{}mutations".format(options.Options().qout))
		utils.cmd('Rscript',
				  'pipeline/meta/pancancer_mutated_features_analysis.R')

	def pancancerStudyMutualExclusion(self):
		utils.cmd("mkdir","-p","{}mutations".format(options.Options().qout))
		utils.cmd('Rscript',
				  'pipeline/meta/pancancer_mutual_exclusion_analysis.R')

	def pancancerStudyPPIMutualExclusion(self):
		m = me_ppi.MEPPI(False,False)
		m.clean()
		m.run()

	def pancancerStudyCoocurrence(self):
		utils.cmd('Rscript', 'pipeline/meta/pancancer_coocurrence_analysis.R')

if __name__ == '__main__':

	if not os.path.exists("{}".format(options.Options().qout)):
		utils.cmd("mkdir","-p","{}logs".format(options.Options().qout))

	logging.basicConfig(level=logging.DEBUG,
						format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
						datefmt='%m-%d %H:%M',
					   	filename='{}logs/{}_{}{}.log'.format(options.Options().qout,
					   		options.Options().tag, options.Options().task,
					   		options.Options().filetag), filemode='w')

	console = logging.StreamHandler()
	console.setLevel(logging.INFO)

	formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
	console.setFormatter(formatter)
	logging.getLogger().addHandler(console)

	M = metaspada()

	if options.Options().task == "get-switches":
	 	M.pancancerGetSwitches()
	elif options.Options().task == "recurrence-analysis":
		M.pancancerRecurrenceAnalysis()
	elif options.Options().task == "wes-mutations-feature-overlap":
		M.pancancerStudyWESMutationsFeatureOverlap()
	elif options.Options().task == "me-analysis":
		M.pancancerStudyMutualExclusion()
	elif options.Options().task == "me-ppi":
		M.pancancerStudyPPIMutualExclusion()
	elif options.Options().task == "co-occurence":
		M.pancancerStudyCoocurrence()
	else:
		M.logger.error("Unrecognized task {}".format(options.Options().task))

	M.logger.info("METASPADA will close.")

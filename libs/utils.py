#!/soft/devel/python-2.7/bin/python

from libs import options

import subprocess
from rpy2.robjects import r
import os
import logging

logger = logging.getLogger("utils")

def cmd(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)

	logger.debug(command)
	subprocess.call(command, shell=True)

def cmdOut(*args):
	command = [ str(x) for x in args ]

	logger.debug(command)
	return subprocess.Popen(command, stdout=subprocess.PIPE)

def readTable(path, sep="\t", header=True):
	"""Read a table in a file, and generate a list of strings per row. 
	Skips rows starting with "#".

	sep (str): field separator.
	header (bool): presence of a header, to discard the first row.
	"""
	counter = 0
	with open(path) as FILE:
		for line in FILE:
			if line[0]=="#": 
				continue
			elif header and counter is 0:
				counter = 1
				continue
			
			yield line.strip().split(sep)

def setEnvironment():

	o = options.Options()
	o.printToFile()

	logger.info("Preparing the environment.")

	if not os.path.exists(".testOld/" + o.out):
		cmd("mkdir","-p",".testOld/" + o.out)
	if not os.path.exists(".testOld2/" + o.out):
		cmd("mkdir","-p",".testOld2/" + o.out)

	cmd("rm -r .testOld2/" + o.out )
	cmd("mv", ".testOld/" + o.out, ".testOld2/" + o.out )
	cmd("mv", o.qout, ".testOld/" + o.out )
	cmd("mkdir -p", "testResults/" + o.out + "RWorkspaces")
	cmd("mkdir", o.qout + "DataExploration")
	cmd("mkdir -p", o.qout + "iLoops/" + o.iLoopsVersion)
	cmd("mkdir", o.qout + "GUILD_experimental")
	cmd("mkdir", o.qout + "GUILD_enriched")
	cmd("mkdir", o.qout + "structural_analysis")
	cmd("mkdir", o.qout + "neighborhood_analysis")
	cmd("mkdir", o.qout + "result_summary")

	if o.initialStep == "get-switches":
		
		#Set R workspace
		r('wd <- \"{0}\"'.format(options.Options().wd))
		r('out <- \"{0}\"'.format(options.Options().qout))
		r('inputData <- list()')	
		r('inputData[["Conditions"]] <- c("N", "T")')
		
		reps = ",".join(set([ "\"" + x + "\"" for x in options.Options().replicates ]))
		ureps = ",".join(set([ "\"" + x + "\"" for x in options.Options().unpairedReplicates ]))

		r('inputData[["Replicates"]] <- c({0})'.format (reps))
		r('inputData[["unpairedReplicates"]] <- c({0})'.format (ureps))

		r('inputData[["minExpression"]] <- {0}'.format(options.Options().minExpression))
		r('save.image("' + options.Options().qout + 'RWorkspaces/0_InitialEnvironment.RData")')

	if o.initialStep not in ["import-data"]:
		cmd("cp -r", ".testOld/" + o.out + "DataExploration", o.qout)
		cmd("cp", ".testOld/" + o.out + "RWorkspaces/*.RData", o.qout + "/RWorkspaces")
		
		cmd("cp", ".testOld/" + o.out + "expression_normal.tsv", ".testOld/" + o.out + "expression_tumor.tsv", o.qout)
		cmd("cp", ".testOld/" + o.out + "candidateList.tsv", ".testOld/" + o.out + "candidateList_v2.tsv",".testOld/" + o.out + "candidateList_v3.tsv",".testOld/" + o.out + "expressedGenes.lst", o.qout)
		cmd("cp", ".testOld/" + o.out + "candidates_normal.gtf", ".testOld/" + o.out + "candidates_tumor.gtf", o.qout)
		cmd("cp", ".testOld/" + o.out + "geneNetwork*.pkl", ".testOld/" + o.out + "txNetwork*.pkl", o.qout)
		cmd("cp", ".testOld/" + o.out + "msInput.txt", o.qout)
			
	if o.initialStep != "launch-iloops" and o.initialStep not in ["import-data","get-switches"]:
		cmd("cp", ".testOld/" + o.out + "candidatesGaudi.lst", o.qout)
	if o.initialStep != "neighborhood-analysis" and o.initialStep not in ["import-data","get-switches"]:
		cmd("cp -r", ".testOld/{0}/GUILD_experimental".format(o.out), o.qout)
		cmd("cp", ".testOld/{0}/geneSubnetwork.pkl".format(o.out), o.qout)
	if o.initialStep != "experimental-network-analysis" and o.initialStep not in ["import-data","get-switches"]:
		cmd("cp -r", ".testOld/{0}/neighborhood_analysis".format(o.out), o.qout)
	if o.initialStep != "predicted-network-analysis" and o.initialStep not in ["import-data","get-switches"]:
		cmd("cp -r", ".testOld/{0}/GUILD_enriched".format(o.out), o.qout)
		cmd("cp -r", ".testOld/{0}/iLoops/{1}".format(o.out, o.iLoopsVersion), o.qout)
	if o.initialStep != "get-relevant-switches" and o.initialStep not in ["import-data","get-switches"]:
		cmd("cp -r", ".testOld/{0}/structural_analysis".format(options.Options().out), o.qout)
	if o.initialStep != "summary" and o.initialStep not in ["import-data","get-switches"]:
		cmd("cp -r", ".testOld/{0}/result_summary".format(o.out), o.qout)

def geneclusterLaunch(tag,base,*args):
	command = base
	for argument in args:
		command += " " + argument

	with open(tag+".sh","w") as configFile:
		configFile.write('#!/bin/sh\n')
		configFile.write('#$ -q normal\n')
		configFile.write('#$ -cwd\n')
		configFile.write("#$ -e /data/users/hector/{0}.log\n".format(tag))
		configFile.write("#$ -o /data/users/hector/{0}.log\n".format(tag))

		configFile.write("source ~/.bashrc\n")
		configFile.write(command+"\n")

	cmd("qsub","-N",tag,tag+".sh")

def selectIloopsSwitches(tx_network,gn_network,prop):

	loopFamilies = {}

	for tx,family in [ (x[0],x[1]["iLoopsFamily"]) for x in tx_network.nodes(data=True) if x[1]["iLoopsFamily"] is not None ]:
		if family in loopFamilies:
			loopFamilies[family].add(tx)
		else:
			loopFamilies[family] = set(tx)

	analyzedLoops = {}

	with open(options.Options().qout + "candidatesGaudi.lst", "w") as CANDIDATES_GAUDI:
		sortedNodes = sorted(gn_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene,gInfo in sortedNodes:

			# genes that dont have a property or are related to one who
			#  does e.g. driver, are not considered
			if not gInfo[prop] or not [ x for x,y in gn_network._net.neighbors(gene) if y[prop] ]:
				break

			for switch in gInfo["isoformSwitches"]:
				nIso = switch.nTx
				tIso = switch.tTx
				nInfo = tx_network._net.node[nIso]
				tInfo = tx_network._net.node[tIso]

				analyze = 0
				comment = "To analyze."

				if switch.score < 0.1: 
					continue
				elif abs(gInfo["diffExpression_logFC"]) > 0.5 or gInfo["diffExpression_p"] < 0.05:
					analyze = -1
					comment = "Gene differentially expressed between conditions."
				elif not switch.cds_diff: 
					analyze = -1
					comment = "No CDS change."
				elif not switch.cds_overlap: 
					analyze = -1
					comment = "No overlap between CDS."
				elif not nInfo["iLoopsFamily"] or not tInfo["iLoopsFamily"]:
					analyze = -1
					comment = "No loops mapped by {0}.".format(options.Options().iLoopsVersion)
				elif nInfo["iLoopsFamily"] == tInfo["iLoopsFamily"]:
					analyze = -1
					comment = "No different loops mapped by {0}.".format(options.Options().iLoopsVersion)

				if analyze < 0:
					CANDIDATES_GAUDI.write("{0}\t{1}\t{2}\n".format(nIso, analyze, comment))
					CANDIDATES_GAUDI.write("{0}\t{1}\t{2}\n".format(tIso, analyze, comment))
					continue

				for isoform,thisLoopPattern in zip([nIso,tIso],[nInfo["iLoopsFamily"],tInfo["iLoopsFamily"]]):

					if os.path.isfile("iLoops/{0}/{1}/{2}.tar.gz".format(
											options.Options().inputType,
											options.Options().iLoopsVersion,
											isoform) ):
							analyze = 1
							comment = "Already analyzed."
					elif thisLoopPattern in analyzedLoops:
							analyze = 2
							if isoform in analyzedLoops[thisLoopPattern]:
								comment = "Already being analyzed in this batch."
							else:
								comment = "Analyzing relative {0}.".format(analyzedLoops[thisLoopPattern])
					else:
						for iso in loopFamilies[thisLoopPattern]:
							if os.path.isfile("{0}iLoops/TCGA/{1}/{2}.tar.gz".format(options.Options().wd, options.Options().iLoopsVersion,iso) ):
								analyze = 2
								comment = "Analyzed relative {0}.".format(iso)
								break

					CANDIDATES_GAUDI.write( "{0}\t{1}\t{2}\n".format(isoform, analyze, comment))

					if analyze == 0: 
						analyzedLoops[thisLoopPattern] = isoform

	cmd("ssh","hectorc@gaudi", "'rm -r {0}'".format(options.Options().gout))
	cmd("ssh","hectorc@gaudi","'mkdir -p {0}Output'".format(options.Options().gout))
	cmd("ssh","hectorc@gaudi","'mkdir {0}Input'".format(options.Options().gout))
	cmd("ssh","hectorc@gaudi","'mkdir {0}logs'".format(options.Options().gout))
	cmd("scp", "-r", options.Options().qout + "candidatesGaudi.lst",
		"hectorc@gaudi.imim.es:" + options.Options().gout)

#!/soft/devel/python-2.7/bin/python

from libs import options

from subprocess import call,Popen,PIPE
from rpy2.robjects import r
import urllib2
import os
import logging

import pdb

logger = logging.getLogger("utils")

def cmd(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)

	logger.debug(command)
	call(command, shell=True)

def cmdOut(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)

	logger.debug(command)
	return Popen(command, shell=True, stdout=PIPE)

def readTable(path, sep="\t", header=True):
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
	cmd("rm -r old2/" + o.out )
	cmd("mv", "old/" + o.out, "old2/" + o.out )
	cmd("mv", o.qout, "old/" + o.out )
	cmd("mkdir -p", "Results/" + o.out + "RWorkspaces")
	cmd("mkdir", o.qout + "DataExploration")
	cmd("mkdir -p", o.qout + "iLoops/" + o.iLoopsVersion)
	cmd("mkdir", o.qout + "GUILD_experimental")
	cmd("mkdir", o.qout + "GUILD_enriched")

	if not o.external:
		if o.initialStep <= 1:
			
			#Set R workspace
			r("wd <- \"" + options.Options().wd + "\"")
			r("out <- \"Results/" + options.Options().out + "\"")
			r("inputData <- list()")	
			r('inputData[["Conditions"]] <- c("N", "T")')
			r('inputData[["Replicates"]] <- ' + str(options.Options().replicates))
			r('save.image("' + options.Options().qout + 'RWorkspaces/0_InitialEnvironment.RData")')
			
			#getDB()

		if o.initialStep > 1:

			cmd("cp -r", "old/" + o.out + "/DataExploration", o.qout)
			cmd("cp", "old/" + o.out + "/RWorkspaces/1_ExploreData.RData", o.qout + "/RWorkspaces")
			cmd("cp", "old/" + o.out + "smartAS.log", o.qout + "/RWorkspaces")
			cmd("cp -r", "old/{0}/iLoops".format(o.out), o.qout)
			#cmd("rm -r", "{0}iLoops/{1}".format(o.qout, o.iLoopsVersion), o.qout)

		if o.initialStep > 2:
			cmd("cp", "old/" + o.out + "/RWorkspaces/2_GetCandidates.RData", o.qout + "/RWorkspaces")
			cmd("cp", "old/" + o.out + "/candidateList.tsv", "old/" + o.out + "/candidateList_v2.tsv","old/" + o.out + "/expressedGenes.lst", o.qout)
			cmd("cp", "old/" + o.out + "/candidates_normal.gtf", "old/" + o.out + "/candidates_tumor.gtf", o.qout)
			cmd("cp", "old/" + o.out + "/geneNetwork*.pkl", "old/" + o.out + "/txNetwork*.pkl", o.qout)
			cmd("cp", "old/" + o.out + "/expression_normal.tsv", "old/" + o.out + "/expression_tumor.tsv", o.qout)
	else:
		
		if o.initialStep == 2:
			oriOut = "_".join(o.out.split("_")[:-1])
			cmd("cp", "Results/" + oriOut + "/RWorkspaces/1_ExploreData.RData", "Results/" + o.out + "/RWorkspaces")
			cmd("Pipeline/scripts/InputUnpaired.r", oriOut, o.unpairedReplicates, "Data/Input/" + o.out + "/")
		if o.initialStep == 3:
			cmd("cp", o.external + ".tsv" , "Results/" + o.out + "/candidateList.tsv")
			cmd("cp", o.external + "_expressedGenes.lst", "Results/" + o.out + "/expressedGenes.lst")

	if o.initialStep > 3:
		cmd("cp -r", "old/{0}/GUILD_experimental".format(o.out), o.qout)
		cmd("cp", "old/{0}/geneSubnetwork.pkl".format(o.out), o.qout)
	if o.initialStep > 4:
		cmd("cp", "old/" + o.out + "/candidatesGaudi.lst", o.qout)
	if o.initialStep > 5:
		pass
		#cmd("cp", "old/{0}".format(options.Options().out), o.qout)
	if o.initialStep > 6:
		cmd("cp -r", "old/{0}/GUILD_enriched".format(o.out), o.qout)
		cmd("cp -r", "old/{0}/iLoops/{1}".format(o.out, o.iLoopsVersion), o.qout)

def getDB():
	with open("Data/Databases/Intogen.tsv", "w") as Intogen:
	
		query="""
		DEFINE
			intogen='/data/project/gene',
			genes='https://bitbucket.org/intogen/intogen-sources.git?ensembl/hsa/genes',
			projects='/data/projects'
		ON
			'https://bitbucket.org/intogen/intogen-mutations.git'
		SELECT
			genes (GENE_ID, SYMBOL),
			projects (PROJECT_NAME),
			intogen (SAMPLE_PROP, FM_QVALUE, CLUST_QVALUE)
		FROM
			intogen
		WHERE
			(
				intogen.FM_QVALUE < '0.05'
				OR
				intogen.CLUST_QVALUE < '0.05'
			)
			"""
			
		req = urllib2.Request("http://www.intogen.org/oql")
		res = urllib2.urlopen(req, query)
		Intogen.write(res.read())

def finish(opt):
	
	print("* Moving files to the Results directory and creating a summary tar file.")
	    
	outFolder = "/home/hector/Results/" + opt["out"]
	outTag = opt["tag1"] + "_mE" + str(opt["minExpression"])

	if not cmdOut("mkdir -p", outFolder).stdout.read().strip():
		overwrite = raw_input("\tDirectory exists. Do you want to overwrite it? (y/n)")
		if not overwrite == "y":
			return

	cmd("cp", "Results/" + opt["out"] + "/candidates_normal.gtf", outFolder + "/" + "candidates_normal." + outTag + ".gtf")
	cmd("cp", "Results/" + opt["out"] + "/candidates_tumor.gtf", outFolder + "/" + "candidates_tumor." + outTag + ".gtf")
	cmd("cp", "Results/" + opt["out"] + "/candidateList.top.tsv", outFolder + "/" + "candidateList." + outTag + ".tsv")

	os.chdir(outFolder)
	cmd("cp -r", "../SmartAS/Results/" + opt["out"] + "/* .")

	cmd("tar -czvf", outTag + ".tar.gz candidates_normal." + outTag + ".gtf candidates_tumor." + outTag + ".gtf candidateList." + outTag + ".tsv")
	cmd("rm", "*" + opt["out"] + "*gff", "*" + opt["out"] + "*tsv")
	os.chdir("/home/hector/SmartAS")

def pickUniqPatterns(tx_network, gn_network):

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
			for switch in properties["isoformSwitches"]:
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
				elif switch.cds_diff == "False": 
					analyze = -1
					comment = "No CDS change."
				elif switch.cds_overlap == "False": 
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
	sshTransfer("candidatesGaudi.lst")

def sshTransfer(*args):
	files = ""
	for aFile in args:
		files += " {0}{1}".format( options.Options().qout, aFile )
	cmd("scp","-r", files,"hectorc@gaudi.imim.es:" + options.Options().gout )
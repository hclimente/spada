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

	strCommand = ''
	for arg in command:
		strCommand += " " + arg

	logger.debug(strCommand)
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

	if o.initialStep in ["import-data","get-switches"]:
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

		cmd("cp -r", ".testOld/" + o.out + "DataExploration", o.qout)
		cmd("mv", ".testOld/" + o.out + "RWorkspaces/*.RData", o.qout + "/RWorkspaces")
		
		cmd("cp", ".testOld/{0}expression_normal.tsv".format(o.out), ".testOld/{0}expression_tumor.tsv".format(o.out), o.qout)
		cmd("cp", ".testOld/{0}candidateList.tsv".format(o.out), ".testOld/{0}candidateList_v2.tsv".format(o.out),".testOld/{0}candidateList_v3.tsv".format(o.out),".testOld/{0}expressedGenes.lst".format(o.out), o.qout)
		cmd("cp", ".testOld/{0}candidates_normal.gtf".format(o.out), ".testOld/{0}candidates_tumor.gtf".format(o.out), o.qout)
		cmd("cp", ".testOld/{0}geneNetwork*.pkl".format(o.out), ".testOld/{0}txNetwork*.pkl".format(o.out), o.qout)
		cmd("cp", ".testOld/{0}msInput.txt".format(o.out), o.qout)
			
	if o.initialStep == "get-relevant-switches":
		cmd("rm","-r",".testOld2/{0}/structural_analysis".format(o.out))
		cmd("mv",".testOld/{0}/structural_analysis".format(o.out),".testOld2/{0}/structural_analysis".format(o.out))
		cmd("mv","{0}/structural_analysis".format(o.qout),".testOld/{0}/structural_analysis".format(o.out))
		cmd("mkdir", "{0}/structural_analysis".format(o.qout))
	elif o.initialStep == "neighborhood-analysis":
		cmd("rm","-r",".testOld2/{0}/neighborhood_analysis".format(o.out))
		cmd("mv",".testOld/{0}/neighborhood_analysis".format(o.out),".testOld2/{0}/neighborhood_analysis".format(o.out))
		cmd("mv","{0}/neighborhood_analysis".format(o.qout),".testOld/{0}/neighborhood_analysis".format(o.out))
		cmd("mkdir", "{0}/neighborhood_analysis".format(o.qout))
	elif o.initialStep == "launch-iloops":
		cmd("rm",".testOld2/{0}/candidatesGaudi.lst".format(o.out))
		cmd("mv",".testOld/{0}/candidatesGaudi.lst".format(o.out),".testOld2/{0}/".format(o.out))
		cmd("mv","{0}/candidatesGaudi.lst".format(o.qout),".testOld/{0}/".format(o.out))
	elif o.initialStep == "get-interaction-changes":
		cmd("rm","-r",".testOld2/{0}/interaction_changes".format(o.out))
		cmd("mv",".testOld/{0}/interaction_changes".format(o.out),".testOld2/{0}/interaction_changes".format(o.out))
		cmd("mv","{0}/interaction_changes".format(o.qout),".testOld/{0}/interaction_changes".format(o.out))
		cmd("mkdir", "{0}/interaction_changes".format(o.qout))
	elif o.initialStep == "experimental-network-analysis":
		cmd("rm","-r",".testOld2/{0}/GUILD_experimental".format(o.out))
		cmd("mv",".testOld/{0}/GUILD_experimental".format(o.out),".testOld2/{0}/GUILD_experimental".format(o.out))
		cmd("mv","{0}/GUILD_experimental".format(o.qout),".testOld/{0}/GUILD_experimental".format(o.out))
		cmd("mkdir", "{0}/GUILD_experimental".format(o.qout))

		cmd("rm",".testOld2/{0}/geneSubnetwork.pkl".format(o.out))
		cmd("mv",".testOld/{0}/geneSubnetwork.pkl".format(o.out),".testOld2/{0}/".format(o.out))
		cmd("mv","{0}/geneSubnetwork.pkl".format(o.qout),".testOld/{0}/".format(o.out))
	elif o.initialStep == "predicted-network-analysis":
		cmd("rm","-r",".testOld2/{0}/GUILD_enriched".format(o.out))
		cmd("mv",".testOld/{0}/GUILD_enriched".format(o.out),".testOld2/{0}/GUILD_enriched".format(o.out))
		cmd("mv","{0}/GUILD_enriched".format(o.qout),".testOld/{0}/GUILD_enriched".format(o.out))
		cmd("mkdir", "{0}/GUILD_enriched".format(o.qout))

	elif o.initialStep == "summary":
		cmd("rm","-r",".testOld2/{0}/result_summary".format(o.out))
		cmd("mv",".testOld/{0}/result_summary".format(o.out),".testOld2/{0}/result_summary".format(o.out))
		cmd("mv","{0}/result_summary".format(o.qout),".testOld/{0}/result_summary".format(o.out))
		cmd("mkdir", "{0}/result_summary".format(o.qout))

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

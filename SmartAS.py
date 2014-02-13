#!/soft/devel/python-2.7/bin/python

import sys, getopt
from shutil import copy, copytree
from os import path, chdir
from libsmartas import *

def main(argv):

	print("""
#######################################
#                                     #
#               SmartAS               #
#    Finding significant AS events    #
#                                     #
#      Hector Climente-GRIB 2014      #
#                                     #
#######################################
	""")

	#Set variables
	opts, args = getopt.getopt(argv, "f:")

	cfgFile = ""

	for opt, arg in opts:
		if opt == "-f":
			cfgFile = arg
		
	opt = setEnvironment(cfgFile)

	if opt["initialStep"] <= 1:
		exploreData(opt)
	if opt["initialStep"] <= 2:
		getCandidates(opt)
	if opt["initialStep"] <= 3:
		candidatePrioritization(opt)
	if opt["initialStep"] <= 4:
		prepareILoopsInput(opt)
	if opt["initialStep"] <= 5:
		launchILoops(opt)
	if opt["initialStep"] <= 6:
		pass
		#exloreILoopsResults(opt)
	
	finish(opt)
	#copytree("Results", "../Dropbox/SmartAS")
	#copy("SmartAS.RData", "../Dropbox/SmartAS")

def exploreData(opt):
	
	print("* Reading and summarizing input files: computing PSI values and plotting correlations between replicates.")
	cmd("Pipeline/ExploreData.r")
	copy("SmartAS.RData", "Results/" + opt["out"] + "/RWorkspaces/1_ExploreData.RData")

def getCandidates(opt):

	print("* Extracting transcripts with high variance and high expression.")
	cmd("Pipeline/GetCandidates.r", opt["minExpression"], opt["minCandidateExpression"])

	copy("SmartAS.RData", "Results/" + opt["out"] + "/RWorkspaces/2_GetCandidates.RData")
	cmd("sort", "Results/" + opt["out"] + "/expressedGenes.lst", ">" + "Results/" + opt["out"] + "/expressedGenes.tmp.lst")
	cmd("mv", "Results/" + opt["out"] + "/expressedGenes.tmp.lst", "Results/" + opt["out"] + "/expressedGenes.lst")

	cmd("Pipeline/OutputCandidates.py", "Results/" + opt["out"] + "/candidateList.tsv", opt["out"])
	
def candidatePrioritization(opt):

	print("* Prioritizing candidates.")
	cmd("Pipeline/CandidatePrioritization.py " + opt["out"])

def prepareILoopsInput(opt):

	print("* Retrieving protein sequences for transcripts and printing to multiFASTA file.")
	cmd("Pipeline/GetiLoopsInput.py", opt["out"], "/expressedGenes.lst", "/candidateList.top.tsv")

def launchILoops(opt):

	print("* Launching iLoops jobs.")
	
	cmd("ssh hectorc@gaudi 'rm -r", opt["gOut"] + "/iLoops; mkdir -p", opt["gOut"] + "'")
	cmd("scp -r " + "Results/" + opt["out"] + "/iLoops hectorc@gaudi.imim.es:" + opt["gOut"])
	cmd("ssh hectorc@gaudi '" + opt["gaudiWd"] + "/Pipeline/launchILoops.py", opt["gOut"] + "/iLoops'")

def exloreILoopsResults(opt):

	print("* Examining iLoops results.")
	cmd("Pipeline/exploreiLoopsOutput.py")

main(sys.argv[1:])
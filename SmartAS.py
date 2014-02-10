#!/soft/devel/python-2.7/bin/python

import sys, getopt
from shutil import copy, copytree
from os import path, chdir

#Custom library
from sh import *

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
		exploreData()
	if opt["initialStep"] <= 2:
		getCandidates(opt["minExpression"], opt["minCandidateExpression"], opt["minPSI"])
	if opt["initialStep"] <= 3:
		candidatePrioritization(opt["top"])
	if opt["initialStep"] <= 4:
		prepareILoopsInput()
	if opt["initialStep"] <= 5:
		pass
		#launchILoops()
	if opt["initialStep"] <= 6:
		pass
		#exloreILoopsResults()
	
	finish(opt)
	#copytree("Results", "../Dropbox/SmartAS")
	#copy("SmartAS.RData", "../Dropbox/SmartAS")

def exploreData():
	
	print("* Reading and summarizing input files: computing PSI values and plotting correlations between replicates.")
	cmd("Pipeline/ExploreData.r")
	copy("SmartAS.RData", opt["out"] + "/RWorkspaces/1_ExploreData.RData")

def getCandidates(minExpression, minCandidateExpression, minPSI):

	print("* Extracting transcripts with high variance and high expression.")
	cmd("Pipeline/GetCandidates.r", minExpression, minCandidateExpression, minPSI)

	copy("SmartAS.RData", "Results/RWorkspaces/2_GetCandidates.RData")
	cmd("sort", opt["out"] + "/expressedGenes.lst", ">" + opt["out"] + "/expressedGenes.tmp.lst")
	cmd("mv", opt["out"] + "/expressedGenes.tmp.lst", opt["out"] + "/expressedGenes.lst")

	cmd("Pipeline/OutputCandidates.py", opt["out"] + "/candidateList.tsv")
	
def candidatePrioritization(top):

	print("* Prioritizing candidates.")
	cmd("Pipeline/CandidatePrioritization.py", top)

def prepareILoopsInput():

	getExpressedGenes = 1

	print("* Retrieving protein sequences for transcripts and printing to multiFASTA file.")
	
	diff = cmdOut("diff old/expressedGenes.lst", opt["out"] + "/expressedGenes.lst 2>&1")
	
	if not diff and path.exists("old/iLoops/ExpressedTranscripts.fasta"):
		getExpressedGenes = 0

	cmd("Pipeline/getiLoopsInput.py", opt["out"] + "/expressedGenes.lst", opt["out"] + "/candidateList.top.tsv", getExpressedGenes)

def launchILoops():

	print("* Launching iLoops jobs.")
	
	cmd("ssh hectorc@gaudi 'mv" + opt["out"] + "/iLoops ~/SmartAS/old; rm -r ~/SmartAS/old'")
	cmd("scp -r " + opt["out"] + "/iLoops hectorc@gaudi.imim.es:" + opt["gOut"])
	cmd("ssh hectorc@gaudi '" + opt["gaudiWd"] + "/Pipeline/launchILoops.py " + opt["gOut"] + "/iLoops'")

	print("\t* Waiting...")

def exloreILoopsResults():

	print("* Examining iLoops results.")
	cmd("Pipeline/exploreiLoopsOutput.py")

main(sys.argv[1:])
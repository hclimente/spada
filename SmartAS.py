#!/usr/bin/python

from subprocess import Popen,PIPE
import sys, getopt
from shutil import copy, copytree
import os

import sh

def main(argv):

	print """
#######################################
#                                     #
#               SmartAS               #
#    Finding significant AS events    #
#                                     #
#      Hector Climente-GRIB 2014      #
#                                     #
#######################################
	"""

	#Set variables
	currentStep = 0
	wd = "/home/hector/SmartAS/"
	gaudiWd = "/sbi/users/hectorc/SmartAS/Results/iLoops"
	minExpression = 0
	minCandidateExpression = 4
	minPSI = 0.25

	opts, args = getopt.getopt(argv, "s:wd:")

	for opt, arg in opts:
		if opt == "-s":
			currentStep = int(arg)
		if opt == "-wd":
			wd = arg

	sh.setEnvironment(wd, currentStep)
	
	if currentStep <= 1:
		exploreData()
	if currentStep <= 2:
		getCandidates(minExpression, minCandidateExpression, minPSI)
	if currentStep <= 3:
		prepareILoopsInput()
	if currentStep <= 4:
		launchILoops()
	if currentStep <= 5:
		exloreILoopsResults()
	
	#copytree("Results", "../Dropbox/SmartAS")
	#copy("SmartAS.RData", "../Dropbox/SmartAS")

def exploreData():
	
	print "* Reading and summarizing input files: computing PSI values and plotting correlations between replicates."
	sh.cmd("Pipeline/PSICalculation.r")
	#subprocess.call("Pipeline/PSICalculation.r", shell=True)
	copy("SmartAS.RData", "Results/RWorkspaces/1_ExploreData.RData")

def getCandidates(minExpression, minCandidateExpression, minPSI):

	print "* Extracting transcripts with high variance and high expression."
	sh.cmd("Pipeline/GetCandidates.r", minExpression, minCandidateExpression, minPSI)
	#subprocess.call("Pipeline/GetCandidates.r " + str(minExpression) + " " + str(minCandidateExpression) + " " + str(minPSI), shell=True)
	copy("SmartAS.RData", "Results/RWorkspaces/2_GetCandidates.RData")

	sh.cmd("sort Results/expressedGenes.lst >Results/expressedGenes.tmp.lst")
	#subprocess.call("sort Results/expressedGenes.lst >Results/expressedGenes.tmp.lst", shell=True)
	sh.cmd("mv Results/expressedGenes.tmp.lst Results/expressedGenes.lst")
	#subprocess.call("mv Results/expressedGenes.tmp.lst Results/expressedGenes.lst", shell=True)

def prepareILoopsInput():

	getExpressedGenes = 1

	print "* Retrieving protein sequences for transcripts and printing to multiFASTA file."
	
	diff = Popen('diff old/expressedGenes.lst Results/expressedGenes.lst 2>&1', shell=True, stdout=PIPE)

	if not diff.stdout.read().strip() and os.path.exists("old/iLoops/input/ExpressedTranscripts.fasta"):
		getExpressedGenes = 0

	#subprocess.call("Pipeline/getiLoopsInput.pl Results/expressedGenes.lst Results/candidateList.lst " + getExpressedGenes, shell=True)
	sh.cmd("Pipeline/getiLoopsInput.pl Results/expressedGenes.lst Results/candidateList.lst ", getExpressedGenes)

def launchILoops():

	print "* Launching iLoops jobs."
	
	#subprocess.call("scp -r Results/iLoops hectorc@gaudi.imim.es:~/SmartAS/Results", shell=True)
	sh.cmd("scp -r Results/iLoops hectorc@gaudi.imim.es:~/SmartAS/Results")
	#subprocess.call("ssh hectorc@gaudi.imim.es '~/SmartAS/Pipeline/launchILoops.sh /sbi/users/hectorc/SmartAS/Results/iLoops", shell=True)
	sh.cmd("ssh hectorc@gaudi.imim.es '~/SmartAS/Pipeline/launchILoops.sh /sbi/users/hectorc/SmartAS/Results/iLoops")

	print "\t* Waiting..."

def exloreILoopsResults():

	print "* Examining iLoops results."
	#subprocess.call("Pipeline/exploreiLoopsOutput.py", shell=True)
	sh.cmd("Pipeline/exploreiLoopsOutput.py")

main(sys.argv[1:])

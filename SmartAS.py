#!/soft/devel/python-2.7/bin/python

from subprocess import Popen,PIPE
import sys, getopt
from shutil import copy, copytree
import os

#Custom library
from sh import *

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
	initialStep = 0
	wd = "/home/hector/SmartAS/"
	gaudiWd = "/sbi/users/hectorc/SmartAS/Results/iLoops"
	minExpression = 0
	minCandidateExpression = 4
	minPSI = 0.25

	Conditions = ["10", "7"]
	Compartments = ["C"]
	Replicates = ["1", "2"]
	Kmer = ["30"]

	top = 2

	opts, args = getopt.getopt(argv, "s:wd:cp:ce:me:")

	for opt, arg in opts:
		if opt == "-s":
			initialStep = int(arg)
		elif opt == "-wd":
			wd = arg
		elif opt == "-cp":
			minPSI = arg
		elif opt == "-ce":
			minCandidateExpression = arg
		elif opt == "-me":
			minExpression = arg

	setEnvironment(wd, initialStep, Conditions, Compartments, Replicates, Kmer)
	
	if initialStep <= 1:
		exploreData()
	if initialStep <= 2:
		getCandidates(minExpression, minCandidateExpression, minPSI)
	if initialStep <= 3:
		bianaInteractions(top)
	if initialStep <= 4:
		prepareILoopsInput()
	if initialStep <= 5:
		launchILoops()
	if initialStep <= 6:
		exloreILoopsResults()
	
	#copytree("Results", "../Dropbox/SmartAS")
	#copy("SmartAS.RData", "../Dropbox/SmartAS")

def exploreData():
	
	print "* Reading and summarizing input files: computing PSI values and plotting correlations between replicates."
	cmd("Pipeline/ExploreData.r")
	copy("SmartAS.RData", "Results/RWorkspaces/1_ExploreData.RData")

def getCandidates(minExpression, minCandidateExpression, minPSI):

	print "* Extracting transcripts with high variance and high expression."
	cmd("Pipeline/GetCandidates.r", minExpression, minCandidateExpression, minPSI)
	copy("SmartAS.RData", "Results/RWorkspaces/2_GetCandidates.RData")

	cmd("sort Results/expressedGenes.lst >Results/expressedGenes.tmp.lst")
	cmd("mv Results/expressedGenes.tmp.lst Results/expressedGenes.lst")
	
def bianaInteractions(top):

	print "* Querying BIANA for known interactions of the candidates."
	cmd("Pipeline/bianaInteractions.py", top)

def prepareILoopsInput():

	getExpressedGenes = 1

	print "* Retrieving protein sequences for transcripts and printing to multiFASTA file."
	
	diff = Popen('diff old/expressedGenes.lst Results/expressedGenes.lst 2>&1', shell=True, stdout=PIPE)

	if not diff.stdout.read().strip() and os.path.exists("old/iLoops/input/ExpressedTranscripts.fasta"):
		getExpressedGenes = 0

	cmd("Pipeline/getiLoopsInput.pl Results/expressedGenes.lst Results/Results/candidateList.top.lst ", getExpressedGenes)

def launchILoops():

	print "* Launching iLoops jobs."
	
	cmd("scp -r Results/iLoops hectorc@gaudi.imim.es:~/SmartAS/Results")
	cmd("ssh hectorc@gaudi '~/SmartAS/Pipeline/launchILoops.sh /sbi/users/hectorc/SmartAS/Results/iLoops'")

	print "\t* Waiting..."

def exloreILoopsResults():

	print "* Examining iLoops results."
	cmd("Pipeline/exploreiLoopsOutput.py")

main(sys.argv[1:])
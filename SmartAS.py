#!/soft/devel/python-2.7/bin/python

#from subprocess import Popen,PIPE
import sys, getopt
from shutil import copy, copytree
from os import path

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
		elif opt == "-w":
			wd = arg
		elif opt == "-p":
			minPSI = arg
		elif opt == "-e":
			minCandidateExpression = arg
		elif opt == "-m":
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
	
	print("* Reading and summarizing input files: computing PSI values and plotting correlations between replicates.")
	cmd("Pipeline/ExploreData.r")
	copy("SmartAS.RData", "Results/RWorkspaces/1_ExploreData.RData")

def getCandidates(minExpression, minCandidateExpression, minPSI):

	print("* Extracting transcripts with high variance and high expression.")
	cmd("Pipeline/GetCandidates.r", minExpression, minCandidateExpression, minPSI)
	cmd("Pipeline/getGenenames.py", "Results/candidateList.lst")
	copy("SmartAS.RData", "Results/RWorkspaces/2_GetCandidates.RData")

	cmd("sort Results/expressedGenes.lst >Results/expressedGenes.tmp.lst")
	cmd("mv Results/expressedGenes.tmp.lst Results/expressedGenes.lst")
	
def bianaInteractions(top):

	print("* Querying BIANA for known interactions of the candidates.")
	cmd("Pipeline/bianaInteractions.py", top)

def prepareILoopsInput():

	getExpressedGenes = 1

	print("* Retrieving protein sequences for transcripts and printing to multiFASTA file.")
	
	diff = cmdOut('diff old/expressedGenes.lst Results/expressedGenes.lst 2>&1')
	
	if not diff and path.exists("old/iLoops/ExpressedTranscripts.fasta"):
		getExpressedGenes = 0

	cmd("Pipeline/getiLoopsInput.py Results/expressedGenes.lst Results/candidateList.top.lst", getExpressedGenes)

def launchILoops():

	print("* Launching iLoops jobs.")
	
	cmd("ssh hectorc@gaudi 'mv ~/SmartAS/Results ~/SmartAS/old; rm -r ~/SmartAS/old'")
	cmd("scp -r Results/iLoops hectorc@gaudi.imim.es:~/SmartAS/Results")
	cmd("ssh hectorc@gaudi '~/SmartAS/Pipeline/launchILoops.py /sbi/users/hectorc/SmartAS/Results/iLoops'")

	print("\t* Waiting...")

def exloreILoopsResults():

	print("* Examining iLoops results.")
	cmd("Pipeline/exploreiLoopsOutput.py")

main(sys.argv[1:])

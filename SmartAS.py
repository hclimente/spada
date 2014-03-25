#!/soft/devel/python-2.7/bin/python

import sys, getopt
from os import path, chdir
from include.libsmartas import cmd, setEnvironment, finish, outputCandidates

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
		if not opt["external"]:
			exit()
	if opt["initialStep"] <= 4:
		launchiLoops(opt)
	if opt["initialStep"] <= 5:
		analyzeInteractions(opt)
	
	#finish(opt)

def exploreData(opt):
	
	print("* Reading and summarizing input files: computing PSI values and intereplicate agreement.")
	cmd("Pipeline/ExploreData.r", opt["out"], "Data/Input/" + opt["inputType"] + "/" + opt["tag1"] + "/")

def getCandidates(opt):

	print("* Extracting transcripts with high variance and high expression.")
	cmd("Pipeline/GetCandidates.r", opt["minExpression"], opt["out"])
	
	cmd("sort", "Results/" + opt["out"] + "/expressedGenes.lst", ">" + "Results/" + opt["out"] + "/expressedGenes.tmp.lst")
	cmd("mv", "Results/" + opt["out"] + "/expressedGenes.tmp.lst", "Results/" + opt["out"] + "/expressedGenes.lst")

	outputCandidates(opt["out"], opt["inputType"])
	
def candidatePrioritization(opt):

	print("* Prioritizing candidates.")
	cmd("Pipeline/CandidatePrioritization.py", opt["out"], opt["inputType"])

def launchiLoops(opt):

	print("* Sending files to Gaudi and performing the iLoops analysis.")
	cmd("ssh hectorc@gaudi 'rm -r", opt["gOut"] + "'")
	cmd("ssh hectorc@gaudi 'mkdir -p", opt["gOut"] + "/iLoops/Output; mkdir -p", opt["gOut"] + "/iLoops/Input; mkdir -p", opt["gOut"] + "/iLoops/logs'")
	cmd("scp -r " + "Results/" + opt["out"] + "/expressedGenes.lst Results/" + opt["out"] + "/candidateList.top.tsv hectorc@gaudi.imim.es:" + opt["gOut"])

	cmd("ssh hectorc@gaudi '" + opt["gaudiWd"] + "/Pipeline/CalculateInteractions.py", opt["gaudiWd"], opt["out"], opt["inputType"] + "'")

def analyzeInteractions(opt):

	print("* Examining iLoops results.")
	cmd("Pipeline/AnalyzeInteractions.py", opt["out"])

main(sys.argv[1:])
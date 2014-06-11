#!/soft/devel/python-2.7/bin/python

import sys, getopt
from os import path, chdir
from libs.utils import cmd, cmdOut, setEnvironment, finish, outputCandidates, pickUniqPatterns
from libs import options

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
    options.Options()
    options.Options().printToFile()

	cfgFile = ""

	for opt, arg in opts:
		if opt == "-f":
			cfgFile = arg
		else:
			print("No configuration file found.")
			exit()
		
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
	cmd("Pipeline/methods/ExploreData.r", opt["out"], "Data/Input/" + opt["inputType"] + "/" + opt["tag1"] + "/")

def getCandidates(opt):

	print("* Extracting transcripts with high variance and high expression.")
	cmd("Pipeline/methods/GetCandidates.r", opt["minExpression"], opt["out"], opt["unpairedReplicates"])
	
	cmd("sort", "Results/" + opt["out"] + "/expressedGenes.lst", ">" + "Results/" + opt["out"] + "/expressedGenes.tmp.lst")
	cmd("mv", "Results/" + opt["out"] + "/expressedGenes.tmp.lst", "Results/" + opt["out"] + "/expressedGenes.lst")

	outputCandidates(opt["out"], opt["inputType"])
	
def candidatePrioritization(opt):

	print("* Prioritizing candidates.")
	cmd("Pipeline/methods/CandidatePrioritization.py", opt["out"], opt["inputType"])

def launchiLoops(opt):

	print("* Selecting isoforms suitable for " + opt["iLoopsVersion"])
	pickUniqPatterns(opt["gOut"], opt["out"], opt["inputType"], opt["iLoopsVersion"], opt["Replicates"] * 0.1)

	print("* Sending list to Gaudi and performing the iLoops analysis.")
	gaudiThread = cmdOut(
							"ssh", "hectorc@gaudi", \
							"'" + opt["gaudiWd"] + "/Pipeline/methods/CalculateInteractions.py " + opt["gaudiWd"] + " " + opt["out"] + " " + opt["inputType"] + " " + opt["iLoopsVersion"] + "'", \
							">Results/" + opt["out"] + "/calculateInteractions.log"
						)

def analyzeInteractions(opt):

	print("* Examining iLoops results.")
	cmd("Pipeline/methods/AnalyzeInteractions.py", opt["out"], opt["inputType"], opt["iLoopsVersion"], opt["Replicates"] * 0.1)

main(sys.argv[1:])
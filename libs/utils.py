#!/soft/devel/python-2.7/bin/python

from subprocess import call,Popen,PIPE
from rpy2.robjects import r
from time import sleep
import urllib2
from os import chdir
from os.path import isfile
import logging

from libs import options

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
			if header and counter is 0:
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
	cmd("mkdir", o.qout + o.iLoopsVersion)

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

			#currentInitialState = o.initialStep
			#opt = parseParam("old/" + o.out + "/Parameters.cfg")
			#o.initialStep = currentInitialState
			#cmd("cp", "old/" + o.out + "/Parameters.cfg", "Results/" + o.out)
			cmd("cp -r", "old/" + o.out + "/DataExploration", o.qout)
			cmd("cp", "old/" + o.out + "/RWorkspaces/1_ExploreData.RData", o.qout + "/RWorkspaces")

		if o.initialStep > 2:
			cmd("cp", "old/" + o.out + "/RWorkspaces/2_GetCandidates.RData", o.qout + "/RWorkspaces")
			cmd("cp", "old/" + o.out + "/candidateList.tsv", "old/" + o.out + "/expressedGenes.lst", o.qout)
			cmd("cp", "old/" + o.out + "/candidates_normal.gtf", "old/" + o.out + "/candidates_tumor.gtf", o.qout)
	else:
		
		if o.initialStep == 2:
			oriOut = "_".join(o.out.split("_")[:-1])
			cmd("cp", "Results/" + oriOut + "/RWorkspaces/1_ExploreData.RData", "Results/" + o.out + "/RWorkspaces")
			cmd("Pipeline/scripts/InputUnpaired.r", oriOut, o.unpairedReplicates, "Data/Input/TCGA/" + o.tag + "/")
		if o.initialStep == 3:
			cmd("cp", o.external + ".tsv" , "Results/" + o.out + "/candidateList.tsv")
			cmd("cp", o.external + "_expressedGenes.lst", "Results/" + o.out + "/expressedGenes.lst")

	if o.initialStep > 3:
		cmd("cp", "old/" + o.out + "/candidateInteractions.tsv", "old/" + o.out + "/candidateList.top.tsv", o.qout)
	if o.initialStep > 4:
		cmd("cp", "old/" + o.out + "/candidatesGaudi.lst", o.qout)
	if o.initialStep > 5:
		cmd("cp", "old/" + o.out + "/iLoops", o.qout)

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

# def printParam(opt):
# 	with open("Results/" + opt["out"] + "/Parameters.cfg", "w") as paramFile:
# 		for aKey in opt.keys():
# 			if not aKey == "Conditions":
# 				paramFile.write(aKey + "=" + str(opt[aKey]) + "\n")

# def parseParam(cfgFile):

# 	opt = { "initialStep" : 0, "wd" : "/home/hector/SmartAS/", "gaudiWd" : "/sbi/users/hectorc/SmartAS", "minExpression" : 0, 
# 			"inputType" : "GENCODE" , "Conditions" : ["N", "T"], "tag1" : "20", "Replicates" : 0, "external" : "", 
# 			"unpairedReplicates" : 0, "iLoopsVersion" : "iLoops_devel"}

# 	with open(cfgFile, "r") as PARAMETERS:
# 		for line in PARAMETERS:
# 			elements = line.strip().split("=")

# 			if elements[0] in ["wd", "gaudiWd", "inputType", "tag1", "out", "gOut", "external", "iLoopsVersion"]:
# 				opt[elements[0]] = elements[1]
# 			elif elements[0] in ["Replicates", "initialStep", "unpairedReplicates"]:
# 				opt[elements[0]] = int(elements[1])
# 			elif elements[0] in ["minExpression"]:
# 				opt[elements[0]] = float(elements[1])
# 			else:
# 				print("Unrecognized option:" + line.strip() )

# 	if opt["external"]:
# 		if opt["unpairedReplicates"] == 0:
# 			opt["out"] = opt["inputType"] + "/" + opt["tag1"]
# 	else:
# 		opt["out"] = opt["inputType"] + "/" + opt["tag1"] + "_mE" + str(opt["minExpression"])

# 	opt["gOut"] = opt["gaudiWd"] + "/Results/" + opt["out"]

# 	return opt

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

	chdir(outFolder)
	cmd("cp -r", "../SmartAS/Results/" + opt["out"] + "/* .")

	cmd("tar -czvf", outTag + ".tar.gz candidates_normal." + outTag + ".gtf candidates_tumor." + outTag + ".gtf candidateList." + outTag + ".tsv")
	cmd("rm", "*" + opt["out"] + "*gff", "*" + opt["out"] + "*tsv")
	chdir("/home/hector/SmartAS")

# def outputCandidates(out, inputType):
# 	with open('Results/' + out + "/candidateList.tsv", "r") as CANDIDATES, \
# 		 open('Results/' + out + "/candidates_normal.gtf", 'w') as GTFn, \
# 		 open("Results/" + out + "/candidates_tumor.gtf", 'w') as GTFt, \
# 		 open("Data/" + inputType + "/annotation.gtf", "r") as ALLTRANSCRIPTS:
# 			candTnt = []
# 			for line in CANDIDATES:
# 				ids = line.strip().split("\t")
# 				candTnt.append([ids[2], ids[3]])
			
# 			for line in ALLTRANSCRIPTS:
# 				for pair in candTnt:
# 					if line.find(pair[0]) != -1:
# 						GTFn.write(line)
# 					elif line.find(pair[1]) != -1:
# 						GTFt.write(line)

def pickUniqPatterns(gOut, out, inputType, iLoopsVersion, minReplicates):

	loopFamilies = {}

	with open("Data/TCGA/UnifiedFasta_" + iLoopsVersion + "_loopFamilies.txt", "r") as LOOP_FAMILIES:
		currentLoopFamily = ""
		for line in LOOP_FAMILIES:
			if ">" in line:
				currentLoopFamily = line[1:].strip().split("\t")[0]
				loopFamilies[currentLoopFamily] = set()
				loopFamilies[currentLoopFamily].add(line.strip().split("\t")[1])
			else:
				loopFamilies[currentLoopFamily].add(line.strip())

	with open("Results/" + out + "/candidateList.top.tsv", "r") as CANDIDATES, \
		 open("Results/" + out + "/candidatesGaudi.lst", "w") as CANDIDATES_GAUDI:

		CANDIDATES.readline()
		analyzedLoops = {}

		for line in CANDIDATES:
			elements = line.split("\t")
			isoN = elements[2]
			isoT = elements[3]
			reps = int(elements[4])
			analyze = 0
			comment = "To analyze."

			if reps < minReplicates:
				break

			with open("Data/" + inputType + "/UnifiedFasta_" + iLoopsVersion + ".fa", "r") as MULTIFASTA:
				loopsN = ""
				loopsT = ""
				for line in MULTIFASTA:
					if isoN in line:
						loopsN = line.strip().split("#")[3]
					elif isoT in line:
						loopsT = line.strip().split("#")[3]

			if elements[9] == "No":
				analyze = -1
				comment = "No CDS change."
			elif not loopsN or not loopsT:
				analyze = -1
				comment = "No loops mapped with " + iLoopsVersion + "."
			elif loopsN == loopsT:
				analyze = -1
				comment = "No difference in loops mapped with " + iLoopsVersion + "."

			if analyze < 0:
				CANDIDATES_GAUDI.write( isoN + "\t" + str(analyze) + "\t" + comment + "\n")
				CANDIDATES_GAUDI.write( isoT + "\t" + str(analyze) + "\t" + comment + "\n")
				continue

			for candidate, thisLoopPattern in zip([isoN, isoT], [loopsN, loopsT]):
		
				analyze = 0
				comment = "To analyze"

				if isfile("iLoops/TCGA/" + iLoopsVersion + "/" + candidate + ".tar.gz"):
					analyze = 1
					comment = "Already analyzed."
				elif thisLoopPattern in analyzedLoops.keys():
					analyze = 2
					comment = "Analyzing relative " + analyzedLoops[thisLoopPattern] + "."
				else:
					for isoform in loopFamilies[thisLoopPattern]:
						if isfile("iLoops/TCGA/" + iLoopsVersion + "/" + isoform + ".tar.gz"):
							analyze = 2
							comment = "Analyzed relative " + isoform + "."
							break

				CANDIDATES_GAUDI.write( candidate + "\t" + str(analyze) + "\t" + comment + "\n")

				if analyze == 0:
					analyzedLoops[thisLoopPattern] = candidate

	cmd("ssh hectorc@gaudi 'rm -r", gOut + "'")
	cmd("ssh hectorc@gaudi 'mkdir -p", gOut + "/Output; mkdir -p", gOut + "/Input; mkdir -p", gOut + "/logs'")
	cmd("scp -r " + "Results/" + out + "/candidatesGaudi.lst hectorc@gaudi.imim.es:" + gOut)

# def importSequences():
# 	kkota = {}
# 	with open("/home/hector/SmartAS/Data/TCGA/UnifiedFasta_iLoops13.fa") as KK:
# 		currentTx = ""
# 		seq = ""
# 		for line in KK:
# 			if ">" in line:
# 				if currentTx:
# 					kkota[currentTx]["Sequence"] = seq
# 					seq=""
# 				elements=line[1:].strip().split("#")
# 				currentTx = elements[0]
# 				kkota.setdefault(currentTx, {})
# 				kkota[currentTx]["Gene"] = elements[1]
# 				kkota[currentTx]["Uniprot"] = elements[2]
# 				kkota[currentTx]["iLoopsFamily"] = elements[1]
# 			else:
# 				seq += line.strip()

# 	return kkota

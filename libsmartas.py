#!/soft/devel/python-2.7/bin/python

from subprocess import call,Popen,PIPE
from rpy2.robjects import r
from time import sleep
import urllib2
from os import chdir, listdir
from fnmatch import filter

def cmd(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)

	call(command, shell=True)

def cmdOut(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)
	
	return Popen(command, shell=True, stdout=PIPE).stdout.read().strip()

def setRWorkspace(wd, out, Replicates):

	r("wd <- \"" + wd + "\"")
	r("out <- \"Results/" + out + "/\"")
	r("inputData <- list()")	
	r('inputData[["Conditions"]] <- c("N", "T")')
	r('inputData[["Replicates"]] <- ' + Replicates)
	r('save.image("SmartAS.RData")')

def setEnvironment(cfgFile):

	opt = parseParam(cfgFile)

	opt["out"] = opt["inputType"] + "/mE" + str(opt["minExpression"]) + "_mCE" + str(opt["minCandidateExpression"]) + "_mP" + str(opt["minPSI"])
	opt["gOut"] = opt["gaudiWd"] + "/Results/" + opt["out"]

	print("* Preparing the environment")
	cmd("rm -r old2/" + opt["out"] + "; mv", "old/" + opt["out"], "old2/"  + opt["out"])
	cmd("mv", "Results/" + opt["out"], "old/" + opt["out"] + "; mv SmartAS.RData ", "old/" + opt["out"] + "/SmartAS.old.RData")
	cmd("mkdir -p", "Results/" + opt["out"] + "/iLoops/Output/Mapping")
	cmd("mkdir", "Results/" + opt["out"] + "/iLoops/Input")
	cmd("mkdir", "Results/" + opt["out"] + "/RWorkspaces")
	cmd("mkdir", "Results/" + opt["out"] + "/DataExploration")
	cmd("mkdir", "Results/" + opt["out"] + "/Input")

	if opt["initialStep"] <= 1:
		opt["Replicates"] = standarizeInput(opt)
		printParam(opt)
		setRWorkspace(opt["wd"], opt["out"], opt["Replicates"])
		getDB()

	if opt["initialStep"] > 1:

		print opt
		currentInitialState = opt["initialStep"]
		opt = parseParam("old/" + opt["out"] + "/Parameters.cfg")
		print opt
		opt["initialStep"] = currentInitialState

		cmd("cp -r", "old/" + opt["out"] + "/DataExploration", "Results/" + opt["out"])
		cmd("cp -r", "old/" + opt["out"] + "/RWorkspaces/1_ExploreData.RData", "Results/" + opt["out"] + "/RWorkspaces")
		cmd("cp -r", "old/" + opt["out"] + "/RWorkspaces/1_ExploreData.RData SmartAS.RData")

		for rep in range(1, int(opt["Replicates"]) + 1):
			for cond in opt["Conditions"]:
				cmd("cp -r", "old/" + opt["out"] + "/" + str(rep) + "_" + cond + ".tsv", "Results/" + opt["out"])
				cmd("cp -r", "old/" + opt["out"] + "/IntraReplicate_" + str(rep) + ".tsv", "Results/" + opt["out"])

	if opt["initialStep"] > 2:
		cmd("cp -r", "old/" + opt["out"] + "/RWorkspaces/2_GetCandidates.RData", "Results/" + opt["out"] + "/RWorkspaces")
		cmd("cp -r", "old/" + opt["out"] + "/RWorkspaces/2_GetCandidates.RData SmartAS.RData")
		cmd("cp", "old/" + opt["out"] + "/candidateList.tsv", "old/" + opt["out"] + "/expressedGenes.lst", "Results/" + opt["out"])
		cmd("cp", "old/" + opt["out"] + "/candidates_normal.gff", "old/" + opt["out"] + "/candidates_tumor.gff", "Results/" + opt["out"])
	if opt["initialStep"] > 3:
		cmd("cp", "old/" + opt["out"] + "/candidateInteractions.tsv", "old/" + opt["out"] + "/candidateList.top.tsv", "Results/" + opt["out"])
	if opt["initialStep"] > 4:
		cmd("cp -r", "old/" + opt["out"] + "/iLoops/Input", "Results/" + opt["out"] + "/iLoops/")
		cmd("cp -r", "old/" + opt["out"] + "/iLoops/ExpressedTranscripts.fasta", "Results/" + opt["out"] + "/iLoops/")
		cmd("cp", "old/" + opt["out"] + "/candidates_normal.top.gff", "old/" + opt["out"] + "/candidates_tumor.top.gff", "Results/" + opt["out"])
	if opt["initialStep"] > 5:
		cmd("cp -r", "old/" + opt["out"] + "/iLoops/Output", "Results/" + opt["out"] + "/iLoops")

	return opt

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

def waitPID(pidQueue):
	for job in pidQueue:
		while True:
			if cmdOut("ps --pid", job, " | grep -v", job, "| wc -l") == "0":
				break
			else:
				print("Awaiting for completion of iLoops jobs.")
				sleep(900)

def printParam(opt):
	with open("Results/" + opt["out"] + "/Parameters.cfg", "w") as paramFile:
		for aKey in opt.keys():
			if not aKey == "Conditions":
				paramFile.write(aKey + "=" + str(opt[aKey]) + "\n")

def parseParam(cfgFile):

	opt = { "initialStep" : 0, "wd" : "/home/hector/SmartAS/", "gaudiWd" : "/sbi/users/hectorc/SmartAS",
		    "minExpression" : 0, "minCandidateExpression" : 4, "minPSI" : 0.25, "inputType" : "GENCODE" , "Conditions" : ["N", "T"],
		    "compartment" : "C", "kmer" : "20", "Replicates" : 0
	}

	with open(cfgFile, "r") as PARAMETERS:
		for line in PARAMETERS:
			elements = line.strip().split("=")
			if elements[0] == "Replicates":
				Replicates = int(elements[1])
			elif elements[0] == "initialStep":
				opt["initialStep"] = int(elements[1])
			elif elements[0] == "wd":
				opt["wd"] = elements[1]
			elif elements[0] == "gaudiWd":
				opt["gaudiWd"] = elements[1]
			elif elements[0] == "minExpression":
				opt["minExpression"] = float(elements[1])
			elif elements[0] == "minCandidateExpression":
				opt["minCandidateExpression"] = float(elements[1])
			elif elements[0] == "minPSI":
				opt["minPSI"] = float(elements[1])
			elif elements[0] == "inputType":
				opt["inputType"] = elements[1]
			elif elements[0] == "kmer":
				opt["kmer"] = elements[1]
			elif elements[0] == "out":
				opt["out"] = elements[1]
			elif elements[0] == "gOut":
				opt["gOut"] = elements[1]
			elif elements[0] == "compartment":
				opt["compartment"] = elements[1]
			else:
				print("Unrecognized option:" + line)

	return opt

def standarizeInput(opt):

	reps = []

	if(opt["inputType"] == "GENCODE"):
		Conditions = {"10": "N", "7": "T"}
		for condition in Conditions.keys():
			replicateCounter = 1
			for sample in filter(listdir("Data/GENCODE/Rawdata/" + opt["kmer"] + "-kmer-length/"), condition + "C*_*"):
				with open("Data/GENCODE/Rawdata/" + opt["kmer"] + "-kmer-length/" + sample + "/quant_bias_corrected.sf", "r") as FILE, \
					 open("Results/" + opt["out"] + "/Input/" + str(replicateCounter) + "_" + Conditions[condition] + ".tsv", "w") as FILTERED:
					for line in FILE:
						if line.find("#") == -1:
							tableValues=line.strip().split("\t")
							splitIds=tableValues[0].split("|")
							FILTERED.write(splitIds[1].split(".")[0] + "\t" + splitIds[0].split(".")[0] + "\t" + splitIds[5] + "\t" + tableValues[2] + "\n")
				reps.append(str(replicateCounter))
				replicateCounter += 1
	
	elif(opt["inputType"] == "TCGA"):
		patients = []
		with open("Data/TCGA/Rawdata/ucec_iso_tpm_paired.txt", "r") as FILE:
			firstLine = FILE.readline().strip().split("\t")
			for sampleType in opt["Conditions"]:
				replicateCounter = 1
				for patient in range(0, len(firstLine)/2):
					patients.append(open("Results/" + opt["out"] + "/Input/" + str(replicateCounter) + "_" + sampleType + ".tsv", "w"))
					reps.append(str(replicateCounter))
					replicateCounter += 1
			
			for line in FILE:
							
				splitted = line.strip().split("\t")
				seqTags = splitted[0].split(",")
				
				gene = seqTags[0]
				transcript = seqTags[1]
				name = gene.split("|")[0]
				
				currentCol = 1
				for patient in patients:
					patient.write( gene + "\t" + transcript + "\t" + name + "\t" + splitted[currentCol] + "\n")
					currentCol += 1
		for patient in patients:
			patient.close()
	
	return str(max(reps))

def finish(opt):
	
	print("* Moving files to the Results directory and creating a summary tar file.")
	minExpresion = ""
	minPSI = ""

	with open("Results/" + opt["out"] + "/Parameters.cfg", "r") as paramFile:
		for line in paramFile:
			elements = line.split("=")

			if elements[0] == "minCandidateExpression":
				minExpresion = elements[1].strip()
			elif elements[0] == "minPSI":
				minPSI = elements[1].strip()
	    
	outFolder = "/home/hector/Results/" + opt["inputType"] + "/" + opt["out"]

	if not cmdOut("mkdir -p", outFolder):
		overwrite = raw_input("\tDirectory exists. Do you want to overwrite it? (y/n)")
		if not overwrite == "y":
			return

	cmd("cp", "Results/" + opt["out"] + "/candidates_normal.gff", outFolder + "/" + "candidates_normal." + opt["out"] + ".gff")
	cmd("cp", "Results/" + opt["out"] + "/candidates_tumor.gff", outFolder + "/" + "candidates_tumor." + opt["out"] + ".gff")
	cmd("cp", "Results/" + opt["out"] + "/candidateList.top.tsv", outFolder + "/" + "candidateList." + opt["out"] + ".tsv")

	chdir(outFolder)
	cmd("cp -r", "../SmartAS/Results/" + opt["out"] + "/* .")

	cmd("tar -czvf", opt["out"] + ".tar.gz candidates_normal." + opt["out"] + ".gff candidates_tumor." + opt["out"] + ".gff candidateList." + opt["out"] + ".tsv")
	cmd("rm", "*" + opt["out"] + "*gff", "*" + opt["out"] + "*tsv")
	chdir("/home/hector/SmartAS")
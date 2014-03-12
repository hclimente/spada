#!/soft/devel/python-2.7/bin/python

from subprocess import call,Popen,PIPE
from rpy2.robjects import r
from time import sleep
import urllib2
from os import chdir

def cmd(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)

	call(command, shell=True)

def cmdOut(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)
	
	return Popen(command, shell=True, stdout=PIPE)

def setRWorkspace(wd, out, Replicates):

	r("wd <- \"" + wd + "\"")
	r("out <- \"Results/" + out + "/\"")
	r("inputData <- list()")	
	r('inputData[["Conditions"]] <- c("N", "T")')
	r('inputData[["Replicates"]] <- ' + str(Replicates))
	r('save.image("Results/' + out + '/RWorkspaces/0_InitialEnvironment.RData")')

def setEnvironment(cfgFile):

	opt = parseParam(cfgFile)

	print("* Preparing the environment")
	cmd("rm -r old2/" + opt["out"] + "; mv", "old/" + opt["out"], "old2/"  + opt["out"])
	cmd("mv", "Results/" + opt["out"], "old/" + opt["out"])
	cmd("mkdir -p", "Results/" + opt["out"] + "/RWorkspaces")
	cmd("mkdir", "Results/" + opt["out"] + "/DataExploration")

	if opt["initialStep"] <= 1 and not opt["External"]:
		printParam(opt)
		setRWorkspace(opt["wd"], opt["out"], opt["Replicates"])
		getDB()

	if opt["initialStep"] > 1 and not opt["External"]:

		currentInitialState = opt["initialStep"]
		opt = parseParam("old/" + opt["out"] + "/Parameters.cfg")
		opt["initialStep"] = currentInitialState

		cmd("cp", "old/" + opt["out"] + "/Parameters.cfg", "Results/" + opt["out"])
		cmd("cp -r", "old/" + opt["out"] + "/DataExploration", "Results/" + opt["out"])
		cmd("cp -r", "old/" + opt["out"] + "/RWorkspaces/1_ExploreData.RData", "Results/" + opt["out"] + "/RWorkspaces")

	if opt["initialStep"] > 2 and not opt["External"]:
		cmd("cp -r", "old/" + opt["out"] + "/RWorkspaces/2_GetCandidates.RData", "Results/" + opt["out"] + "/RWorkspaces")
		cmd("cp", "old/" + opt["out"] + "/candidateList.tsv", "old/" + opt["out"] + "/expressedGenes.lst", "Results/" + opt["out"])
		cmd("cp", "old/" + opt["out"] + "/candidates_normal.gtf", "old/" + opt["out"] + "/candidates_tumor.gtf", "Results/" + opt["out"])
	if opt["initialStep"] > 3:
		cmd("cp", "old/" + opt["out"] + "/candidateInteractions.tsv", "old/" + opt["out"] + "/candidateList.top.tsv", "Results/" + opt["out"])
	if opt["initialStep"] > 4:
		cmd("mv", "old/" + opt["out"] + "/iLoops", "Results/" + opt["out"])
	else:
		cmd("mkdir -p", "Results/" + opt["out"] + "/iLoops/Output/Mapping")
		cmd("mkdir", "Results/" + opt["out"] + "/iLoops/Input")
	if opt["initialStep"] > 5:
		pass
		#cmd("cp -r", "old/" + opt["out"] + "/iLoops/Output", "Results/" + opt["out"] + "/iLoops")

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

def printParam(opt):
	with open("Results/" + opt["out"] + "/Parameters.cfg", "w") as paramFile:
		for aKey in opt.keys():
			if not aKey == "Conditions":
				paramFile.write(aKey + "=" + str(opt[aKey]) + "\n")

def parseParam(cfgFile):

	opt = { "initialStep" : 0, "wd" : "/home/hector/SmartAS/", "gaudiWd" : "/sbi/users/hectorc/SmartAS", "minExpression" : 0, 
			"inputType" : "GENCODE" , "Conditions" : ["N", "T"], "tag1" : "20", "Replicates" : 0, "External" : False }

	with open(cfgFile, "r") as PARAMETERS:
		for line in PARAMETERS:
			elements = line.strip().split("=")
			if elements[0] == "Replicates":
				opt["Replicates"] = int(elements[1])
			elif elements[0] == "initialStep":
				opt["initialStep"] = int(elements[1])
			elif elements[0] == "wd":
				opt["wd"] = elements[1]
			elif elements[0] == "gaudiWd":
				opt["gaudiWd"] = elements[1]
			elif elements[0] == "minExpression":
				opt["minExpression"] = float(elements[1])
			elif elements[0] == "inputType":
				opt["inputType"] = elements[1]
			elif elements[0] == "tag1":
				opt["tag1"] = elements[1]
			elif elements[0] == "out":
				opt["out"] = elements[1]
			elif elements[0] == "gOut":
				opt["gOut"] = elements[1]
			elif elements[0] == "External":
				if elements[1] == "True":
					opt["External"] = True
				elif elements[1] == "False":
					opt["External"] = False
			else:
				print("Unrecognized option:" + line.strip() )

	opt["out"] = opt["inputType"] + "/" + opt["tag1"] + "_mE" + str(opt["minExpression"])
	if opt["External"]:
		opt["out"] = opt["inputType"] + "/" + opt["tag1"]

	opt["gOut"] = opt["gaudiWd"] + "/Results/" + opt["out"]

	return opt

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

def outputCandidates(out, inputType):
	with open('Results/' + out + "/candidateList.tsv", "r") as CANDIDATES, \
		 open('Results/' + out + "/candidates_normal.gtf", 'w') as GTFn, \
		 open("Results/" + out + "/candidates_tumor.gtf", 'w') as GTFt, \
		 open("Data/" + inputType + "/annotation.gtf", "r") as ALLTRANSCRIPTS:
			candTnt = []
			for line in CANDIDATES:
				ids = line.strip().split("\t")
				candTnt.append([ids[2], ids[3]])
			
			for line in ALLTRANSCRIPTS:
				for pair in candTnt:
					if line.find(pair[0]) != -1:
						GTFn.write(line)
					elif line.find(pair[1]) != -1:
						GTFt.write(line)
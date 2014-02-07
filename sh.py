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
	
	return Popen(command, shell=True, stdout=PIPE).stdout.read().strip()

def setRWorkspace(wd, Conditions, Compartments, Replicates, Kmer):

	r("wd <- \"" + wd + "\"")
	r("inputData <- list()")
	
	condStat = 'inputData[["Conditions"]] <- c('
	for condition in Conditions:
		condStat = condStat + "\"" + condition + "\","
	condStat = condStat[:-1] + ")"
	
	compStat = 'inputData[["Compartments"]] <- c('
	for compartment in Compartments:
		compStat = compStat + "\"" + compartment + "\","
	compStat = compStat[:-1] + ")"
	
	replStat = 'inputData[["Replicates"]] <- c('
	for replicate in Replicates:
		replStat = replStat + "\"" + replicate + "\","
	replStat = replStat[:-1] + ")"
	
	kmerStat = 'inputData[["K-mer"]] <- c('
	for kmer in Kmer:
		kmerStat = kmerStat + "\"" + kmer + "\","
	kmerStat = kmerStat[:-1] + ")"
	
	r(condStat)
	r(compStat)
	r(replStat)
	r(kmerStat)
	
	r('save.image("SmartAS.RData")')

def setEnvironment(wd, initialStep, Conditions, Compartments, Replicates, Kmer, inputType):

	print("* Preparing the environment")
	cmd("cd " + wd)
	cmd("rm -r old2; mv old old2")
	cmd("mv Results old; mv SmartAS.RData old/SmartAS.old.RData")
	cmd("mkdir -p Results/iLoops/Output/Mapping")
	cmd("mkdir Results/iLoops/Input")
	cmd("mkdir Results/RWorkspaces")
	cmd("mkdir Results/DataExploration")

	if initialStep <= 1:
		setRWorkspace(wd, Conditions, Compartments, Replicates, Kmer)
		getDB()

		if(inputType == "GENCODE"):
			for kmer in Kmer:
  				for replicate in Replicates:
  					for condition in Conditions:
						tag = condition + "C" + replicate + "_" + kmer
						with open("Data/GENCODE/Rawdata/" + kmer + "-kmer-length/" + tag + "/quant_bias_corrected.sf", "r") as FILE, \
							 open("Data/GENCODE/Filtered/" + tag + ".filtered.sf", "w") as FILTERED:
							for line in FILE:
								if line.find("#") == -1:
									tableValues=line.split("\t")
									splitIds=tableValues[0].split("|")
									FILTERED.write(splitIds[1].split(".")[0] + "\t" + splitIds[0].split(".")[0] + "\t" + splitIds[5] + "\t" + tableValues[2] + "\n")
		elif(inputType == "TCGA"):
			patients = []

			with open("Data/TCGA/Rawdata/ucec_iso_tpm_paired.txt", "r") as FILE:
				firstLine = FILE.readline().split("\t")
				for patient in firstLine:
					patients.append(patient[0:len(patient) - 1])
				
				with open("Data/TCGA/Filtered/" + tag + ".filtered.sf", "w") as FILTERED:
				FILTERED.write( gene + "\t" + transcript + "\t" + name + "\t" + splitted[1] + "\n")

				for line in FILE:
					transcriptInfo = {}
					
					splitted = line.split("\t")
					seqTags = splitted[1].split(",")
					
					transcriptInfo["Gene"] = seqTags[0]
					transcriptInfo["Transcript"] = seqTags[1]
					transcriptInfo["Name"] = gene.split("|")[0]
					transcriptInfo["TPMs"] = {}
					
					currentCol = 1

					for oneCondition in ["N", "T"]:
						transcriptInfo["TPMs"][oneCondition] = {}
						for patient in patients:
							transcriptInfo["TPMs"][oneCondition][str(i)] = splitted[i]
							currentCol += 1

			for oneCondition in ["N", "T"]:
				for patient in patients:										
					tag = patient + "_" + oneCondition

			with open("Data/TCGA/Filtered/" + tag + ".filtered.sf", "w") as FILTERED:
				FILTERED.write( gene + "\t" + transcript + "\t" + name + "\t" + splitted[1] + "\n")

	if initialStep > 1:
		cmd("cp -r old/DataExploration Results")
		cmd("cp -r old/RWorkspaces/1_ExploreData.RData Results/RWorkspaces")
		cmd("cp -r old/RWorkspaces/1_ExploreData.RData SmartAS.RData")
		for comp in Compartments:
			for rep in Replicates:
				for kmer in Kmer:
					for cond in Conditions:
						cmd("cp -r old/" + cond + comp + rep + "_" + kmer + ".tsv Results")
					cmd("cp -r old/IntraReplicate" + comp + rep + "_" + kmer + ".tsv Results")
		
		#If any kind of data is recycled, check that the parameters didn't change between runs.
		diff = cmdOut('diff Results/Parameters.cfg old/Parameters.cfg 2>&1')
	
		if diff:
			print("WARNING: parameter files don't match.")
			with open("Results/Parameters.cfg", "w") as paramFile:
				paramFile.write("\n***Did not match previous run***")
	if initialStep > 2:
		cmd("cp -r old/RWorkspaces/2_GetCandidates.RData Results/RWorkspaces")
		cmd("cp -r old/RWorkspaces/2_GetCandidates.RData SmartAS.RData")
		cmd("cp old/candidateList.tsv old/expressedGenes.lst Results")
		cmd("cp  old/candidates_normal.gff old/candidates_tumor.gff Results")
	if initialStep > 3:
		cmd("cp old/candidateInteractions.tsv old/candidateList.top.tsv Results")
	if initialStep > 4:
		cmd("cp -r old/iLoops/Input Results/iLoops/")
		cmd("cp -r old/iLoops/ExpressedTranscripts.fasta Results/iLoops/")
		cmd("cp old/candidates_normal.top.gff old/candidates_tumor.top.gff Results")
	if initialStep > 5:
		cmd("cp -r old/iLoops/Output Results/iLoops")

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

def printParam(initialStep, wd, gaudiWd, minExpression, minCandidateExpression, minPSI, Conditions, Compartments, Replicates, Kmer, top):
	with open("Results/Parameters.cfg", "w") as paramFile:
		paramFile.write("initialStep=" + str(initialStep) + "\n")
		paramFile.write("wd=" + wd + "\n")
		paramFile.write("gaudiWd=" + gaudiWd + "\n")
		paramFile.write("minExpression=" + str(minExpression) + "\n")
		paramFile.write("minCandidateExpression=" + str(minCandidateExpression) + "\n")
		paramFile.write("minPSI=" + str(minPSI) + "\n")
		paramFile.write("Conditions=" + str(Conditions) + "\n")
		paramFile.write("Compartments=" + str(Compartments) + "\n")
		paramFile.write("Replicates=" + str(Replicates) + "\n")
		paramFile.write("Kmer=" + str(Kmer) + "\n")
		paramFile.write("top=" + str(top) + "\n")

def finish():
	
	print("* Moving files to the Results directory and creating a summary tar file.")
	minExpresion = ""
	minPSI = ""

	with open("Results/Parameters.cfg", "r") as paramFile:
		for line in paramFile:
			elements = line.split("=")

			if elements[0] == "minCandidateExpression":
				minExpresion = elements[1].strip()
			elif elements[0] == "minPSI":
				minPSI = elements[1].strip()
	    
	resultsDir="/home/hector/Results/GENECODE19"
	outFolder="e" + minExpresion + "p" + minPSI

	if not cmdOut("mkdir", resultsDir + "/" + outFolder):
		overwrite = raw_input("\tDirectory exists. Do you want to overwrite it? (y/n)")
		if overwrite == "n":
			return

	cmd("cp", "Results/candidates_normal.gff", resultsDir + "/" + outFolder + "/" + "candidates_normal." + outFolder + ".gff")
	cmd("cp", "Results/candidates_tumor.gff", resultsDir + "/" + outFolder + "/" + "candidates_tumor." + outFolder + ".gff")
	cmd("cp", "Results/candidateList.top.tsv", resultsDir + "/" + outFolder + "/" + "candidateList." + outFolder + ".tsv")

	chdir(resultsDir + "/" + outFolder)
	cmd("cp -r ~/SmartAS/Results/* .")

	cmd("tar -czvf", outFolder + ".tar.gz candidates_normal." + outFolder + ".gff candidates_tumor." + outFolder + ".gff candidateList." + outFolder + ".tsv")
	cmd("rm", "*" + outFolder + "*gff", "*" + outFolder + "*tsv")
	chdir("/home/hector/SmartAS")
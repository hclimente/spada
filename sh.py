#!/soft/devel/python-2.7/bin/python

from subprocess import call,Popen,PIPE
from rpy2.robjects import r
import httplib2
import json
from time import sleep

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

def setEnvironment(wd, initialStep, Conditions, Compartments, Replicates, Kmer):

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
	if initialStep > 1:
		cmd("cp -r old/DataExploration Results")
		cmd("cp -r old/RWorkspaces/1_ExploreData.RData Results/RWorkspaces")
		cmd("cp -r old/RWorkspaces/1_ExploreData.RData SmartAS.RData")
		cmd("cp -r old/10C1_30.tsv old/7C1_30.tsv old/10C2_30.tsv old/7C2_30.tsv old/IntraReplicateC1_30.tsv old/IntraReplicateC2_30.tsv Results")
	if initialStep > 2:
		cmd("cp -r old/RWorkspaces/2_GetCandidates.RData Results/RWorkspaces")
		cmd("cp -r old/RWorkspaces/2_GetCandidates.RData SmartAS.RData")
		cmd("cp -r old/candidateList.lst old/candidateList_withGenenames.tsv old/expressedGenes.lst Results")
	if initialStep > 3:
		cmd("cp old/candidateInteractions.tsv old/candidateList.top.lst Results")
	if initialStep > 4:
		cmd("cp -r old/iLoops/Input Results/iLoops/")
		cmd("cp -r old/iLoops/ExpressedTranscripts.fasta Results/iLoops/")
		cmd("cp -r old/candidates.v3.gff old/candidates_normal.v2.gff old/candidates_tumor.v2.gff Results")
	if initialStep > 5:
		cmd("cp -r old/iLoops/Output Results/iLoops")

def waitPID(pidQueue):
	for job in pidQueue:
		while True:
			if cmdOut("ps --pid", job, " | grep -v", job, "| wc -l") == "0":
				break
			else:
				print("Awaiting for completion of iLoops jobs.")
				sleep(900)

def getFeature(query, ensemblId, question):

	http = httplib2.Http(".cache")
	server = "http://beta.rest.ensembl.org"

	ext = "/" + query + "/id/" + ensemblId + "?" + question
	resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
			
	if not resp.status == 200:
		errorCode = {400: "Bad Request (id not found)", 404: "Not Found (incorrect format)", 429: "Too Many Requests", 503: "Service Unavailable"}
		print("\t\tError with " + ensemblId + " (" + server + ext + "): " + errorCode[resp.status])
		return {}
	
	#Ensembl REST API doesn't accept more than 3 queries/second.
	sleep(0.33)

	return json.loads(content)
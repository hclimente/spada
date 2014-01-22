#!/soft/devel/python-2.7/bin/python

from subprocess import call
from rpy2.robjects import r

def cmd(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)
	call(command, shell=True)

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

def setEnvironment(wd, currentStep, Conditions, Compartments, Replicates, Kmer):

	print "* Preparing the environment"
	cmd("cd " + wd)

	if currentStep <= 5:
		cmd("rm -r old/iLoops/output; mkdir -p old/iLoops/output")
		cmd("mv Results/iLoops/output/* old/iLoops/output")
		cmd("mkdir -p Results/iLoops/output")

	if currentStep <= 3:
		cmd("rm -r old/iLoops/input/ENST*; mkdir -p old/iLoops/input")
		cmd("mv Results/iLoops/input/* old/iLoops/input")
		cmd("mv Results/candidates.gff old")
		cmd("mkdir -p Results/iLoops/input")
		cmd("cp Results/RWorkspaces/2_GetCandidates.RData SmartAS.RData")
	if currentStep <= 2:
		cmd("mv Results/expressedGenes.lst old")
		cmd("mv Results/candidateList.lst old")
		cmd("mv Results/RWorkspaces/2_GetCandidates.RData old/RWorkspaces")
		cmd("cp Results/RWorkspaces/1_ExploreData.RData SmartAS.RData")
	if currentStep <= 1:
		cmd("rm -r old")
		cmd("mkdir -p old/DataExploration")
		cmd("mkdir -p old/RWorkspaces")
		cmd("mv Results/DataExploration/* old/DataExploration")
		cmd("mv Results/* old")
		cmd("mv Results/RWorkspaces/* old/RWorkspaces")
		cmd("mv SmartAS.RData old/SmartAS.old.RData")
		cmd("mkdir -p Results/DataExploration")
		cmd("mkdir -p Results/RWorkspaces")
		setRWorkspace(wd, Conditions, Compartments, Replicates, Kmer)
		#cmd("Pipeline/SetWorkspace.r " + wd)
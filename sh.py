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

def setEnvironment(wd, initialStep, Conditions, Compartments, Replicates, Kmer):

	print "* Preparing the environment"
	cmd("cd " + wd)
	cmd("mv Results old; mv SmartAS.RData old/SmartAS.old.RData")
	cmd("mkdir -p Results/iLoops/output")
	cmd("mkdir Results/iLoops/input")
	cmd("mkdir Results/RWorkspaces")
	cmd("mkdir Results/DataExploration")

	if initialStep <= 1:
		setRWorkspace(wd, Conditions, Compartments, Replicates, Kmer)
	if initialStep > 1:
		cmd("mkdir Results/RWorkspaces")
		cmd("cp old/DataExploration Results")
		cmd("cp old/RWorkspaces/1_ExploreData.RData Results/RWorkspaces")
		cmd("cp old/RWorkspaces/1_ExploreData.RData SmartAS.RData")
		cmd("cp old/10C1_30.tsv old/7C1_30.tsv old/10C2_30.tsv old/7C2_30.tsv old/IntraReplicateC1_30.tsv old/IntraReplicateC2_30.tsv Results")
	if initialStep > 2:
		cmd("cp old/RWorkspaces/2_GetCandidates.RData Results/RWorkspaces")
		cmd("cp old/RWorkspaces/2_GetCandidates.RData SmartAS.RData")
		cmd("cp old/candidateList.lst old/expressedGenes.lst Results")
	if initialStep > 3:
		cmd("cp old/iLoops/input Results/iLoops/")
		cmd("cp old/candidates.gff Results")
	if initialStep > 5:
		cmd("cp old/iLoops/output Results/iLoops")
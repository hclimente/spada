#!/soft/devel/python-2.7/bin/python

from subprocess import call

def cmd(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)
	call(command, shell=True)

def setEnvironment(wd, currentStep):
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
		cmd("Pipeline/SetWorkspace.r " + wd)
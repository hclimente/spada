#!/bin/sh

#Set variables
currentStep=$1
wd="/home/hector/SmartAS/"
gaudiWd="/sbi/users/hectorc/SmartAS/Results/iLoops"

echo "#######################################"
echo "#                                     #"
echo "#               SmartAS               #"
echo "#    Finding significant AS events    #"
echo "#                                     #"
echo "#      Hector Climente-GRIB 2014      #"
echo "#                                     #"
echo "#######################################"

main(){

	setEnvironment

	if [[ currentStep -le 1 ]]; then
		exploreData
	fi
	if [[ currentStep -le 2 ]]; then
		getCandidates
	fi
	if [[ currentStep -le 3 ]]; then
		prepareILoopsInput
	fi
	if [[ currentStep -le 4 ]]; then
		launchILoops
	fi
	if [[ currentStep -le 5 ]]; then
		exloreILoopsResults
	fi

	cp -r Results ~/Dropbox/SmartAS
	cp SmartAS.RData ~/Dropbox/SmartAS

}

function setEnvironment(){
	echo "* Preparing the environment"
	cd $wd
	
	if [[ currentStep -le 5 ]]; then
		rm -r old/iLoops/output; mkdir -p old/iLoops/output
		mv Results/iLoops/output/* old/iLoops/output
		mkdir -p Results/iLoops/output
	fi
	if [[ currentStep -le 3 ]]; then
		rm -r old/iLoops/input/ENST*; mkdir -p old/iLoops/input
		mv Results/iLoops/input/* old/iLoops/input
		mv Results/candidates.gff old
		mkdir -p Results/iLoops/input
		cp Results/RWorkspaces/2_GetCandidates.RData SmartAS.RData
	fi
	if [[ currentStep -le 2 ]]; then
		mv Results/expressedGenes.lst old
		mv Results/candidateList.lst old
		mv Results/RWorkspaces/2_GetCandidates.RData old/RWorkspaces
		cp Results/RWorkspaces/1_ExploreData.RData SmartAS.RData
	fi
	if [[ currentStep -le 1 ]]; then
		rm -r old
		mkdir -p old/DataExploration
		mkdir -p old/RWorkspaces
		mv Results/DataExploration/* old/DataExploration
		mv Results/* old
		mv Results/RWorkspaces/* old/RWorkspaces
		mv SmartAS.RData old/SmartAS.old.RData
		mkdir -p Results/DataExploration
		mkdir -p Results/RWorkspaces
		Pipeline/SetWorkspace.r $wd
	fi
	
}

function exploreData(){
	
	echo "* Reading and summarizing input files: computing PSI values and plotting correlations between replicates."
	Pipeline/PSICalculation.r
	cp SmartAS.RData Results/RWorkspaces/1_ExploreData.RData
	
}

function getCandidates(){

	minExpression="0"
	minCandidateExpression="4"
	minPSI="0.25"

	echo "* Extracting transcripts with high variance and high expression."
	Pipeline/GetCandidates.r $minExpression $minCandidateExpression $minPSI
	cp SmartAS.RData Results/RWorkspaces/2_GetCandidates.RData

	sort Results/expressedGenes.lst >Results/expressedGenes.tmp.lst
	mv Results/expressedGenes.tmp.lst Results/expressedGenes.lst

}

function prepareILoopsInput(){

	echo "* Retrieving protein sequences for transcripts and printing to multiFASTA file."
	
	if [[ `diff old/expressedGenes.lst Results/expressedGenes.lst 2>&1` == "" && -e old/iLoops/input/ExpressedTranscripts.fasta ]]
		then 
		getExpressedGenes="0"
	else
		getExpressedGenes="1"
	fi

	Pipeline/getiLoopsInput.pl Results/expressedGenes.lst Results/candidateList.lst $getExpressedGenes

}

function launchILoops(){

	echo "* Launching iLoops jobs."
	
	scp -r Results/iLoops hectorc@gaudi.imim.es:~/SmartAS/Results
	ssh hectorc@gaudi.imim.es '~/SmartAS/Pipeline/launchILoops.sh /sbi/users/hectorc/SmartAS/Results/iLoops' #$gaudiWd'

	echo -e "\t* Waiting..."

}

function exloreILoopsResults(){
	
	echo "* Examining iLoops results."
	Pipeline/exploreiLoopsOutput.py

}

main
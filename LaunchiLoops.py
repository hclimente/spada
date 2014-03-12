#!/soft/devel/python-2.7/bin/python

import sys, os
from include.libsmartas import *
from fnmatch import filter
import include.custom_iLoops_xml_parser as parser
import pdb

def writeFasta(basename, inputType, expressedTranscripts):
	wannaWrite = False
	fileCounter = 1
	transcriptCounter = 0

	expressedTranscriptsSet = set()
	if isinstance(expressedTranscripts, basestring):
		with open(expressedTranscripts, "r") as EXPRESSED:
			for line in EXPRESSED:
				expressedTranscriptsSet.add(line.strip().split("\t")[0])
	elif isinstance(expressedTranscripts, set):
		expressedTranscriptsSet = expressedTranscripts

	MULTIFASTA = open(basename + "_" + str(fileCounter) + ".fasta", "w")

	with open("Data/" + inputType + "/proteins.fa", "r") as gcMULTIFASTA:
		for line in gcMULTIFASTA:
			if line.find(">") != -1:

				if transcriptCounter >= 9999:
					fileCounter += 1
					MULTIFASTA.close()
					MULTIFASTA = open(basename + "_" + str(fileCounter) + ".fasta", "w")
					transcriptCounter = 0

				identifiers = ((line[1:].split("|"))[0].split("."))[0]
				if identifiers in expressedTranscriptsSet:
					MULTIFASTA.write(">" + identifiers + "\n")
					wannaWrite = True
					transcriptCounter += 1
				else:
					wannaWrite = False
			elif wannaWrite:
				MULTIFASTA.write(line)

def parseMapping(iLoopsFolder, tag):
	loopFamilies = {}
	noLoops = set()
	loopModels = set()

	myParser = parser.iLoopsParser()

	for mappingBatch in filter(os.listdir(iLoopsFolder + "Output/"), "Mapping_" + tag + "_*" ):
		xmlFile = iLoopsFolder + "Output/" + mappingBatch + "/" + mappingBatch + ".xml"
		errFile = iLoopsFolder + "Output/" + mappingBatch + "/" + mappingBatch + ".err.xml"

		newLoops = myParser.parseLoops(xmlOutput					 = xmlFile,
									   output_proteins               = True, 
									   output_alignments             = False,
									   output_domain_mappings        = False,
									   output_protein_features       = True,
									   output_domain_assignations    = False,
									   output_interactions           = False,
									   output_interaction_signatures = False,
									   output_RF_results             = False,
									   output_RF_precisions          = False
									  )
		
		for aTranscript, loops in sorted(newLoops.iteritems()):
			if loops not in loopFamilies.keys():
				loopFamilies[loops] = []

			loopFamilies[loops].append(aTranscript)

		newNoLoops = myParser.parseNoLoops(xmlOutput                     = errFile,
										   iLoopsPath                    = iLoopsFolder,
										   output_proteins               = True, 
										   output_alignments             = False,
										   output_domain_mappings        = False,
										   output_protein_features       = False,
										   output_domain_assignations    = False,
										   output_interactions           = False,
										   output_interaction_signatures = False,
										   output_RF_results             = False,
										   output_RF_precisions          = False
										  )
		for aTranscript in newNoLoops:
			noLoops.add(aTranscript)

	with open(iLoopsFolder + tag + "_loopFamilies.txt", "w") as loopFamiliesList:
		 for loopPattern in loopFamilies.keys():
		 	loopModels.add(loopFamilies[loopPattern][0])

		 	loopFamiliesList.write(">" + loopPattern + "\t" + loopFamilies[loopPattern][0] + "\n")
		 	for transcript in loopFamilies[loopPattern][1:]:
		 		loopFamiliesList.write(transcript + "\n")
	
	return (loopFamilies, loopModels, noLoops)

def getFASTAInput(iLoopsFolder, tag, inputType, transcripts):

	writeFasta(iLoopsFolder + tag, inputType, transcripts)
	for expressedFasta in filter(os.listdir(iLoopsFolder), tag + "_*.fasta"):
		
		assignationBatch = expressedFasta.split(".")[0].split("_")[1]
		cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
			"-f " + iLoopsFolder + expressedFasta,
			"-j " + iLoopsFolder + "Output/Mapping_" + tag + "_" + assignationBatch,
			"-x Mapping_" + tag + "_" + assignationBatch + ".xml",
			"-v",
			"-m",
			"-n 25"
		   )

	loopFamilies, loopModels, noLoops = parseMapping(iLoopsFolder, tag)
	writeFasta(iLoopsFolder + tag + "_uniqLoops", inputType, loopModels)

	with open(iLoopsFolder + tag + "_noLoops.li", "a") as NOLOOPS:
		for aTranscript in noLoops:
			NOLOOPS.write(aTranscript + "\n")

	return (loopFamilies, noLoops)

def getPairsInput(iLoopsFolder, goodCandidates):
	for aPair in goodCandidates:
		for aCandidate in aPair:

			cmd("mkdir " + iLoopsFolder + "Input/" + aCandidate)
			for expressedFasta in filter(os.listdir(iLoopsFolder), "Expressed_uniqLoops_*.fasta"):
				cmd("cp", iLoopsFolder + expressedFasta, iLoopsFolder + "Input/" + aCandidate)
				with open(iLoopsFolder + "Input/" + aCandidate + "/" + expressedFasta, "a") as exprFast, \
					 open("Data/" + inputType + "/proteins.fa", "r") as motherFast:
					fileNumber = (expressedFasta.split(".")[0]).split("_")[2]
					writeSeq = False
					
					for line in motherFast:
						if line.find(">") != -1:
							writeSeq = False
							if line.find(aCandidate) != -1:
								writeSeq = True
								exprFast.write(">" + aCandidate + "\n")
						elif writeSeq:
							exprFast.write(line)
				with open(iLoopsFolder + "Input/" + aCandidate + '/' + aCandidate + '_' + fileNumber + '.net', "w") as PAIRS, \
					 open(iLoopsFolder + "Input/" + aCandidate + "/" + expressedFasta, "r") as exprFast:
					for rawExpressed in exprFast:
						if rawExpressed.find(">") == -1:
							continue
						expressedTranscript = rawExpressed.strip()[1:]
						PAIRS.write(aCandidate + "\t" + expressedTranscript + "\n")
		
def getFASTAandPairs(iLoopsFolder, inputType, transcripts):
	candidates = set()
	candidatePairs = []

	with open(transcripts, "r") as CANDIDATES:

		CANDIDATES.readline()
		top = 0

		for line in CANDIDATES:
			elements = line.split("\t")
			candidates.add(elements[2])
			candidates.add(elements[3])
			candidatePairs.append( (elements[2], elements[3]) )
			
			top += 1
			if top >= 30:
				break

	loopFamilies, noLoops = getFASTAInput(iLoopsFolder, "Candidate", inputType, candidates)

	toDelete = set()

	for candidatePair in candidatePairs:
		for candidate in candidatePair:
			if candidate in noLoops:
				toDelete.add(candidatePair)

		for loop in loopFamilies.keys():
			if candidatePair[0] in loopFamilies[loop] and candidatePair[1] in loopFamilies[loop]:
				toDelete.add(candidatePair)

	for failCandidate in toDelete:
		candidatePairs.remove(failCandidate)

	getPairsInput(iLoopsFolder, candidatePairs)

	return candidatePairs

os.chdir(sys.argv[1])
out = "Results/" + sys.argv[2]
inputType = sys.argv[3]
iLoopsFolder = out + "/iLoops/"
expressedTranscripts = out + "/expressedGenes.lst"
candidateTranscripts = out + "/candidateList.top.tsv"

print("\t* Preparing FASTA files for all transcripts.")
getFASTAInput(iLoopsFolder, "Expressed", inputType, expressedTranscripts)

print("\t* Checking loop mapping for top candidates and preparing input.")
goodCandidates = getFASTAandPairs(iLoopsFolder, inputType, candidateTranscripts)

print("\t* Launching iLoops.")

for transcriptPair in goodCandidates:
	for transcript in transcriptPair:
		for configFile in filter(os.listdir(iLoopsFolder + "Input/" + transcript), "*net"):
			batch = (configFile.split(".")[0]).split("_")[1]
			cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
				"-f " + iLoopsFolder + "/Input/" + transcript + "/" + "Expressed_uniqLoops_" + batch + ".fasta",
				"-q " + iLoopsFolder + "/Input/" + transcript + "/" + configFile,
				"-j " + iLoopsFolder + "/Output/" + configFile,
				"-x " + configFile + ".xml",
				"-g all",
				"-n 25",
				"-Q sbi",
				"-c 1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50",
				"-v"
			   )

cmd("scp -r " + iLoopsFolder + " hector@feynman.imim.es:~/SmartAS/" + out)
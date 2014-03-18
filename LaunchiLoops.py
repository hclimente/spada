#!/soft/devel/python-2.7/bin/python

import sys
from os import listdir,chdir
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
				if inputType == "TCGA":
					identifiers = line[1:].strip()
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

	for mappingBatch in filter(listdir(iLoopsFolder + "Output/"), "Mapping_" + tag + "_*" ):
		for mappingBatch_nodeResult in filter(listdir(iLoopsFolder + "Output/" + mappingBatch + "/sge_output"), "*.assignation.[012][0-9].xml" ):
			xmlFile = iLoopsFolder + "Output/" + mappingBatch + "/sge_output/" + mappingBatch_nodeResult

			newLoops = myParser.parseLoops(
											xmlOutput					  = xmlFile,
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
		
		for mappingBatch_nodeErr in filter(listdir(iLoopsFolder + "Output/" + mappingBatch + "/sge_output"), "*.assignation.[012][0-9].err.xml" ):
			errFile = iLoopsFolder + "Output/" + mappingBatch + "/sge_output/" + mappingBatch_nodeErr

			newNoLoops = myParser.parseNoLoops(
												xmlOutput                     = errFile,
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
	for expressedFasta in filter(listdir(iLoopsFolder), tag + "_*.fasta"):
		
		assignationBatch = expressedFasta.split(".")[0].split("_")[1]
		cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
			"-f " + iLoopsFolder + expressedFasta,
			"-j " + iLoopsFolder + "Output/Mapping_" + tag + "_" + assignationBatch,
			"-x Mapping_" + tag + "_" + assignationBatch + ".xml",
			"-g all",
			"-v",
			"-m",
			"-n 25",
			"-Q bigmem",
			"2>&1 >" + iLoopsFolder + "logs/Mapping_" + tag + "_" + assignationBatch + ".log"
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
			for expressedFasta in filter(listdir(iLoopsFolder), "Expressed_uniqLoops_*.fasta"):
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
			top += 1

			elements = line.split("\t")
			candidates.add(elements[2])
			candidates.add(elements[3])
			candidatePairs.append( (elements[2], elements[3]) )
			
			if top >= 30:
				break

	loopFamilies, noLoops = getFASTAInput(iLoopsFolder, "Candidate", inputType, candidates)

	toDelete = set()

	for candidatePair in candidatePairs:
		#delete if none has a loop
		if candidateiPair[0] in noLoops and candidateiPair[1] in noLoops:
			toDelete.add(candidatePair)
		
		#delete if both have the same loops
		for loop in loopFamilies.keys():
			if candidatePair[0] in loopFamilies[loop] and candidatePair[1] in loopFamilies[loop]:
				toDelete.add(candidatePair)

	for failCandidate in toDelete:
		candidatePairs.remove(failCandidate)

	getPairsInput(iLoopsFolder, candidatePairs)

	return candidatePairs

chdir(sys.argv[1])
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
 		for configFile in filter(listdir(iLoopsFolder + "Input/" + transcript), "*net"):
 			batch = (configFile.split(".")[1]).split("_")[1]
 			
			cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
 				"-f " + iLoopsFolder + "Input/" + transcript + "/Expressed_uniqLoops_" + batch + ".fasta",
 				"-q " + iLoopsFolder + "Input/" + transcript + "/" + configFile,
 				"-j " + iLoopsFolder + "Output/" + configFile,
 				"-x " + configFile + ".xml",
 				"-v",
 				"-g all",
 				"-n 25",
 				"-Q bigmem",
 				"-c 1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50",
 				"2>&1 >" + iLoopsFolder + "logs/" + configFile + ".log"
 			   )

		with open(iLoopsFolder + "Output/" + transcript + ".ips", "w") as ISO_OUTPUT:
			ISO_OUTPUT.write("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
			ISO_OUTPUT.write("<xml>\n")
			
			for candidate in filter(listdir(iLoopsFolder + "Output"), transcript + "_*"):
				for tens in range(2):
					for units in range(10):
						number = str(tens) + str(units)
						if number == "00": continue

						
				for nodeMap in filter(listdir(iLoopsFolder + "Output/" + candidate + "/sge_output"), "*.assignation." + number + ".xml" ):
					mapFile = iLoopsFolder + "Output/" + candidate + "/sge_output/" + nodeMap
					with open(mapFile, "r") as MAPPED:
						for line in MAPPED:
							if line.strip() != "<?xml version=\"1.0\" encoding=\"utf-8\"?>" and line.strip() != "<xml>" and line.strip() != "</xml>":
								ISO_OUTPUT.write(line)
			
			for candidate in filter(listdir(iLoopsFolder + "Output"), transcript + "_*"):
				for nodeAss in filter(listdir(iLoopsFolder + "Output/" + candidate + "/sge_output"), "*.scoring.[012][0-9].xml" ):
					assFile = iLoopsFolder + "Output/" + candidate + "/sge_output/" + nodeAss
					with open(assFile, "r") as ASSIGNED:
						for line in ASSIGNED:
							if line.strip() != "<?xml version=\"1.0\" encoding=\"utf-8\"?>" and line.strip() != "<xml>" and line.strip() != "</xml>":
								ISO_OUTPUT.write(line)

			ISO_OUTPUT.write("</xml>\n")

cmd("scp -r " + iLoopsFolder + "Output/*.ips hector@feynman.imim.es:~/SmartAS/" + iLoopsFolder + "Output")

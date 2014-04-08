#!/soft/devel/python-2.7/bin/python

import sys
from os import listdir,chdir
from include.libsmartas import *
from fnmatch import filter
import include.custom_iLoops_xml_parser as parser

def getFASTAInput(basename, inputType, expressedTranscripts):
	wannaWrite = False
	fileCounter = 1
	transcriptCounter = 0

	noLoops = set()

	expressedTranscriptsSet = set()
	if isinstance(expressedTranscripts, basestring):
		with open(expressedTranscripts, "r") as EXPRESSED:
			for line in EXPRESSED:
				expressedTranscriptsSet.add(line.strip().split("\t")[0])
	elif isinstance(expressedTranscripts, set):
		expressedTranscriptsSet = expressedTranscripts

	MULTIFASTA = open(basename + "_" + str(fileCounter) + ".fasta", "w")

	with open("Data/" + inputType + "/UnifiedFasta.fa", "r") as gcMULTIFASTA:
		familiesAdded = {}
		for line in gcMULTIFASTA:
			if ">" in line:

				loopFamily = line.strip().split("#")[3]

				if transcriptCounter >= 1499:
					fileCounter += 1
					MULTIFASTA.close()
					MULTIFASTA = open(basename + "_" + str(fileCounter) + ".fasta", "w")
					transcriptCounter = 0

				identifier = line[1:].strip().split("#")[0]
				if inputType == "GENCODE":
					identifier = ((line[1:].split("|"))[0].split("."))[0]
				
				if identifier in expressedTranscriptsSet and loopFamily != "" and loopFamily not in familiesAdded.keys():
					MULTIFASTA.write(line)
					familiesAdded[loopFamily] = [identifier]
					wannaWrite = True
					transcriptCounter += 1
				else:
					wannaWrite = False
					if loopFamily == "":
						noLoops.add(identifier)
					else:
						if loopFamily not in familiesAdded.keys():
							familiesAdded[loopFamily] = [identifier]
						else:
							familiesAdded[loopFamily].append(identifier)
					
			elif wannaWrite:
				MULTIFASTA.write(line)

	return (familiesAdded, noLoops)

# def parseMapping(iLoopsFolder, tag):
# 	loopFamilies = {}
# 	noLoops = set()
# 	loopModels = set()

# 	myParser = parser.iLoopsParser()

# 	for mappingBatch in filter(listdir(iLoopsFolder + "Output/"), "Mapping_" + tag + "_*" ):
# 		for tens in ["", range(3) ]:
# 			for units in range(10):
# 				number = str(tens) + str(units)
# 				if number == "00": continue
				
# 				for mappingBatch_nodeResult in filter(listdir(iLoopsFolder + "Output/" + mappingBatch + "/sge_output"), "*.assignation." + number + ".xml" ):
# 					xmlFile = iLoopsFolder + "Output/" + mappingBatch + "/sge_output/" + mappingBatch_nodeResult

# 					newLoops = myParser.parseLoops(
# 													xmlOutput					  = xmlFile,
# 													output_proteins               = True, 
# 													output_alignments             = False,
# 													output_domain_mappings        = False,
# 													output_protein_features       = True,
# 													output_domain_assignations    = False,
# 													output_interactions           = False,
# 													output_interaction_signatures = False,
# 													output_RF_results             = False,
# 													output_RF_precisions          = False
# 												  )
					
# 					for aTranscript, loops in sorted(newLoops.iteritems()):
# 						if loops not in loopFamilies.keys():
# 							loopFamilies[loops] = []

# 						loopFamilies[loops].append(aTranscript)
		
# 		for tens in ["", range(3) ]:
# 			for units in range(10):
# 				number = str(tens) + str(units)
# 				if number == "00": continue

# 				for mappingBatch_nodeErr in filter(listdir(iLoopsFolder + "Output/" + mappingBatch + "/sge_output"), "*.assignation." + number + ".err.xml" ):
# 					errFile = iLoopsFolder + "Output/" + mappingBatch + "/sge_output/" + mappingBatch_nodeErr

# 					newNoLoops = myParser.parseNoLoops(
# 														xmlOutput                     = errFile,
# 														iLoopsPath                    = iLoopsFolder,
# 														output_proteins               = True, 
# 														output_alignments             = False,
# 														output_domain_mappings        = False,
# 														output_protein_features       = False,
# 														output_domain_assignations    = False,
# 														output_interactions           = False,
# 														output_interaction_signatures = False,
# 														output_RF_results             = False,
# 														output_RF_precisions          = False
# 													  )
					
# 					for aTranscript in newNoLoops:
# 						noLoops.add(aTranscript)

# 	with open(iLoopsFolder + tag + "_loopFamilies.txt", "w") as loopFamiliesList:
# 		 for loopPattern in loopFamilies.keys():
# 		 	loopModels.add(loopFamilies[loopPattern][0])

# 		 	loopFamiliesList.write(">" + loopPattern + "\t" + loopFamilies[loopPattern][0] + "\n")
# 		 	for transcript in loopFamilies[loopPattern][1:]:
# 		 		loopFamiliesList.write(transcript + "\n")
	
# 	return (loopFamilies, loopModels, noLoops)

# def getFASTAInput(iLoopsFolder, tag, inputType, transcripts):

# 	writeFasta(iLoopsFolder + tag, inputType, transcripts)
# 	for expressedFasta in filter(listdir(iLoopsFolder), tag + "_*.fasta"):
		
# 		assignationBatch = expressedFasta.split(".")[0].split("_")[1]
# 		cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
# 			"-f " + iLoopsFolder + expressedFasta,
# 			"-j " + iLoopsFolder + "Output/Mapping_" + tag + "_" + assignationBatch,
# 			"-x Mapping_" + tag + "_" + assignationBatch + ".xml",
# 			"-g all",
# 			"-v",
# 			"-m",
# 			"-n 25",
# 			"-Q sbi",
# 			"2>&1 >" + iLoopsFolder + "logs/Mapping_" + tag + "_" + assignationBatch + ".log"
# 		   )

# 	loopFamilies, loopModels, noLoops = parseMapping(iLoopsFolder, tag)
# 	writeFasta(iLoopsFolder + tag + "_uniqLoops", inputType, loopModels)

# 	with open(iLoopsFolder + tag + "_noLoops.li", "a") as NOLOOPS:
# 		for aTranscript in noLoops:
# 			NOLOOPS.write(aTranscript + "\n")

# 	return (loopFamilies, noLoops)

def getPairsInput(iLoopsFolder, goodCandidates):
	for aPair in goodCandidates:
		for aCandidate in aPair:

			cmd("mkdir " + iLoopsFolder + "Input/" + aCandidate)
			for expressedFasta in filter(listdir(iLoopsFolder), "Expressed_*.fasta"):
				cmd("cp", iLoopsFolder + expressedFasta, iLoopsFolder + "Input/" + aCandidate)
				with open(iLoopsFolder + "Input/" + aCandidate + "/" + expressedFasta, "a") as exprFast, \
					 open("Data/" + inputType + "/proteins.fa", "r") as motherFast:
					fileNumber = expressedFasta.split(".")[0].split("_")[1]
					writeSeq = False
					
					for line in motherFast:
						if ">" in line:
							writeSeq = False
							if aCandidate in line:
								writeSeq = True
								exprFast.write(">" + aCandidate + "\n")
						elif writeSeq:
							exprFast.write(line)
				with open(iLoopsFolder + "Input/" + aCandidate + '/' + aCandidate + '_' + fileNumber + '.net', "w") as PAIRS, \
					 open(iLoopsFolder + "Input/" + aCandidate + "/" + expressedFasta, "r") as exprFast:
					for rawExpressed in exprFast:
						if not ">" in rawExpressed:
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
			if elements[9] == "No":
				continue
			candidates.add(elements[2])
			candidates.add(elements[3])
			candidatePairs.append( (elements[2], elements[3]) )
			
			if top >= 30:
				break

	loopFamilies, noLoops = getFASTAInput(iLoopsFolder + "Candidate", inputType, candidates)

	toDelete = set()

	for candidatePair in candidatePairs:
		#delete if none has a loop
		if candidatePair[0] in noLoops and candidatePair[1] in noLoops:
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
out = "Results/" + sys.argv[2] + "/"
inputType = sys.argv[3]
iLoopsFolder = out + "iLoops/"
expressedTranscripts = out + "expressedGenes.lst"
candidateTranscripts = out + "candidateList.top.tsv"

print("\t* Preparing FASTA files for all transcripts.")
getFASTAInput(iLoopsFolder + "Expressed", inputType, expressedTranscripts)

print("\t* Checking loop mapping for top candidates and preparing input.")
goodCandidates = getFASTAandPairs(iLoopsFolder, inputType, candidateTranscripts)

print("\t* Launching iLoops.")

for transcriptPair in goodCandidates:
	for transcript in transcriptPair:
 		for configFile in filter(listdir(iLoopsFolder + "Input/" + transcript), "*net"):
 			batch = (configFile.split(".")[1]).split("_")[1]
 			
			cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
 				"-f " + iLoopsFolder + "Input/" + transcript + "/Expressed_" + batch + ".fasta",
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

			maxBatch = len( filter( listdir(iLoopsFolder + "Output"), transcript + "_*") )
			
			for batch in range(1, maxBatch + 1):
				candidateOut = iLoopsFolder + "Output/" + transcript + "_" + str(batch) + "/sge_output"
				for tens in ["", range(3) ]:
					for units in range(10):
						number = str(tens) + str(units)
						if number == "00": continue

						for nodeMap in filter(listdir( candidateOut ), "*.assignation." + number + ".xml" ):
							mapFile =  candidateOut + "/" + nodeMap
							with open(mapFile, "r") as MAPPED:
								for line in MAPPED:
									if line.strip() != "<?xml version=\"1.0\" encoding=\"utf-8\"?>" and line.strip() != "<xml>" and line.strip() != "</xml>":
										ISO_OUTPUT.write(line)

			for batch in range(1, maxBatch + 1):
				candidateOut = iLoopsFolder + "Output/" + transcript + "_" + str(batch) + "/sge_output"
				for tens in ["", range(3) ]:
					for units in range(10):
						number = str(tens) + str(units)
						if number == "00": continue

						for nodeAss in filter(listdir( candidateOut ), "*.scoring." + number + ".xml" ):
							assFile = candidateOut + "/" + nodeAss
							with open(assFile, "r") as ASSIGNED:
								for line in ASSIGNED:
									if line.strip() != "<?xml version=\"1.0\" encoding=\"utf-8\"?>" and line.strip() != "<xml>" and line.strip() != "</xml>":
										ISO_OUTPUT.write(line)

			ISO_OUTPUT.write("</xml>\n")

cmd("scp -r " + iLoopsFolder + "Output/*.ips hector@feynman.imim.es:~/SmartAS/" + iLoopsFolder + "Output")

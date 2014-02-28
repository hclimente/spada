#!/soft/devel/python-2.7/bin/python

import sys, os
from shutil import copy
from libsmartas import *
from fnmatch import filter
import iLoops_xml_parser as parser

class iLoopsParser(parser.ILXMLParser):
	def custom_protein_output(self, protein_object, **kwds): 
		return protein_object
	
	def parseLoops(self, xmlOutput, **kwds):
		parsedLoops = {}

		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds): 
			if isinstance(resultItem, parser.ILXMLProtein):
				if resultItem.get_loops():
					loopList = []
					for aLoop in resultItem.get_loops():
						loopList.append(aLoop.get_code())
					loopList.sort()
					parsedLoops[resultItem.get_name()] = ";".join(loopList)

		return parsedLoops

	def parseNoLoops(self, iLoopsPath, xmlOutput, **kwds):
		with open(iLoopsPath + "noLoops.li", "a") as noLoops:
			for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds): 
				if isinstance(resultItem, parser.ILXMLProtein):
					noLoops.write(resultItem.get_name() + "\n")


def writeFasta(basename, inputType, expressedTranscripts):
	wannaWrite = False
	fileCounter = 1
	transcriptCounter = 0

	expressedTranscriptsSet = set()
	if isinstance(expressedTranscripts, basestring):
		with open(expressedTranscripts, "r") as EXPRESSED:
			for line in EXPRESSED:
				expressedTranscriptsSet.add(line.strip())
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

def parseMapping(iLoopsFolder):
	loopFamilies = {}

	for mappingBatch in filter(os.listdir(iLoopsFolder + "Output/"), "Mapping_*"):
		xmlFile = iLoopsFolder + "Output/" + mappingBatch + "/" + mappingBatch + ".xml"
		errFile = iLoopsFolder + "Output/" + mappingBatch + "/" + mappingBatch + ".err.xml"

		newLoops = myParser.parseLoops(xmlOutput					   = xmlFile,
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

		myParser.parseNoLoops(xmlOutput                     = errFile,
							  iLoopsPath 					= iLoopsFolder,
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

	loopPatternModel = set()
	with open(iLoopsFolder + "loopFamilies.txt", "w") as loopFamiliesList:
		 for loopPattern in loopFamilies.keys():
		 	loopPatternModel.add(loopFamilies[loopPattern][0])

		 	loopFamiliesList.write(">" + loopPattern + "\t" + loopFamilies[loopPattern][0] + "\n")
		 	for transcript in loopFamilies[loopPattern][1:]:
		 		loopFamiliesList.write(transcript + "\n")

	return loopPatternModel

os.chdir(sys.argv[1])
out = "Results/" + sys.argv[2]
expressedTranscripts = out + sys.argv[3]
candidateTranscripts = out + sys.argv[4]
inputType = sys.argv[5]
iLoopsFolder = out + "/iLoops/"
myParser = iLoopsParser()

print("\t* Preparing FASTA files.")

writeFasta(iLoopsFolder + "ExpressedTranscripts", inputType, expressedTranscripts)

for expressedFasta in filter(os.listdir(iLoopsFolder), "ExpressedTranscripts_*.fasta"):
	
	assignationBatch = expressedFasta.split("_")[1].split(".")[0]
	cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
		"-f " + iLoopsFolder + expressedFasta,
		"-j " + iLoopsFolder + "Output/Mapping_" + assignationBatch,
		"-x Mapping_" + assignationBatch + ".xml",
		"-v",
		"-m",
		"-n 25"
	   )

loopPatternReps = parseMapping(iLoopsFolder)
writeFasta(iLoopsFolder + "ExpressedTranscripts.uniqLoops", inputType, loopPatternReps)

print("\t* Writing the pairs files.")
top = 0

with open(candidateTranscripts, "r") as CANDIDATES:
	CANDIDATES.readline()
	for line in CANDIDATES:
		elements = line.split("\t")
		
		for aCandidate in [elements[2], elements[3]]:
	
			cmd("mkdir " + out + "/iLoops/Input/" + aCandidate)

			for expressedFasta in filter(os.listdir(iLoopsFolder), "ExpressedTranscripts.uniqLoops_*.fasta"):
				cmd("cp", iLoopsFolder + expressedFasta, iLoopsFolder + "Input/" + aCandidate)
				with open(iLoopsFolder + "Input/" + aCandidate + "/" + expressedFasta, "a") as exprFast, \
					 open("Data/" + inputType + "/proteins.fa", "r") as motherFast:

					 fileNumber = (expressedFasta.split(".")[1]).split("_")[1]

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
						
		top +=1
		if top >=30:
			break

print("\t* Launching iLoops.")
for transcript in filter(os.listdir(iLoopsFolder + "/Input"), "ENST*"):
	for configFile in filter(os.listdir(iLoopsFolder + "/Input/" + transcript), "*net"):
		batch = (configFile.split(".")[0]).split("_")[1]
		cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
			"-f " + iLoopsFolder + "/Input/" + transcript + "/" + "ExpressedTranscripts.uniqLoops_" + batch + ".fasta",
			"-q " + iLoopsFolder + "/Input/" + transcript + "/" + configFile,
			"-j " + iLoopsFolder + "/Output/" + configFile,
			"-x " + configFile + ".xml",
			"-g all",
			"-n 25",
			"-Q sbi",
			"-c 1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50",
			"-v"
			)

cmd("scp -r " + iLoopsFolder + " hector@feynman.imim.es:~/SmartAS/Results/")

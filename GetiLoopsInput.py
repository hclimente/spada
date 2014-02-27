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
	with open(expressedTranscripts, "r") as EXPRESSED:
		for line in EXPRESSED:
			expressedTranscriptsSet.add(line.strip())

	MULTIFASTA = open(basename + "_" + str(fileCounter) + ".fasta", "w")

	with open("Data/" + inputType + "/proteins.fa", "r") as gcMULTIFASTA:
		for line in gcMULTIFASTA:
			if line.find(">") != -1:

				if transcriptCounter >= 100:
					break
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

out = sys.argv[1]
expressedTranscripts = out + sys.argv[2]
candidateTranscripts = out + sys.argv[3]
inputType = sys.argv[4]
iLoopsFolder = out + "/iLoops/"
myParser = iLoopsParser()
os.chdir("/sbi/users/hectorc/SmartAS")

print("\t* Preparing original FASTA file.")

writeFasta(iLoopsFolder + "ExpressedTranscripts", inputType, expressedTranscripts)

isoformSeq = {}
for expressedFasta in filter(os.listdir(iLoopsFolder), "ExpressedTranscripts_*.fasta"):
	with open(iLoopsFolder + "/" + expressedFasta, "r") as expFasta:
		currentTranscript = ""
		for line in expFasta:
			trueLine = line.strip()
			if trueLine.find(">") != -1:
				currentTranscript = trueLine[1:]
				isoformSeq[currentTranscript] = ""
			else:
				isoformSeq[currentTranscript] += trueLine
	
	assignationBatch = expressedFasta.split("_")[1].split(".")[0]
	cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
		"-f " + iLoopsFolder + expressedFasta,
		"-j " + iLoopsFolder + "Output/Mapping_" + assignationBatch,
		"-x Mapping_" + assignationBatch + ".xml",
		"-v",
		"-m",
		"-n 25"
	   )

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

with open(iLoopsFolder + "expressedTranscripts.uniqLoops.li", "w") as uniqLoops, \
	 open(iLoopsFolder + "loopFamilies.txt", "w") as loopFamiliesList:
	 for loopPattern in loopFamilies.keys():
	 	uniqLoops.write(loopFamilies[loopPattern][0] + "\n")

	 	loopFamiliesList.write(">" + loopPattern + "\n")
	 	for transcript in loopFamilies[loopPattern]:
	 		loopFamiliesList.write(transcript + "\n")

writeFasta(iLoopsFolder + "ExpressedTranscripts.uniqLoops", inputType, iLoopsFolder + "expressedTranscripts.uniqLoops.li")

print("\t* Writing the pairs files.")

with open(candidateTranscripts, "r") as CANDIDATES:
	CANDIDATES.readline()
	for line in CANDIDATES:
		elements = line.split("\t")
		
		for aCandidate in [elements[2], elements[3]]:
	
			fileNumber = 1
			numberOfCandidates = 0
			cmd("mkdir " + out + "/iLoops/Input/" + aCandidate)

			print aCandidate

			for expressedFasta in filter(os.listdir(iLoopsFolder), "ExpressedTranscripts.uniqLoops_*.fasta"):
				cmd("cp", iLoopsFolder + expressedFasta, iLoopsFolder + "/Input/" + aCandidate)
				with open(iLoopsFolder + "/Input/" + aCandidate + "/" + expressedFasta, "a") as exprFast, \
					 open("Data/" + inputType + "/proteins.fa", "r") as motherFast:
					 writeSeq = False
					 for line in motherFast:
					 	if line.find(">") != -1:
					 		writeSeq = False
					 		if line.find(aCandidate) != -1:
					 			writeSeq = True
						 		exprFast.write(">" + aCandidate + "\n")
						elif writeSeq:
							exprFast.write(line)

				PAIRS = open(out + "/iLoops/Input/" + aCandidate + '/' + aCandidate + '_' + str(fileNumber) + '.net', "w")
				with open(iLoopsFolder + "/Input/" + aCandidate + "/" + expressedFasta, "r") as exprFast:
					for rawExpressed in exprFast:
						expressedTranscript = rawExpressed.strip()
						PAIRS.write(aCandidate + "\t" + expressedTranscript + "\n")
						numberOfCandidates += 1
						if(numberOfCandidates >= 10000):
							fileNumber += 1
							numberOfCandidates = 0
							PAIRS.close()
							PAIRS = open(out + "/iLoops/Input/" + aCandidate + '/' + aCandidate + '_' + str(fileNumber) + '.net', "w")
					PAIRS.close()

		break

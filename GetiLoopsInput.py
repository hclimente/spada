#!/soft/devel/python-2.7/bin/python

import sys
from shutil import copy
from libsmartas import *
import os

out = sys.argv[1]
expressedTranscripts = "Results/" + out + sys.argv[2]
candidateTranscripts = "Results/" + out + sys.argv[3]

expressedTranscriptsSet = set()
with open(expressedTranscripts, "r") as EXPRESSED:
	for line in EXPRESSED:
		expressedTranscriptsSet.add(line.strip())

wannaWrite = False
expressedTranscriptsCounter = len(expressedTranscriptsSet)

with open("Data/GENCODE/proteins.fa", "r") as gcMULTIFASTA:
	with open("Results/" + out + "/iLoops/ExpressedTranscripts.fasta", "w") as MULTIFASTA:
		for line in gcMULTIFASTA:
			if line.find(">") != -1:
				identifiers = ((line.split("|"))[0].split("."))[0]
				if identifiers in expressedTranscriptsSet:
					MULTIFASTA.write(identifiers + "\n")
					wannaWrite = True
				else:
					wannaWrite = False
			elif wannaWrite:
				MULTIFASTA.write(line)

print("\t* Writing the pairs files.")

with open(candidateTranscripts, "r") as CANDIDATES:
	CANDIDATES.readline()
	for line in CANDIDATES:
		elements = line.split("\t")
		if len(elements) == 4:
			break
		
		for aCandidate in [elements[2], elements[3]]:
	
			cmd("mkdir Results/" + out + "/iLoops/Input/" + aCandidate)
			fileNumber = 1
			numberOfCandidates = 0
			with open(expressedTranscripts, "r") as EXPRESSED:
				PAIRS = open("Results/" + out + "/iLoops/Input/" + aCandidate + '/' + aCandidate + '_' + str(fileNumber) + '.net', "w")
				for rawExpressed in EXPRESSED:
					expressedTranscript = rawExpressed.strip()
					PAIRS.write(aCandidate + "\t" + expressedTranscript + "\n")
					numberOfCandidates += 1
					if(numberOfCandidates >= 10000):
						fileNumber += 1
						numberOfCandidates = 0
						PAIRS.close()
						PAIRS = open("Results/" + out + "/iLoops/Input/" + aCandidate + '/' + aCandidate + '_' + str(fileNumber) + '.net', "w")
				PAIRS.close()
#!/soft/devel/python-2.7/bin/python

from interface import iLoops_parser as parser
from interface import iLoops_outputPruner as pruner
from libs import utils

import sys
import os
import fnmatch

def splitFASTA(basename, inputType, iLoopsVersion):
	wannaWrite = False
	fileCounter = 1
	transcriptCounter = 0

	familiesAdded = set()

	MULTIFASTA = open(basename + str(fileCounter) + ".fasta", "w")

	with open("Data/" + inputType + "/UnifiedFasta_" + iLoopsVersion + ".fa", "r") as gcMULTIFASTA:
		
		for line in gcMULTIFASTA:
			if ">" in line:

				wannaWrite = False
				loopFamily = line.strip().split("#")[3]
				identifier = line[1:].strip().split("#")[0]
				if inputType == "GENCODE":
					identifier = ((line[1:].split("|"))[0].split("."))[0]

				if transcriptCounter >= 1499:
					fileCounter += 1
					MULTIFASTA.close()
					MULTIFASTA = open(basename + str(fileCounter) + ".fasta", "w")
					transcriptCounter = 0
				
				if loopFamily and loopFamily not in familiesAdded:
					familiesAdded.add(loopFamily)
					wannaWrite = True
					transcriptCounter += 1
					MULTIFASTA.write(">" + identifier + "\n")
					
			elif wannaWrite:
				MULTIFASTA.write(line)

def getFASTAandPairs(transcript, batch, inputType, out, iLoopsVersion):
	tag = transcript + "_" + batch
	utils.cmd("cp", out + "allProteome_" + batch + ".fasta", out + "Input/" + tag + ".fasta")

 	with open(out + "Input/" + tag + ".fasta", "r") as FASTUKI, open(out + "Input/" + tag + ".net", "w") as PAIRSUKI:
 		for line in FASTUKI:
 			if ">" in line:
 				pair = line[1:].strip()

 				PAIRSUKI.write(transcript + "\t" + pair + "\n")

 	with open("Data/" + inputType + "/UnifiedFasta_" + iLoopsVersion + ".fa", "r") as gcMULTIFASTA, \
 		 open(out + "Input/" + tag + ".fasta", "a") as THIS_FASTA:
 		wannaWrite = False
 		for line in gcMULTIFASTA:
 			if ">" in line:
 				wannaWrite = False
 				if transcript in line:
 					THIS_FASTA.write(">" + transcript + "\n")
 					wannaWrite = True

 			elif wannaWrite:
 				THIS_FASTA.write(line)

os.chdir(sys.argv[1])
out = "Results/" + sys.argv[2] + "/"
inputType = sys.argv[3]
iLoopsVersion = sys.argv[4]

print("\t* Preparing FASTA files for all transcripts.")
splitFASTA(out + "allProteome_", inputType, iLoopsVersion)

print("\t* Launching iLoops.")

for line in utils.readTable(out + "candidatesGaudi.lst", header=False):
	transcript = line[0]
	analysisCode = line[1]
	
	if analysisCode != "0":
		continue
	
 	for fastaFile in fnmatch.filter(os.listdir(out), "allProteome_*.fasta"):
 		batch = (fastaFile.split(".")[0]).split("_")[1]
 		tag = transcript + "_" + batch
 		getFASTAandPairs(transcript, batch, inputType, out, iLoopsVersion)
 			
		utils.cmd("/soft/devel/python-2.7/bin/python /sbi/programs/" + iLoopsVersion + "/iLoops.py",
 			"-f " + out + "Input/" + tag + ".fasta",
 			"-q " + out + "Input/" + tag + ".net",
 			"-j " + out + "Output/" + tag,
 			"-x " + tag + ".xml",
 			"-v",
 			"-g all",
 			"-n 25",
 			"-Q sbi",
 			"-c 1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50",
 			"2>&1 >" + out + "logs/" + tag + ".log"
 		   )

	p = pruner.iLoopsOutput_pruner(transcript, out + "Output/")
	p.joinFiles()
	if p.makeLiteVersion():
		utils.cmd("scp","-r", out + "Output/" + transcript + ".tar.gz", "hector@einstein.imim.es:~/SmartAS/iLoops/" + inputType + "/" + iLoopsVersion)
	else:
		print("Error in generation of file.")
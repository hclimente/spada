#!/soft/devel/python-2.7/bin/python

import sys
from os import listdir,chdir
from include.libsmartas import *
from fnmatch import filter
import include.custom_iLoops_xml_parser as parser

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
	cmd("cp", out + "allProteome_" + batch + ".fasta", out + "Input/" + tag + ".fasta")

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

chdir(sys.argv[1])
out = "Results/" + sys.argv[2] + "/"
inputType = sys.argv[3]
iLoopsVersion = sys.argv[4]

print("\t* Preparing FASTA files for all transcripts.")
splitFASTA(out + "allProteome_", inputType, iLoopsVersion)

print("\t* Launching iLoops.")

with open(out + "candidatesGaudi.lst", "r") as CANDIDATES:
	for line in CANDIDATES:
		elements = line.strip().split("\t")
		transcript = elements[0]
		analysisCode = elements[1]
		
		if analysisCode != "0":
			continue
		
	 	for fastaFile in filter(listdir(out), "allProteome_*.fasta"):
	 		batch = (fastaFile.split(".")[0]).split("_")[1]
	 		tag = transcript + "_" + batch
	 		getFASTAandPairs(transcript, batch, inputType, out, iLoopsVersion)
	 			
			cmd("/soft/devel/python-2.7/bin/python /sbi/programs/" + iLoopsVersion + "/iLoops.py",
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

		with open(out + "Output/" + transcript + ".ips", "w") as ISO_OUTPUT:
			ISO_OUTPUT.write("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
			ISO_OUTPUT.write("<xml>\n")

			maxBatch = len( filter( listdir(out + "Output"), transcript + "_*") )
				
			for batch in range(1, maxBatch + 1):
				candidateOut = out + "Output/" + transcript + "_" + str(batch) + "/sge_output"
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
				candidateOut = out + "Output/" + transcript + "_" + str(batch) + "/sge_output"
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

		cmd("scp","-r", out + "Output/" + transcript + ".ips", "hector@feynman.imim.es:~/SmartAS/iLoops/" + inputType + "/" + iLoopsVersion)
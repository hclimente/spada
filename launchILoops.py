#!/soft/devel/python-2.7/bin/python

import sys, os
from libsmartas import *
from fnmatch import filter
import iLoopsXMLParser as parser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

if(len(sys.argv) != 2):
	print("No arguments passed.")
	exit()

class iLoopsParser(parser.ILXMLParser):
	def custom_protein_output(self, protein_object, **kwds): 
		return protein_object
	
	def parseResults(self, xmlOutput, **kwds):
		parsedLoops = {}

		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=1, **kwds): 
			if isinstance(resultItem, parser.ILXMLProtein):
				print(resultItem.get_name())
				loopList = []
				for aLoop in resultItem.get_loops():
					loopList.append(aLoop.get_code())

				if loopList:
					parsedLoops[resultItem.get_name()] = ";".join(loopList.sort())
				
		print parsedLoops

		return parsedLoops

os.chdir(sys.argv[1])

pidQueue = []

isoformSeq = {}

for expressedFasta in filter(os.listdir("."), "ExpressedTranscripts_*.fasta"):
	with open(expressedFasta, "r") as expFasta:
		currentTranscript = ""
		for line in expFasta:
			trueLine = line.strip()
			if trueline.find(">"):
				currentTranscript = trueline[1:]
				isoformSeq[currentTranscript] = ""
			else:
				isoformSeq[currentTranscript] += trueLine
	
	assignationBatch = expressedFasta.split("_")[1].split(".")[0]
	# cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
	# 	"-f " + expressedFasta,
	# 	"-j Output/Mapping_" + assignationBatch,
	# 	"-x Mapping_" + assignationBatch + ".xml",
	# 	"-v",
	# 	"-m",
	# 	"-n 25"
	# 	)

myParser = iLoopsParser()
allTranscripts = {}

for mappingBatch in filter(os.listdir("Output/"), "Mapping_*"):
	for xmlFile in filter(os.listdir("Output/" + mappingBatch + "/sge_output"), "*assignation.??.xml"):
		newLoops = myParser.parseResults(xmlOutput=xmlFile, outputInteraction_signatures=True, outputRFPrecisions=True)
		allTranscripts = dict(allTranscripts.items() + newLoops.items())

print(allTranscripts)

loopPatterns = []
loopFamilies = {}

for transcript, loops in sorted(allTranscripts.iteritems()):
	if loops in loopPatterns:
		loopFamilies[aTranscript].append(aTranscript)
	else:
		loopPatterns.append(loops)
		loopFamilies[aTranscript] = []

with open("ExpressedTranscripts.loopFiltered.fasta", "w") as FILTERED:
	for fasta in isoformSeq:
		for representative in loopPatterns.keys():
			if representative != fasta.id:
				break
			FILTERED.write(">" + fasta.id + "\n" + fasta.seq.tostring() + "\n")

# for transcript in filter(os.listdir("Input"), "ENST*"):
# 	for configFile in filter(os.listdir("Input/" + transcript), "*net"):
# 		#Map the loops for the query sequence.
# 		print transcript
# 		print configFile
# 		cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
# 			"-f ExpressedTranscripts.fasta",
# 			"-q Input/" + transcript + "/" + configFile,
# 			"-j Output/" + configFile,
# 			"-x " + configFile + ".xml",
# 			"-g all",
# 			"-n 25",
# 			"-Q sbi",
# 			"-c 1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50",
# 			"-v",
# 			"-m &"
# 			)
# 	 	pidQueue.append(cmdOut("$!"))
# 		exit()

# waitPID(pidQueue)

cmd("scp -r Results/iLoops/Output hectorc@gaudi.imim.es:~/SmartAS/Results/iLoops")

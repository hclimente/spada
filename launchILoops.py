#!/soft/devel/python-2.7/bin/python

import sys, os
from libsmartas import *
from fnmatch import filter
import iLoops_xml_parser as parser

class iLoopsParser(parser.ILXMLParser):
	def custom_protein_output(self, protein_object, **kwds): 
		return protein_object
	
	def parseResults(self, xmlOutput, **kwds):
		parsedLoops = {}

		with open("noLoops.li", "w") as noLoops:
			for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds): 
				if isinstance(resultItem, parser.ILXMLProtein):
					if resultItem.get_loops():
						loopList = []
						for aLoop in resultItem.get_loops():
							loopList.append(aLoop.get_code())
						loopList.sort()
						parsedLoops[resultItem.get_name()] = ";".join(loopList)
					else:
						noLoops.write(resultItem.get_name() + "\n")


		return parsedLoops

if(len(sys.argv) != 2):
	print("No arguments passed.")
	exit()

os.chdir(sys.argv[1])

isoformSeq = {}
for expressedFasta in filter(os.listdir("."), "ExpressedTranscripts_*.fasta"):
	with open(expressedFasta, "r") as expFasta:
		currentTranscript = ""
		for line in expFasta:
			trueLine = line.strip()
			if trueLine.find(">") != -1:
				currentTranscript = trueLine[1:]
				isoformSeq[currentTranscript] = ""
			else:
				isoformSeq[currentTranscript] += trueLine
	
	assignationBatch = expressedFasta.split("_")[1].split(".")[0]
#	cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
#		"-f " + expressedFasta,
#		"-j Output/Mapping_" + assignationBatch,
#		"-x Mapping_" + assignationBatch + ".xml",
#		"-v",
#		"-m",
#		"-n 25"
#	   )

myParser = iLoopsParser()
allTranscripts = {}

for mappingBatch in filter(os.listdir("Output/"), "Mapping_*"):
	for xmlFile in filter(os.listdir("Output/" + mappingBatch + "/sge_output"), "*assignation.??.xml"):
		filePath = "Output/" + mappingBatch + "/sge_output/" + xmlFile
		newLoops = myParser.parseResults(xmlOutput=filePath, 
										 output_proteins               = True, 
                                         output_alignments             = False,
                                         output_domain_mappings        = False,
                                         output_protein_features       = True,
                                         output_domain_assignations    = False,
                                         output_interactions           = False,
                                         output_interaction_signatures = False,
                                         output_RF_results             = False,
                                         output_RF_precisions          = False)
		allTranscripts = dict(allTranscripts.items() + newLoops.items())

loopFamilies = {}

for aTranscript, loops in sorted(allTranscripts.iteritems()):
	if loops not in loopFamilies.keys():
		loopFamilies[loops] = []

	loopFamilies[loops].append(aTranscript)

a = 0
for loop in loopFamilies.keys():
	a += len(loopFamilies[loop])

print a

with open("ExpressedTranscripts.loopFiltered.fasta", "w") as FILTERED:
	for representative in loopFamilies.keys():
		FILTERED.write(">" + loopFamilies[representative][0] + "\n" + isoformSeq[loopFamilies[representative][0]] + "\n")

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

#!/soft/devel/python-2.7/bin/python

import sys, os
from sh import *
from fnmatch import filter
import iLoopsXMLParser as parser

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

os.chdir(sys.argv[1])

pidQueue = []

cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
	"-f ExpressedTranscripts.fasta",
	"-j Output/" + configFile,
	"-x " + configFile + ".xml",
	"-g all",
	"-n 25",
	"-Q sbi",
	"-c 1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50",
	"-v",
	"-m"
	)

myParser = iLoopsParser()
xmlFile = "output/ENST00000243253_3.net/sge_output/22939.assignation.01.xml"

myParser.parseResults(xmlOutput=xmlFile, outputInteraction_signatures=True, outputRFPrecisions=True)

# for transcript in filter(os.listdir("Input"), "ENST*"):
# 	for configFile in filter(os.listdir("Input/" + transcript), "*net"):
# 		#Map the loops for the query sequence.
# 		print transcript
# 		print configFile
# 		cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
# 			"-f ExpressedTranscripts.fasta",
# 			"-q Input/" + transcript + "/" + configFile,
# 			"-j output/" + configFile,
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
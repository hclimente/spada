#!/soft/devel/python-2.7/bin/python

import sys, os
from sh import *
from fnmatch import filter
import iLoops_xml_parser as parser

if(len(sys.argv) != 2):
	print("No arguments passed.")
	exit()

class iLoopsParser(parser.ILXMLParser):
	def parseResults(self, xmlOutput, **kwds):
		parsedProteins = [] # to store parsed proteins
		reportLevel = 1

		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=reportLevel, **kwds): 
			# process the customized result item
			# if it is a protein, process the object as desired (print to STDOUT the protein and store it in parsedProteins list)
			if isinstance(resultItem, parser.ILXMLLoop): 
				parsedProteins.append(resultItem)
				print resultItem
				
			# else, just print it!
			else: 
				print resultItem

		# demonstrate that the full protein object has been stored into parsedProteins list
		if len(parsedProteins) > 0: 
			print "PARSED PROTEINS:", len(parsedProteins), "\t", ", ".join([ x.get_name() for x in parsedProteins ])

os.chdir(sys.argv[1])

pidQueue = []

# for transcript in filter(os.listdir("input"), "ENST*"):
# 	for configFile in filter(os.listdir("input/" + transcript), "*net"):
# 		#Map the loops for the query sequence.
# 		print transcript
# 		print configFile
# 		cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
# 			"-f ExpressedTranscripts.fasta",
# 			"-q input/" + transcript + "/" + configFile,
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

myParser = iLoopsParser()
xmlFile = "output/ENST00000243253_3.net/sge_output/22939.assignation.01.xml"

myParser.parseResults(xml_file=file, output_interaction_signatures = False, output_RF_precisions= False)



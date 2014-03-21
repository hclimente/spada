#!/soft/devel/python-2.7/bin/python

from include.libsmartas import *
import include.custom_iLoops_xml_parser as parser
from os import listdir
from os.path import isfile
from fnmatch import filter
import sys

import pdb

out = "Results/" + sys.argv[1] + "/"
myParser = parser.iLoopsParser()

InteraX = {}
expressedTranscripts = {}

with open(out + "expressedGenes.lst", "r") as EXPRESSED:
	for line in EXPRESSED:
		elements = line.strip().split("\t")
		expressedTranscripts[elements[0]] = elements[1]

with open(out + "candidateList.top.tsv", "r") as CANDIDATES, open(out + "InteraxChanges.tsv", "w") as INTERAX:

	INTERAX.write("Gene" + "\t" + "Origin" + "\t" + "Transcript" + "\t" + "Interactions" + "\n")
	CANDIDATES.readline()
	top = 0
	
	for line in CANDIDATES:
		elements = line.split("\t")
		interactions = {}
		InteraX = {}
		top += 1
		
		if not isfile(out + "iLoops/Output/" + elements[2] + ".ips") and not isfile(out + "iLoops/Output/" + elements[3] + ".ips"):
			print("\t * Gene " + elements[1] + ": iLoops couldn't process this candidate.")
			continue

		print("\t * Gene " + elements[1])
				
		for iso, ori in zip([elements[2], elements[3]], ["Normal","Tumor"]):

			INTERAX.write(elements[0] + "\t" + ori + "\t" + iso + "\t")

			print("\t\t" + ori + ": " + iso)

			xmlFile = out + "iLoops/Output/" + iso + ".ips"

			if not isfile(xmlFile):
				interactions[ori] = {}
			else:
				interactions[ori] = myParser.parseInteractions(
															 	thisCandidate				  = iso,
															 	xmlOutput					  = xmlFile,
															 	output_proteins               = False, 
															 	output_alignments             = False,
															 	output_domain_mappings        = False,
															 	output_protein_features       = False,
															 	output_domain_assignations    = False,
															 	output_interactions           = True,
															 	output_interaction_signatures = True,
															 	output_RF_results             = True,
															 	output_RF_precisions          = False
															  )

			for partner in interactions[ori].keys():
				if ori == "Normal":
					InteraX[partner] = InteraX.get(partner, 0) + interactions[ori][partner]
				elif ori == "Tumor":
					InteraX[partner] = InteraX.get(partner, 0) - interactions[ori][partner]
				
				INTERAX.write(partner + "(" + str(interactions[ori][partner]) + ")" + ";")

			INTERAX.write("\n")

		with open(out + elements[0] + ".dot", "w") as DOTFile:
			for iso, ori in zip([elements[2], elements[3]], ["Normal","Tumor"]):
		 		DOTFile.write("graph " + iso.split(".")[0] + " {\n")
		 		DOTFile.write("\t" + iso.split(".")[0] + " [label=" + elements[0] + ", shape=polygon,sides=5];\n")
		 		for partner in interactions[ori].keys():
		 			DOTFile.write("\t" + partner.split(".")[0] + " [shape=record,label=\"<f0> "+ partner.split(".")[0] +"|<f1> " + expressedTranscripts[partner] + "\"];\n") 

		 		for partner in interactions[ori].keys():
		 			DOTFile.write("\t" + iso.split(".")[0] + " -- " + partner.split(".")[0] + " [label=\"" + str(interactions[ori][partner]) + "\"") 
		 			if InteraX[partner] >= 6:
		 				DOTFile.write(", style=bold, color=blue, weight=" + str(InteraX[partner]) )
		 			elif InteraX[partner] <= -6:
		 				DOTFile.write(", style=bold, color=blue, weight=" + str(InteraX[partner]) )
		 			DOTFile.write("];\n")
		 		DOTFile.write("}\n")

		cmd("circo", "-Tps", out + elements[0] + ".dot", "-o", out + elements[0] + ".ps")

		if top >= 30:
			break
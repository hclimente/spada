#!/soft/devel/python-2.7/bin/python

from include.libsmartas import *
import include.custom_iLoops_xml_parser as parser
from os import listdir
from fnmatch import filter
import sys

out = "Results/" + sys.argv[1] + "/"
myParser = parser.iLoopsParser()

InteraX = {}
expressedTranscripts = {}

with open(out + "expressedGenes.lst", "r") as EXPRESSED:
	for line in EXPRESSED:
		elements = line.strip().split("\t")
		expressedTranscripts[elements[0]] = elements[1]

with open(out + "candidateList.top.tsv", "r") as CANDIDATES:

	CANDIDATES.readline()
	top = 0
	
	for line in CANDIDATES:
		elements = line.split("\t")
		interactions = {}
		InteraX = {}
		top += 1
		
		if not filter(listdir(out + "/iLoops/Output"), elements[2] + "_*") and not filter(listdir(out + "/iLoops/Output"), elements[3] + "_*"):
			print("\t * Gene " + elements[1] + ": iLoops couldn't process this candidate.")
			continue

		print("\t * Gene " + elements[1])
		
		for iso, ori in zip([elements[2], elements[3]], ["Normal","Tumor"]):

			interactions[ori] = {}
			print("\t\t" + ori + ": " + iso)
			for candidate in filter(listdir(out + "/iLoops/Output"), iso + "_*"):
				#print out + "/iLoops/Output/" + candidate + "/" + candidate + ".xml"

				for nodeResult in filter(os.listdir(iLoopsFolder + "Output/" + candidate + "/sge_output"), "*.assignation.[012][0-9].xml" ):
					xmlFile = iLoopsFolder + "Output/" + candidate + "/sge_output/" + nodeResult

					candidateTag = candidate.split("_")[0]
					newInteractions = myParser.parseInteractions(
																 	thisCandidate				   = candidateTag,
																 	xmlOutput					   = xmlFile,
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

					interactions[ori].update(newInteractions)
					for partner in newInteractions.keys():
						if ori == "Normal":
							InteraX[partner] = InteraX.get(partner, 0) + newInteractions[partner]
						elif ori == "Tumor":
							InteraX[partner] = InteraX.get(partner, 0) - newInteractions[partner]
		
		for iso, ori in zip([elements[2], elements[3]], ["Normal","Tumor"]):
		 	with open(out + elements[0] + ".dot", "w") as DOTFile:
		 		DOTFile.write("graph " + iso + " {\n")
		 		DOTFile.write("\t" + iso + " [label=" + elements[0] + ", shape=polygon,sides=5];\n")
		 		for partner in interactions[ori].keys():
		 			DOTFile.write("\t" + partner + " [shape=record,label=\"<f0> "+ partner +"|<f1> " + expressedTranscripts[partner] + "\"];\n") 

		 		for partner in interactions[ori].keys():
		 			DOTFile.write("\t" + iso + " -- " + partner + " [label=\"" + str(interactions[ori][partner]) + "\"") 
		 			if InteraX[partner] >= 6:
		 				DOTFile.write(", style=bold, color=blue, weight=" + str(InteraX[partner]) )
		 			elif InteraX[partner] <= -6:
		 				DOTFile.write(", style=bold, color=blue, weight=" + str(InteraX[partner]) )
		 			DOTFile.write("];\n")
		 		DOTFile.write("}\n")

		if top >= 3:
			break

with open(out + "InteraxChanges.tsv", "w") as INTERAX:
	INTERAX.write("Gene" + "\t" + "Normal" + "\t" + "Tumor" + "\t" + "Loss in tumor" + "\t" + "Gain in tumor" + "\n")
	for candidate in InteraX.keys():
		elements = candidate.split("_")
		INTERAX.write(elements[0] + "\t" + elements[1] + "\t" + elements[2] + "\t")
		for alteration in ["loss", "gain"]:
			for protein in InteraX[candidate][alteration]:
				INTERAX.write(protein + ";")
			INTERAX.write("\t")

		INTERAX.write("\n")
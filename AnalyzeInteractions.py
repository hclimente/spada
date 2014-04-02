#!/soft/devel/python-2.7/bin/python

from include.libsmartas import *
import include.custom_iLoops_xml_parser as parser
from os.path import isfile
import sys
from math import fabs
import networkx as nx
import matplotlib.pyplot as plt

out = "Results/" + sys.argv[1] + "/"
myParser = parser.iLoopsParser()

minDiff = 6

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
		gene = elements[0]
		normalIso = elements[2]
		tumorIso = elements[3]
		interactions = {}
		InteraX = {}
		top += 1

		if not isfile(out + "iLoops/Output/" + normalIso + ".ips") and not isfile(out + "iLoops/Output/" + tumorIso + ".ips"):
			print("\t * Gene " + elements[1] + ": iLoops couldn't process this candidate.")
			continue

		with open(out + "InteraxChanges_" + gene + ".tsv", "w") as INTERAX:

			INTERAX.write("Transcript\tPartner\tPartner gene\tmaxRC\tDiff\n")
		
			print("\t * Gene " + elements[1])
					
			for iso, ori in zip([normalIso, tumorIso], ["Normal","Tumor"]):

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
					
					INTERAX.write(iso + "\t" + partner + "\t" + expressedTranscripts[partner] + "\t" + str(interactions[ori][partner]) + "\t" + str(InteraX[partner]) + "\n")

			G=nx.Graph()

	 		# DOTFile.write("graph " + gene + " {\n")
	 		# DOTFile.write("\t" + gene + " [label=" + gene + ", shape=polygon,sides=5];\n")
	 		G.add_node(gene)
	 		for partner in InteraX.keys():
	 		 	if fabs(InteraX[partner]) < minDiff:
	 		 		continue

	 		 	partnerCoreName = partner.split(".")[0]
	 			G.add_node(partnerCoreName)
	 			G.add_edge(gene, partnerCoreName)
	 			G[gene][partnerCoreName]['weight'] = str( fabs(InteraX[partner]) )
	 			G[gene][partnerCoreName]['label'] = str(InteraX[partner])
	 			
	 			if InteraX[partner] > 0:
	 				G[gene][partnerCoreName]['color']='red'
	 			elif InteraX[partner] < 0:
	 				G[gene][partnerCoreName]['color']='blue'

			nx.draw(G)
			plt.savefig( out + gene + ".png", dpi=1000)
			nx.write_dot(G, out + gene + '.dot')
			nx.write_graphml(G, out + gene + '.gml')

			if top >= 30:
				break
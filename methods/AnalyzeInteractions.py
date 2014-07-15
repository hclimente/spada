#!/soft/devel/python-2.7/bin/python

from libs.utils import *
from interface import iLoops_parser as parser
from os.path import isfile
import sys
from time import sleep
import tarfile
from os import remove

out = "Results/" + sys.argv[1] + "/"
inputType = sys.argv[2]
iLoopsVersion = sys.argv[3]
minReplicates = float(sys.argv[4])
myParser = parser.iLoopsParser()

InteraX = {}
expressedTranscripts = {}
articleCompilation = {}
loopFamilies = {}
candidatesGaudi = {}

with open(out + "expressedGenes.lst", "r") as EXPRESSED:
	for line in EXPRESSED:
		elements = line.strip().split("\t")
		expressedTranscripts[elements[0]] = elements[1]

with open("Data/Databases/compilationTable.tsv", "r") as compilationTable:
	for line in compilationTable:
		elements = line.strip().split("\t")
		genename = elements[0].split("|")[0]
		driverDb = [elements[8], elements[13]]
		apoptosis = [elements[10]]
		embryo = [elements[9]]
		rbps = [elements[3], elements[4], elements[5], elements[11]]
		epigenet = [elements[6]]
		
		if "yes" in driverDb:
			articleCompilation[genename] = "Driver"
		elif "yes" in rbps:
			articleCompilation[genename] = "RBP"
		elif "yes" in epigenet:
			articleCompilation[genename] = "Epigenetic factor"
		else:
			articleCompilation[genename] = "-"

with open("Data/TCGA/UnifiedFasta_" + iLoopsVersion + "_loopFamilies.txt", "r") as LOOP_FAMILIES:
	currentLoopFamily = ""
	for line in LOOP_FAMILIES:
		if ">" in line:
			currentLoopFamily = line[1:].strip().split("\t")[0]
			loopFamilies[currentLoopFamily] = set()
			loopFamilies[currentLoopFamily].add(line[1:].strip().split("\t")[1])
		else:
			loopFamilies[currentLoopFamily].add(line.strip())

with open( out + "candidatesGaudi.lst", "r") as GAUDI:
	for line in GAUDI:
		element = line.strip().split("\t")
		if int(element[1]) >= 0: candidatesGaudi[element[0]] = int(element[1])

with open(out + "candidateList.top.tsv", "r") as CANDIDATES:

	CANDIDATES.readline()
	
	for line in CANDIDATES:

		elements = line.split("\t")
		gene = elements[0]
		geneEntrez = elements[1]
		normalIso = elements[2]
		tumorIso = elements[3]
		reps = int(elements[4])
		uniprotIds = elements[6].split("#")
		tag = gene + "_" + normalIso + "_" + tumorIso
		interactions = {}

		if reps < minReplicates:
			break

		sys.stdout.write("\t * Gene " + elements[1])

		if not normalIso in candidatesGaudi.keys() and not tumorIso in candidatesGaudi.keys():
			print(" - not processed")
			continue
		elif not normalIso in candidatesGaudi.keys() or not tumorIso in candidatesGaudi.keys():
			print(" - not fully processed")
			continue	
		
		print("")
					
		for iso, ori in zip([normalIso, tumorIso], ["Normal","Tumor"]):
			tarFile = "iLoops/" + inputType + "/" + iLoopsVersion + "/" + iso + ".tar.gz"
			if candidatesGaudi[iso] == 2:
				for family in loopFamilies.keys():
					if iso in loopFamilies[family]:
						for relative in loopFamilies[family]:
							tarFile = "iLoops/" + inputType + "/" + iLoopsVersion + "/" + relative + ".tar.gz"
							if isfile(tarFile):
								break

			sys.stdout.write("\t\t* " + ori + ": " + iso)
			sys.stdout.flush()
			while not isfile(tarFile):
				sleep(900)
			
			tar = tarfile.open(tarFile)
			xmlFile = tar.getmembers()[0].name
			tar.extract(xmlFile)

			interactions[ori] = myParser.parseInteractions(
															thisCandidate				  = iso,
															xmlOutput					  = xmlFile,
															expressedIsoforms			  = expressedTranscripts.keys(),
															output_proteins               = False, 
															output_alignments             = False,
															output_domain_mappings        = False,
															output_protein_features       = False,
															output_domain_assignations    = False,
															output_interactions           = True,
															output_interaction_signatures = False,
															output_RF_results             = True,
															output_RF_precisions          = True
														  )	
			
			print(" - OK")
			remove(xmlFile)

		with open(out + "iLoops/InteraXChanges_" + tag + ".tsv", "w") as INTERAX:

			INTERAX.write("Partner\tPartner_gene\tRC_N\tRC_T\tdRC\tAnnotation\n")

			#Iterate all interactions that are predicted in normal
			for familyRep in interactions["Normal"].keys():
				diff = 100
				if familyRep in interactions["Tumor"].keys():
					diff = interactions["Normal"][familyRep] - interactions["Tumor"][familyRep]
				for family in [f for f in loopFamilies.keys() if familyRep in loopFamilies[f] ]:
					for partner in [p for p in loopFamilies[family] if p in expressedTranscripts.keys() ]:
						INTERAX.write( partner + "\t" + expressedTranscripts[partner] + "\t")
						INTERAX.write( str(interactions["Normal"][familyRep] ) + "\t")
						INTERAX.write( str(interactions["Tumor"].get(familyRep, "") ) + "\t" )
						INTERAX.write( str(diff) + "\t")
						INTERAX.write( articleCompilation.get(expressedTranscripts[partner], "-") + "\n")

			#Interate the rest i.e. the iterations in tumor that are not in normal
			for familyRep in [p for p in interactions["Tumor"].keys() if p not in interactions["Normal"].keys()]:
				for family in [f for f in loopFamilies.keys() if familyRep in loopFamilies[f] ]:
					for partner in [p for p in loopFamilies[family] if p in expressedTranscripts.keys() ]:
						INTERAX.write( partner + "\t" + expressedTranscripts[partner] + "\t")
						INTERAX.write( "\t") #Empty value for normal
						INTERAX.write( str(interactions["Tumor"][familyRep] ) + "\t" )
						INTERAX.write( "-100\t")
						INTERAX.write( articleCompilation.get(expressedTranscripts[partner], "-") + "\n")

		cmd("mkdir", out + "iLoops/" + tag )
		cmd("Pipeline/AnalyzeInteractors.r", out, tag )
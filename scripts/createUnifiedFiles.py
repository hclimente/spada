#!/soft/devel/python-2.7/bin/python

from os import listdir
from libsmartas import *
from fnmatch import filter
import custom_iLoops_xml_parser as parser

def writeFasta(basename, expressedTranscripts):
	wannaWrite = False
	fileCounter = 1
	transcriptCounter = 0

	expressedTranscriptsSet = set()
	if isinstance(expressedTranscripts, basestring):
		with open(expressedTranscripts, "r") as EXPRESSED:
			for line in EXPRESSED:
				expressedTranscriptsSet.add(line.strip())
	elif isinstance(expressedTranscripts, set):
		expressedTranscriptsSet = expressedTranscripts

	MULTIFASTA = open(basename + "_" + str(fileCounter) + ".fasta", "w")

	with open("UnifiedFasta.fa", "r") as gcMULTIFASTA:
		for line in gcMULTIFASTA:
			if ">" in line:

				if transcriptCounter >= 9999:
					fileCounter += 1
					MULTIFASTA.close()
					MULTIFASTA = open(basename + "_" + str(fileCounter) + ".fasta", "w")
					transcriptCounter = 0

				identifiers = line[1:].split("#")[0]
				
				if identifiers in expressedTranscriptsSet:
					MULTIFASTA.write(">" + identifiers + "\n")
					wannaWrite = True
					transcriptCounter += 1
				else:
					wannaWrite = False
			elif wannaWrite:
				MULTIFASTA.write(line)

def parseMapping(basename):
	loopFamilies = {}
	noLoops = set()
	loopModels = set()

	myParser = parser.iLoopsParser()

	for mappingBatch in filter(listdir("out/"), basename + "_*" ):
		for tens in ["", range(3) ]:
			for units in range(10):
				number = str(tens) + str(units)
				if number == "00": continue
				
				for mappingBatch_nodeResult in filter(listdir("out/" + mappingBatch + "/sge_output"), "*.assignation." + number + ".xml" ):
					xmlFile = "out/" + mappingBatch + "/sge_output/" + mappingBatch_nodeResult

					newLoops = myParser.parseLoops(
													xmlOutput					  = xmlFile,
													output_proteins               = True, 
													output_alignments             = False,
													output_domain_mappings        = False,
													output_protein_features       = True,
													output_domain_assignations    = False,
													output_interactions           = False,
													output_interaction_signatures = False,
													output_RF_results             = False,
													output_RF_precisions          = False
												  )
					
					for aTranscript, loops in sorted(newLoops.iteritems()):
						if loops not in loopFamilies.keys():
							loopFamilies[loops] = []

						loopFamilies[loops].append(aTranscript)
		
		for tens in ["", range(3) ]:
			for units in range(10):
				number = str(tens) + str(units)
				if number == "00": continue

				for mappingBatch_nodeErr in filter(listdir("out/" + mappingBatch + "/sge_output"), "*.assignation." + number + ".err.xml" ):
					errFile = "out/" + mappingBatch + "/sge_output/" + mappingBatch_nodeErr

					newNoLoops = myParser.parseNoLoops(
														xmlOutput                     = errFile,
														iLoopsPath                    = "/sbi/users/hectorc/AllMapping",
														output_proteins               = True, 
														output_alignments             = False,
														output_domain_mappings        = False,
														output_protein_features       = False,
														output_domain_assignations    = False,
														output_interactions           = False,
														output_interaction_signatures = False,
														output_RF_results             = False,
														output_RF_precisions          = False
													  )
					
					for aTranscript in newNoLoops:
						noLoops.add(aTranscript)

	with open("allMapping_loopFamilies.txt", "w") as loopFamiliesList:
		 for loopPattern in loopFamilies.keys():
		 	loopModels.add(loopFamilies[loopPattern][0])

		 	loopFamiliesList.write(">" + loopPattern + "\t" + loopFamilies[loopPattern][0] + "\n")
		 	for transcript in loopFamilies[loopPattern][1:]:
		 		loopFamiliesList.write(transcript + "\n")
	
	return (loopFamilies, loopModels, noLoops)

iLoops_ver = "iLoops_devel"
basename = "allMapping_" + iLoops_ver

writeFasta(basename, "transcripts.lst")
for expressedFasta in filter(listdir("."), basename + "_*.fasta"):

	assignationBatch = expressedFasta.split(".")[0].split("_")[1]

	cmd("/soft/devel/python-2.7/bin/python /sbi/programs/" + iLoops_ver + "/iLoops.py",
		"-f " + expressedFasta,
		"-j " + "out/" + basename + "_" + assignationBatch,
		"-x " + basename + "_" + assignationBatch + ".xml",
		"-g all",
		"-v",
		"-m",
		"-n 25",
		"-Q bigmem",
		"2>&1 >logs/" + basename + "_" + assignationBatch + ".log"
	   )

loopFamilies, loopModels, noLoops = parseMapping(basename)
writeFasta( basename + "_uniqLoops", loopModels)
with open( basename + "_noLoops.li", "a") as NOLOOPS:
	for aTranscript in noLoops:
		NOLOOPS.write(aTranscript + "\n")

fams = {}

with open( basename + "_loopFamilies.txt", "r") as LOOP_FAMILIES:
	currentK = ""
	for line in LOOP_FAMILIES:
		if ">" in line:
			elements=line.strip().split("\t")
			currentK = elements[0][1:]
			fams[currentK] = [elements[1]]
		else:
			fams[currentK].append(line.strip())


with open("UnifiedFasta.fa", "r") as gcMULTIFASTA, \
	 open("UnifiedFasta_" + iLoops_ver + ".fa", "w") as LOOP_FAMILIES:

	for line in gcMULTIFASTA:
		if not ">" in line:
			LOOP_FAMILIES.write(line)
			continue

		LOOP_FAMILIES.write(line.strip() + "#")
		ucscId = line.split("#")[0][1:]
		for loop_key, proteins in fams.iteritems():
			if ucscId in proteins:
				LOOP_FAMILIES.write(loop_key)

		LOOP_FAMILIES.write("\n")
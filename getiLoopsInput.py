#!/soft/devel/python-2.7/bin/python

import sys
from shutil import copy
from sh import *

def getGFFTrack(candidate, GFF3_TRACK, GFF2n_TRACK, GFF2t_TRACK):

	gffTrack = {}
	gffTrack["transcript"] = {}
	gffTrack["exon"] = {}

	ensGene = candidate.pop().strip()
	normalIso = candidate[0].strip()
	tumorIso = candidate[1].strip()
	geneName = ""

	gffInfo = getFeature("feature", ensGene, "feature=gene")
	if not gffInfo:
		return False

	for line in gffInfo:
		if line["ID"] == ensGene:
			geneName = line["external_name"]
			
			strand = "+"
			if line["strand"] == -1:
				strand = "-"
			GFF3_TRACK.write("##sequence-region   " + geneName + " " + str(line["start"])  + " " + str(line["end"]) + "\n")
			GFF3_TRACK.write(geneName + "\t.\t" + line["feature_type"] + "\t" + str(line["start"]) + "\t" +\
							 str(line["end"]) + "\t.\t" + strand + "\t.\t" + "ID=" + line["ID"] + "\n")

	for rawString in candidate:
		aCandidate = rawString.strip()	
		
		for feature in ["transcript", "exon"]:
			gffInfo = getFeature("feature", aCandidate, "feature=" + feature)
			if not gffInfo:
				return False
		
			for line in gffInfo:
				if feature == "exon" and line["Parent"] != aCandidate:
					continue
				elif feature == "transcript" and line["ID"] != aCandidate:
					continue

				if not line["ID"] in gffTrack[feature]:
					gffTrack[feature][line["ID"]] = {
						"seqid": "chr" + line["seq_region_name"],
						"type": line["feature_type"],
						"start": str(line["start"]),
						"end": str(line["end"]),
						"strand": "+",
						"Atributes": "ID=" + line["ID"] + ";Parent=" + line["Parent"],
						"Group": [aCandidate]
					}
	
					if line["strand"] == -1:
						gffTrack[feature][line["ID"]]["strand"] = "-"
				else:
					gffTrack[feature][line["ID"]]["Atributes"] += "," + line["Parent"]
					gffTrack[feature][line["ID"]]["Group"].append(line["Parent"])

	for feature in ["transcript", "exon"]:
		for geneId in gffTrack[feature]:
			thisLine = gffTrack[feature][geneId]
			GFF3_TRACK.write(geneName + "\t.\t" + thisLine["type"] + "\t" + thisLine["start"] + "\t" +\
						 thisLine["end"] + "\t.\t" + thisLine["strand"] + "\t.\t" + thisLine["Atributes"] + "\n")
			if feature == "exon":
				if normalIso in thisLine["Group"]:
					GFF2n_TRACK.write(thisLine["seqid"] + "\t.\t" + "exon" + "\t" + thisLine["start"] + "\t" +\
						  			  thisLine["end"] + "\t.\t" + thisLine["strand"] + "\t.\t" + "Transcript " + normalIso + "\n")
				if tumorIso in thisLine["Group"]:
					GFF2t_TRACK.write(thisLine["seqid"] + "\t.\t" + "exon" + "\t" + thisLine["start"] + "\t" +\
						  			  thisLine["end"] + "\t.\t" + thisLine["strand"] + "\t.\t" + "Transcript " + tumorIso + "\n")

	return True

expressedTranscripts = sys.argv[1];
candidateTranscripts = sys.argv[2];
getExpressedGenes = bool(int(sys.argv[3]))

if(getExpressedGenes):

	print("\t* Writing the multiFASTA file with all the expressed transcripts.")

	with open('Results/iLoops/ExpressedTranscripts.fasta', "w") as MULTIFASTA:
		with open(expressedTranscripts, "r") as EXPRESSED:
			for line in EXPRESSED:
				stableId = line.strip()
	
				sequence = (getFeature("sequence", stableId, "type=protein"))["seq"]
				if sequence:
					MULTIFASTA.write(">" + stableId + "\n")
					MULTIFASTA.write(sequence + "\n")

else:
	copy("old/iLoops/ExpressedTranscripts.fasta", "Results/iLoops/ExpressedTranscripts.fasta")

print("\t* Writing the pairs files.")
GFF3_TRACK = open('Results/candidates.v3.gff', 'w')
GFF2n_TRACK = open('Results/candidates_normal.v2.gff', 'w')
GFF2t_TRACK = open('Results/candidates_tumor.v2.gff', 'w')
GFF3_TRACK.write("##gff-version 3" + "\n")

with open(candidateTranscripts, "r") as CANDIDATES:
	for line in CANDIDATES:
		candidates = line.split("\t")
		
		if getGFFTrack(candidates, GFF3_TRACK, GFF2n_TRACK, GFF2t_TRACK):
			for rawCandidate in candidates:
	
				aCandidate = rawCandidate.strip()
				cmd("mkdir Results/iLoops/Input/" + aCandidate)
				fileNumber = 1
				numberOfCandidates = 0
	
				with open(expressedTranscripts, "r") as EXPRESSED:
					PAIRS = open("Results/iLoops/Input/" + aCandidate + '/' + aCandidate + '_' + str(fileNumber) + '.net', "w")
	
					for rawExpressed in EXPRESSED:
						expressedTranscript = rawExpressed.strip()
						PAIRS.write(aCandidate + "\t" + expressedTranscript + "\n")
						numberOfCandidates += 1
	
						if(numberOfCandidates >= 10000):
							fileNumber += 1
							numberOfCandidates = 0
							PAIRS.close()
							PAIRS = open("Results/iLoops/Input/" + aCandidate + '/' + aCandidate + '_' + str(fileNumber) + '.net', "w")
	
					PAIRS.close()
		else:
			for rawCandidate in candidates:
				aCandidate = rawCandidate.strip()
				cmd("rm -r Results/iLoops/Input/" + aCandidate)

GFF3_TRACK.close()
GFF2n_TRACK.close()
GFF2t_TRACK.close()

#!/soft/devel/python-2.7/bin/python

import httplib2, sys
from shutil import copy
from sh import *
from time import sleep
import json

def getGFF3Track(twoCandidates, GFF3_TRACK, GFF2n_TRACK, GFF2t_TRACK):

	gffTrack = {}
	gffTrack["gene"] = {}
	gffTrack["transcript"] = {}
	gffTrack["exon"] = {}

	ensGene = ""
	geneName = ""

	for rawCandidate in twoCandidates:
		aCandidate = rawCandidate.strip()	
		for feature in ["transcript", "exon", "gene"]:
			
			ext = "/feature/id/" + aCandidate + "?feature=" + feature
			resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
			
			if not resp.status == 200:
				print("\t\tCandidate", aCandidate, "(", server + ext, ") not valid. Error", resp.status)
				return True, ""
		
			gffInfo = json.loads(content)
			
			for line in gffInfo:

				if feature == "transcript":
					if not ensGene and line["ID"] == aCandidate:
						ensGene = line["Parent"]
					
					if line["ID"] != aCandidate:
						continue
				elif feature == "gene":
					if line["ID"] != ensGene:
						continue
					geneName = line["external_name"]
				elif feature == "exon":
					if line["Parent"] != aCandidate:
						continue

				if not line["ID"] in gffTrack[feature]:
					gffTrack[feature][line["ID"]] = {}

					gffTrack[feature][line["ID"]]["seqid"] = "chr" + line["seq_region_name"]
					gffTrack[feature][line["ID"]]["type"] = line["feature_type"]
					gffTrack[feature][line["ID"]]["start"] = str(line["start"])
					gffTrack[feature][line["ID"]]["end"] = str(line["end"])
					gffTrack[feature][line["ID"]]["strand"] = "+"
					gffTrack[feature][line["ID"]]["Atributes"] = "ID=" + line["ID"]
	
					if line["strand"] == -1:
						gffTrack[feature][line["ID"]]["strand"] = "-"
	
					if feature == "transcript" or feature == "exon":
						gffTrack[feature][line["ID"]]["Atributes"] += ";Parent=" + line["Parent"]
					else:
						gffTrack[feature][line["ID"]]["Atributes"] += ";Name=" + line["external_name"]
				else:
					if feature == "transcript" or feature == "exon":
						gffTrack[feature][line["ID"]]["Atributes"] += "," + line["Parent"]

	GFF3_TRACK.write("##sequence-region   " + geneName + " " + gffTrack["gene"][ensGene]["start"]  + " " + gffTrack["gene"][ensGene]["end"] + "\n")

	for feature in ["gene", "transcript", "exon"]:
		for geneId in gffTrack[feature]:
			thisLine = gffTrack[feature][geneId]
			GFF3_TRACK.write(geneName + "\t.\t" + thisLine["type"] + "\t" + thisLine["start"] + "\t" +\
						 thisLine["end"] + "\t.\t" + thisLine["strand"] + "\t.\t" + thisLine["Atributes"] + "\n")

	for feature in ["transcript", "exon"]:
		for geneId in gffTrack[feature]:
			type = "exon"
			if feature == "transcript":
				type = "mRNA"
			thisLine = gffTrack[feature][geneId]
			GFF2n_TRACK.write(thisLine["seqid"] + "\t.\t" + type + "\t" + thisLine["start"] + "\t" +\
						   thisLine["end"] + "\t.\t" + thisLine["strand"] + "\t.\t" + thisLine["Atributes"] + "\n")

	return False

http = httplib2.Http(".cache")
server = "http://beta.rest.ensembl.org"

expressedTranscripts = sys.argv[1];
candidateTranscripts = sys.argv[2];
getExpressedGenes = bool(int(sys.argv[3]))

# if(getExpressedGenes):

# 	print("\t* Writing the multiFASTA file with all the expressed transcripts.")

# 	expressedMultiFasta = open('Results/iLoops/ExpressedTranscripts.fasta', "w")
# 	ensemblQueryRestriction = 0

# 	with open(expressedTranscripts, "r") as EXPRESSED:
# 		for line in EXPRESSED:
# 			stableId = line.strip()
# 			ext = "/sequence/id/" + stableId + "?type=protein"
# 			resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"text/plain"})

# 			if not resp.status == 200:
# 				print("\t\tCouldn't retrieve", stableId, "(", server + ext, "). Error", resp.status)
# 				continue

# 			expressedMultiFasta.write(">" + stableId + "\n")
# 			expressedMultiFasta.write(content + "\n")

# 			#Ensembl REST API doesn't accept more than 3 queries/second
# 			ensemblQueryRestriction += 1
# 			if ensemblQueryRestriction == 3:
# 				sleep(1)
# 				ensemblQueryRestriction = 0

# 	expressedMultiFasta.close()
# else:
# 	copy("old/iLoops/ExpressedTranscripts.fasta", "Results/iLoops/ExpressedTranscripts.fasta")

print("\t* Writing the pairs files.")
GFF3_TRACK = open('Results/candidates.v3.gff', 'w')
GFF2n_TRACK = open('Results/candidates_normal.v2.gff', 'w')
GFF2t_TRACK = open('Results/candidates_tumor.v2.gff', 'w')
GFF3_TRACK.write("##gff-version 3" + "\n")

with open(candidateTranscripts, "r") as CANDIDATES:
	for line in CANDIDATES:
		candidates = line.split("\t")
		delete = False

		delete = getGFF3Track(candidates, GFF3_TRACK, GFF2n_TRACK, GFF2t_TRACK)

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

		if delete:
			for rawCandidate in candidates:
				aCandidate = rawCandidate.strip()
				cmd("rm -r Results/iLoops/Input/" + aCandidate)

GFF_TRACK.close()
GFF2n_TRACK.close()
GFF2t_TRACK.close()
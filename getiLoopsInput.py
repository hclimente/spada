#!/soft/devel/python-2.7/bin/python

import sys
from shutil import copy
from sh import *
import os

expressedTranscripts = sys.argv[1];
candidateTranscripts = sys.argv[2];
getExpressedGenes = bool(int(sys.argv[3]))

if(getExpressedGenes):

	print("\t* Writing the multiFASTA file with all the expressed transcripts.")

	with open('Results/iLoops/ExpressedTranscripts.fasta', "w") as MULTIFASTA:
		with open(expressedTranscripts, "r") as EXPRESSED:
			for line in EXPRESSED:
				stableId = line.strip()
	
				proteinFeat = getFeature("sequence", stableId, "type=protein")
				if proteinFeat:
					MULTIFASTA.write(">" + stableId + "\n")
					MULTIFASTA.write(proteinFeat["seq"] + "\n")

else:
	copy("old/iLoops/ExpressedTranscripts.fasta", "Results/iLoops/ExpressedTranscripts.fasta")

print("\t* Writing the pairs files.")
GFF3_TRACK = open('Results/candidates.top.v3.gff', 'w')
GFF2n_TRACK = open('Results/candidates_normal.top.v2.gff', 'w')
GFF2t_TRACK = open('Results/candidates_tumor.top.v2.gff', 'w')
GFF3_TRACK.write("##gff-version 3" + "\n")

with open(candidateTranscripts, "r") as CANDIDATES:
	for line in CANDIDATES:
		elements = line.split("\t")
		candidates = [elements[2], elements[3]]
		
		if getGFFTrack(candidates + [elements[1]], GFF3_TRACK, GFF2n_TRACK, GFF2t_TRACK):
			for aCandidate in candidates:
	
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
			for aCandidate in candidates:
				delFile = "Results/iLoops/Input/" + aCandidate
				if os.path.exists(delFile):
					cmd("rm -r", delFile)

GFF3_TRACK.close()
GFF2n_TRACK.close()
GFF2t_TRACK.close()

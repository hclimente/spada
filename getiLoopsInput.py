#!/soft/devel/python-2.7/bin/python

import httplib2, sys
from shutil import copy
from sh import *
 
http = httplib2.Http(".cache")
server = "http://beta.rest.ensembl.org"

expressedTranscripts = sys.argv[1];
candidateTranscripts = sys.argv[2];
getExpressedGenes = sys.argv[3];

notFound = []

if(getExpressedGenes):

	print "\t* Writing the multiFASTA file with all the expressed transcripts."

	expressedMultiFasta = open('Results/iLoops/ExpressedTranscripts.fasta', "w")
	with open(expressedTranscripts, "r") as EXPRESSED:
		for line in EXPRESSED:
			stableId = line.strip()
			ext = "/sequence/id/" + stableId + "?type=protein"
			resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"text/plain"})
 
			if not resp.status == 200:
				print "\t\tCouldn't retrieve", stableId, "(", server + ext, "). Error", resp.status
				notFound.append(stableId)
				continue
 
			expressedMultiFasta.write(">" + stableId + "\n")
			expressedMultiFasta.write(content + "\n")

	expressedMultiFasta.close()
else:
	copy("old/iLoops/input/ExpressedTranscripts.fasta", "esults/iLoops/ExpressedTranscripts.fasta")

print "\t* Writing the pairs files."
GFF_TRACK = open('Results/candidates.gff', 'w')

with open(candidateTranscripts, "r") as CANDIDATES:
	for line in CANDIDATES:
		candidates = line.split("\t")
		delete = False

		for rawCandidate in candidates:

			aCandidate = rawCandidate.strip()
			ext = "/feature/id/" + aCandidate + "?feature=exon"
			resp, content = http.request(server+ext, method="GET", headers={"text/x-gff3"})

			if not resp.status == 200:
				print "\t\tCandidate ", stableId, " (", server + ext, ") not valid. Error ", resp.status
				delete = True
				break

			GFF_TRACK.write(content)

			cmd("mkdir Results/iLoops/input/" + aCandidate)
			fileNumber = 1
			numberOfCandidates = 0

			with open(expressedTranscripts, "r") as EXPRESSED:
				PAIRS = open("Results/iLoops/input/" + aCandidate + '/' + aCandidate + '_' + str(fileNumber) + '.net', "w")

				for expressedTranscript in EXPRESSED:
					PAIRS.write(aCandidate + "\t" + expressedTranscript + "\n")
					numberOfCandidates += 1

					if(numberOfCandidates >= 10000):
						fileNumber += 1
						numberOfCandidates = 0
						PAIRS.close()
						PAIRS = open("Results/iLoops/input/" + aCandidate + '/' + aCandidate + '_' + str(fileNumber) + '.net', "w")

				PAIRS.close()

		if delete:
			for rawCandidate in candidates:
				aCandidate = rawCandidate.strip()
				cmd("rm -r Results/iLoops/input/" + aCandidate)

GFF_TRACK.close()
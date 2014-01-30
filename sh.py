#!/soft/devel/python-2.7/bin/python

from subprocess import call,Popen,PIPE
from rpy2.robjects import r
import httplib2
import json
from time import sleep
import urllib2

def cmd(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)
	call(command, shell=True)

def cmdOut(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)
	
	return Popen(command, shell=True, stdout=PIPE).stdout.read().strip()

def setRWorkspace(wd, Conditions, Compartments, Replicates, Kmer):

	r("wd <- \"" + wd + "\"")
	r("inputData <- list()")
	
	condStat = 'inputData[["Conditions"]] <- c('
	for condition in Conditions:
		condStat = condStat + "\"" + condition + "\","
	condStat = condStat[:-1] + ")"
	
	compStat = 'inputData[["Compartments"]] <- c('
	for compartment in Compartments:
		compStat = compStat + "\"" + compartment + "\","
	compStat = compStat[:-1] + ")"
	
	replStat = 'inputData[["Replicates"]] <- c('
	for replicate in Replicates:
		replStat = replStat + "\"" + replicate + "\","
	replStat = replStat[:-1] + ")"
	
	kmerStat = 'inputData[["K-mer"]] <- c('
	for kmer in Kmer:
		kmerStat = kmerStat + "\"" + kmer + "\","
	kmerStat = kmerStat[:-1] + ")"
	
	r(condStat)
	r(compStat)
	r(replStat)
	r(kmerStat)
	
	r('save.image("SmartAS.RData")')

def setEnvironment(wd, initialStep, Conditions, Compartments, Replicates, Kmer):

	print("* Preparing the environment")
	cmd("cd " + wd)
	cmd("rm -r old2; mv old old2")
	cmd("mv Results old; mv SmartAS.RData old/SmartAS.old.RData")
	cmd("mkdir -p Results/iLoops/Output/Mapping")
	cmd("mkdir Results/iLoops/Input")
	cmd("mkdir Results/RWorkspaces")
	cmd("mkdir Results/DataExploration")

	if initialStep <= 1:
		setRWorkspace(wd, Conditions, Compartments, Replicates, Kmer)
		getDB()
	if initialStep > 1:
		cmd("cp -r old/DataExploration Results")
		cmd("cp -r old/RWorkspaces/1_ExploreData.RData Results/RWorkspaces")
		cmd("cp -r old/RWorkspaces/1_ExploreData.RData SmartAS.RData")
		cmd("cp -r old/10C1_30.tsv old/7C1_30.tsv old/10C2_30.tsv old/7C2_30.tsv old/IntraReplicateC1_30.tsv old/IntraReplicateC2_30.tsv Results")
		
		#If any kind of data is recycled, check that the parameters didn't change between runs.
		diff = cmdOut('diff Results/Parameters.cfg old/Parameters.cfg 2>&1')
	
		if diff:
			print("WARNING: parameter files don't match.")
			with open("Results/Parameters.cfg", "w") as paramFile:
				paramFile.write("\nDID NOT MATCH")
	if initialStep > 2:
		cmd("cp -r old/RWorkspaces/2_GetCandidates.RData Results/RWorkspaces")
		cmd("cp -r old/RWorkspaces/2_GetCandidates.RData SmartAS.RData")
		cmd("cp old/candidateList.lst old/candidateList_withGenenames.tsv old/expressedGenes.lst Results")
		cmd("cp old/candidates.v3.gff old/candidates_normal.v2.gff old/candidates_tumor.v2.gff Results")
	if initialStep > 3:
		cmd("cp old/candidateInteractions.tsv old/candidateList.top.lst Results")
	if initialStep > 4:
		cmd("cp -r old/iLoops/Input Results/iLoops/")
		cmd("cp -r old/iLoops/ExpressedTranscripts.fasta Results/iLoops/")
		cmd("cp old/candidates.top.v3.gff old/candidates_normal.top.v2.gff old/candidates_tumor.top.v2.gff Results")
	if initialStep > 5:
		cmd("cp -r old/iLoops/Output Results/iLoops")

def getDB():
	with open("Data/Intogen.tsv", "w") as Intogen:
	
		query="""
		DEFINE
			intogen='/data/project/gene',
			genes='https://bitbucket.org/intogen/intogen-sources.git?ensembl/hsa/genes',
			projects='/data/projects'
		ON
			'https://bitbucket.org/intogen/intogen-mutations.git'
		SELECT
			genes (GENE_ID, SYMBOL),
			projects (PROJECT_NAME),
			intogen (SAMPLE_PROP, FM_QVALUE, CLUST_QVALUE)
		FROM
			intogen
		WHERE
			(
				intogen.FM_QVALUE < '0.05'
				OR
				intogen.CLUST_QVALUE < '0.05'
			)
			"""
			
		req = urllib2.Request("http://www.intogen.org/oql")
		res = urllib2.urlopen(req, query)
		Intogen.write(res.read())

def waitPID(pidQueue):
	for job in pidQueue:
		while True:
			if cmdOut("ps --pid", job, " | grep -v", job, "| wc -l") == "0":
				break
			else:
				print("Awaiting for completion of iLoops jobs.")
				sleep(900)

def getFeature(query, ensemblId, question):

	http = httplib2.Http(".cache")
	server = "http://beta.rest.ensembl.org"

	ext = "/" + query + "/id/" + ensemblId + "?" + question
	resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
			
	if not resp.status == 200:
		#errorCode = {400: "Bad Request (id not found)", 404: "Not Found (incorrect format)", 429: "Too Many Requests", 503: "Service Unavailable"}
		#print("\t\tError with " + ensemblId + " (" + server + ext + "): " + errorCode[resp.status])
		print("\t\tError with " + ensemblId + " (" + server + ext + "): " + str(resp.status))
		return {}
	
	#Ensembl REST API doesn't accept more than 3 queries/second.
	sleep(0.33)

	return json.loads(content)

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

def printParam(initialStep, wd, gaudiWd, minExpression, minCandidateExpression, minPSI, Conditions, Compartments, Replicates, Kmer, top):
	with open("Results/Parameters.cfg", "w") as paramFile:
		paramFile.write("initialStep=" + str(initialStep) + "\n")
		paramFile.write("wd=" + wd + "\n")
		paramFile.write("gaudiWd=" + gaudiWd + "\n")
		paramFile.write("minExpression=" + str(minExpression) + "\n")
		paramFile.write("minCandidateExpression=" + str(minCandidateExpression) + "\n")
		paramFile.write("minPSI=" + str(minPSI) + "\n")
		paramFile.write("Conditions=" + str(Conditions) + "\n")
		paramFile.write("Compartments=" + str(Compartments) + "\n")
		paramFile.write("Replicates=" + str(Replicates) + "\n")
		paramFile.write("Kmer=" + str(Kmer) + "\n")
		paramFile.write("top=" + str(top) + "\n")
#!/soft/devel/python-2.7/bin/python

import sys

def doStuff(knsur, iLoopsVersion):

	minReplicates = 0
	out = "testResults/TCGA/" + knsur + "_mE-1.0/"

	with open(out + "Parameters.cfg", "r") as PARAM:
		for line in PARAM:
			if "Replicates" in line and "unpaired" not in line:
				minReplicates = float(line.strip().split("=")[1]) * 0.1

	with open(out + "candidateList.top.tsv", "r") as CANDIDATES, open(knsur + ".stats.tsv", "w") as OUTPUT:
		OUTPUT.write("Gene\tnT\ttT\tDriver\tCDS\tmapLoops\tdiffLoops\n")
		CANDIDATES.readline()
		for line in CANDIDATES:
			elements = line.strip().split("\t")
			hugo = elements[0]
			gn = elements[1]
			nTx = elements[2]
			tTx = elements[3]
			patients = int(elements[4])
			knownPPI = 0
			if elements[5] != "na":
				knownPPI = int(elements[5])
			uniprot = elements[6].split("#")
			CDS = False
			if elements[7] == "Yes":
				CDS = True
			UTR_change = False
			if elements[8] == "Yes": 
				UTR_change = True
			CDS_change = False
			if elements[9] == "Yes":
				CDS_change = True
			driver = False
			if "yes" in [elements[10], elements[16], elements[17], elements[21]]:
				driver = True

			if patients < minReplicates:
				break

			#loops, loopDiff = checkLoopDiff(nTx, tTx, iLoopsVersion)
			#OUTPUT.write(gn + "\t" + nTx + "\t" + tTx + "\t" + str(driver) + "\t" + str(CDS) + "\t" + str(loops) + "\t" + str(loopDiff) + "\n" )

			OUTPUT.write(gn + "\t" + nTx + "\t" + tTx + "\t" + str(driver) + "\n" )

def checkLoopDiff(nTx, tTx,iLoopsVersion):

	nLoops = ""
	tLoops = ""

	loops = False
	loopDiff = False

	with open("Data/TCGA/UnifiedFasta_" + iLoopsVersion + ".fa", "r") as FASTA:
		for line in FASTA:
			if nTx in line:
				nLoops = line.strip().split("#")[3]
			elif tTx in line:
				tLoops = line.strip().split("#")[3]

			if nLoops and tLoops:
				break

	if nLoops and tLoops:
		loops = True
		if nLoops != tLoops:
			loopDiff = True

	return (loops, loopDiff)

knsur = "luad"

doStuff(knsur, "iLoops_devel")
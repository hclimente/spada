#!/soft/devel/python-2.7/bin/python

import sys
from subprocess import call
from os import listdir
from fnmatch import filter

def cmd(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)

	call(command, shell=True)

inputType = sys.argv[1]
tag1 = sys.argv[2]

reps = []

outDir = "Data/Input/" + inputType + "/" + tag1 + "/"

cmd("rm -r", outDir, "; mkdir -p", outDir)

if(inputType == "GENCODE"):
	Conditions = {"10": "N", "7": "T"}
	for condition in Conditions.keys():
		replicateCounter = 1
		for sample in filter(listdir("Data/GENCODE/Rawdata/" + tag1 + "-kmer-length/"), condition + "C*_*"):
			with open("Data/GENCODE/Rawdata/" + tag1 + "-kmer-length/" + sample + "/quant_bias_corrected.sf", "r") as FILE, \
				 open(outDir + str(replicateCounter) + "_" + Conditions[condition] + ".tsv", "w") as FILTERED:
				for line in FILE:
					if line.find("#") == -1:
						tableValues=line.strip().split("\t")
						splitIds=tableValues[0].split("|")
						FILTERED.write(splitIds[1] + "\t" + splitIds[0] + "\t" + splitIds[5] + "\t" + tableValues[2] + "\n")
			reps.append(replicateCounter)
			replicateCounter += 1

elif(inputType == "TCGA"):
	patients = []
	with open("Data/TCGA/Rawdata/" + tag1 + "_iso_tpm_paired-filtered.txt", "r") as FILE:
		firstLine = FILE.readline().strip().split("\t")
		for sampleType in ["N", "T"]:
			replicateCounter = 1
			for patient in range(0, len(firstLine)/2):
				patients.append(open(outDir + str(replicateCounter) + "_" + sampleType + ".tsv", "w"))
				reps.append(replicateCounter)
				replicateCounter += 1
		
		for line in FILE:
						
			splitted = line.strip().split("\t")
			seqTags = splitted[0].split(",")
			
			gene = seqTags[0]
			transcript = seqTags[1]
			name = gene.split("|")[0]
			
			currentCol = 1
			for patient in patients:
				patient.write( gene + "\t" + transcript + "\t" + name + "\t" + splitted[currentCol] + "\n")
				currentCol += 1
	for patient in patients:
		patient.close()
elif(inputType == "TCGA_unpaired"):
	patients = []
	with open("Data/TCGA/Rawdata/" + tag1 + "_iso_tpm_tumor-filtered.txt", "r") as FILE:
		firstLine = FILE.readline().strip().split("\t")
		replicateCounter = 1
		for patient in range(0, len(firstLine) ):
			patients.append(open(outDir + str(replicateCounter) + "_T.tsv", "w"))
			reps.append(replicateCounter)
			replicateCounter += 1
		
		for line in FILE:
						
			splitted = line.strip().split("\t")
			seqTags = splitted[0].split(",")
			
			gene = seqTags[0]
			transcript = seqTags[1]
			name = gene.split("|")[0]
			
			currentCol = 1
			for patient in patients:
				patient.write( gene + "\t" + transcript + "\t" + name + "\t" + splitted[currentCol] + "\n")
				currentCol += 1
	for patient in patients:
		patient.close()

cmd("cp Data/config.cfg", outDir)
with open(outDir + "config.cfg", "a") as CONFIG:
	CONFIG.write("\ntag1=" + tag1)
	CONFIG.write("\nReplicates=" + str(max(reps)))
	CONFIG.write("\ninputType=" + inputType)
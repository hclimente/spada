#!/soft/devel/python-2.7/bin/python

from libs import options
from libs import utils

import sys
import os
import fnmatch

def standarizeInput():
	numberOfPatients = 0
	outDir = "Data/Input/{0}/{1}/".format(options.Options().inputType, options.Options().tag)

	if options.Options().inputType == "GENCODE":
		utils.cmd("rm -r", outDir, "; mkdir -p", outDir)
		Conditions = {"10": "N", "7": "T"}
		for condition in Conditions.keys():
			replicateCounter = 1
			for sample in fnmatch.filter(os.listdir("Data/GENCODE/Rawdata/" + options.Options().tag + "-kmer-length/"), condition + "C*_*"):
				with open("{0}{1}_{2}.tsv".format(outDir,replicateCounter,Conditions[condition]), "w") as FILTERED:
					for line in utils.readTable("Data/GENCODE/Rawdata/{0}-kmer-length/{1}/quant_bias_corrected.sf".format(options.Options().tag, sample)):
						if "#" not in line[0]:
							FILTERED.write("{0}\t{1}\n".format(line[0], line[2]))
				if replicateCounter > numberOfPatients: numberOfPatients = replicateCounter
				replicateCounter += 1


	elif options.Options().inputType == "TCGA":
		utils.cmd("rm -r", outDir, "; mkdir -p", outDir)
		patientFiles = []

		inputFile = "Data/TCGA/Rawdata/{0}_iso_tpm_paired-filtered.txt".format(options.Options().tag)

		if options.Options().unpairedReplicates:
			inputFile = "Data/TCGA/Rawdata/{0}_iso_tpm_tumor-filtered.txt".format(options.Options().tag)

		for line in utils.readTable(inputFile):
			numberOfPatients = (len(line) - 1)/2
			
			if options.Options().unpairedReplicates:
				for patient in range(1, numberOfPatients+1):
					patientFiles.append(open("{0}{1}_T.tsv".format(outDir, patient), "w"))
			else:
				for sampleType in ["N", "T"]:
					for patient in range(1, numberOfPatients+1):
						patientFiles.append(open("{0}{1}_{2}.tsv".format(outDir, patient, sampleType), "w"))
			break

		for line in utils.readTable(inputFile):
			identifiers = line[0].split(",")
				
			gene = identifiers[0]
			transcript = identifiers[1]
				
			currentCol = 1
			for patient in patientFiles:
				patient.write( "{0}\t{1}\t{2}\n".format(gene, transcript, line[currentCol]) )
				currentCol += 1
		
		for patient in patientFiles:
			patient.close()

	options.Options().printToFile(initialStep=1,replicates=numberOfPatients)
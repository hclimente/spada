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
		patientFiles 		 = []
		unpairedPatientFiles = []
		pairedPatients 	 	 = set()
		unpairedPatients 	 = set()
		tag = options.Options().tag 
		if options.Options().unpairedReplicates:
			tag = tag[2:]

		inputFile = "Data/TCGA/Rawdata/{0}_iso_tpm_paired-filtered.txt".format(tag)

		for line in utils.readTable(inputFile,header=False):
			pairedPatients = set([ x[:-1] for x in line ])
			patientFiles = [ open("{0}{1}_{2}.tsv".format(outDir,line[p][:-1],s),"w") for s in ["N", "T"] for p in range(len(pairedPatients)) ]
			break

		for line in utils.readTable(inputFile):
			identifiers = line[0].split(",")
			gene = identifiers[0]
			transcript = identifiers[1]

			for i in range(len(patientFiles)):
				patientFiles[i].write( "{0}\t{1}\t{2}\n".format(gene, transcript, line[i+1]) )
		
		for patient in patientFiles:
			patient.close()

		if options.Options().unpairedReplicates:
			inputFile = "Data/TCGA/Rawdata/{0}_iso_tpm_tumor-filtered.txt".format(tag)

			for line in utils.readTable(inputFile,header=False):
				unpairedPatients = set([ x[:-1] for x in line ])
				unpairedPatientFiles = [ open("{0}{1}_T.tsv".format(outDir,line[p][:-1]),"w") for p in range(len(unpairedPatients)) ]
				break

			for line in utils.readTable(inputFile):
				identifiers = line[0].split(",")
				gene = identifiers[0]
				transcript = identifiers[1]

				for i in range(len(unpairedPatientFiles)):
					unpairedPatientFiles[i].write( "{0}\t{1}\t{2}\n".format(gene, transcript, line[i+1]) )

			for patient in unpairedPatientFiles:
				patient.close()

	options.Options().printToFile(initialStep=1,replicates=pairedPatients,unpairedReplicates=unpairedPatients)
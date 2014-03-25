#!/soft/devel/python-2.7/bin/python

from os import listdir
from libsmartas import cmd

for michalFile in listdir("."):
	core = michalFile.split(".")[0]
	kansurType = michalFile.split("_")[0]
	with open(michalFile, "r") as MICHAL, open("/home/hector/SmartAS/Data/TCGA/External/" + core + ".tsv","w") as SMARTAS:
		for line in MICHAL:
			elements = line.strip().split("\t")
			tumor = elements[0].split(",")[1]
			normal = elements[1].split(",")[1]
			gene = elements[0].split(",")[0]
			genename = gene.split("|")[0]

			SMARTAS.write(genename + "\t" + gene + "\t" + normal + "\t" + tumor + "\n")

	ans = raw_input("Background from " + kansurType + " for " + michalFile + "? (y/kansur type)")

	if ans != "y":
		kansurType = ans

	cmd("ln", "-s", "/home/hector/SmartAS/Results/TCGA/" + kansurType + "_mE-1.0/expressedGenes.lst", 
		"/home/hector/SmartAS/Data/TCGA/External/" + core + "_expressedGenes.lst"	)
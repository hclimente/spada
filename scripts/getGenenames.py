#!/soft/devel/python-2.7/bin/python

from sh import *
import sys

candidates = sys.argv[1];

with open(candidates, "r") as CANDIDATES:
	with open ("Results/candidateList_withGenenames.lst", "w") as GENENAMES:
		GENENAMES.write("Genename\tEnsembl Id\tNormal isoform\tTumor isoform\n")
		for line in CANDIDATES:
			ids = line.split("\t")
			ensGene = ids[2].strip()

			geneInfo = getFeature("feature", ensGene, "feature=gene")

			for line in geneInfo:
				if line["ID"] == ensGene:
					GENENAMES.write(line["external_name"] + "\t" + ensGene + "\t" + ids[0] + "\t" + ids[1] + "\n")
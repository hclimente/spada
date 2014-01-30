#!/soft/devel/python-2.7/bin/python

from sh import *
import sys

candidates = sys.argv[1];

GFF3_TRACK = open('Results/candidates.v3.gff', 'w')
GFF2n_TRACK = open('Results/candidates_normal.v2.gff', 'w')
GFF2t_TRACK = open('Results/candidates_tumor.v2.gff', 'w')
GFF3_TRACK.write("##gff-version 3" + "\n")

with open(candidates, "r") as CANDIDATES:
	with open ("Results/candidateList_withGenenames.tsv", "w") as GENENAMES:
		GENENAMES.write("Genename\tEnsembl Id\tNormal isoform\tTumor isoform\n")
		for line in CANDIDATES:
			ids = line.split("\t")
			ensGene = ids[2].strip()
			geneInfo = getFeature("feature", ensGene, "feature=gene")

			for line in geneInfo:
				if line["ID"] == ensGene:
					GENENAMES.write(line["external_name"] + "\t" + ensGene + "\t" + ids[0] + "\t" + ids[1] + "\n")
			
			getGFFTrack(ids, GFF3_TRACK, GFF2n_TRACK, GFF2t_TRACK)


GFF3_TRACK.close()
GFF2n_TRACK.close()
GFF2t_TRACK.close()
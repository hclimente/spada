#!/soft/devel/python-2.7/bin/python

from sh import *
import sys

candidates = sys.argv[1];

with open(candidates, "r") as CANDIDATES, open('Results/candidates_normal.v2.gff', 'w') as GFF2n_TRACK, open('Results/candidates_tumor.v2.gff', 'w') as GFF2t_TRACK:
	with open("Data/annotation.gtf", "r") as ALLTRANSCRIPTS:
		candTnt = []
		for line in CANDIDATES:
			ids = line.split("\t")
			pair = [ids[0], ids[1]]
			candTnt.append(pair)
		
		for line in ALLTRANSCRIPTS:
			for pair in candTnt:
				if line.find(pair[0]) != -1:
					GFF2n_TRACK.write(line)
				elif line.find(pair[1]) != -1:
					GFF2t_TRACK.write(line)


exit()

candidates = sys.argv[1];

GFF3_TRACK = open('Results/candidates.v3.gff', 'w')
GFF2n_TRACK = open('Results/candidates_normal.v2.gff', 'w')
GFF2t_TRACK = open('Results/candidates_tumor.v2.gff', 'w')
GFF3_TRACK.write("##gff-version 3" + "\n")

with open(candidates, "r") as CANDIDATES:
	with open ("Results/candidateList.tsv", "w") as GENENAMES:
		#GENENAMES.write("Genename\tEnsembl Id\tNormal isoform\tTumor isoform\n")
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

cmd("rm Results/candidateList.lst")
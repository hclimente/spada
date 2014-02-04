#!/soft/devel/python-2.7/bin/python

from sh import *
import sys

candidates = sys.argv[1];

with open(candidates, "r") as CANDIDATES, open('Results/candidates_normal.gff', 'w') as GFF2n_TRACK, open('Results/candidates_tumor.gff', 'w') as GFF2t_TRACK:
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
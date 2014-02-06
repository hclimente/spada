#!/soft/devel/python-2.7/bin/python

from sh import *
import sys

candidates = sys.argv[1];
particle = ""
if candidates.find("top") != -1:
	particle = ".top"

with open(candidates, "r") as CANDIDATES, \
	 open('Results/candidates_normal' + particle + '.gff', 'w') as GFF2n_TRACK, \
	 open('Results/candidates_tumor' + particle + '.gff', 'w') as GFF2t_TRACK, \
	 open("Data/annotation.gtf", "r") as ALLTRANSCRIPTS:
		candTnt = []
		for line in CANDIDATES:
			ids = line.strip().split("\t")
			candTnt.append([ids[2], ids[3]])
		
		for line in ALLTRANSCRIPTS:
			for pair in candTnt:
				if line.find(pair[0]) != -1:
					GFF2n_TRACK.write(line)
				elif line.find(pair[1]) != -1:
					GFF2t_TRACK.write(line)
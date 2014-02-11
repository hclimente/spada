#!/soft/devel/python-2.7/bin/python

from libsmartas import *
import sys

candidates = sys.argv[1]
out = sys.argv[2]

with open(candidates, "r") as CANDIDATES, \
	 open('Results/' + out + "/candidates_normal.gff", 'w') as GFF2n_TRACK, \
	 open("Results/" + out + "/candidates_tumor.gff", 'w') as GFF2t_TRACK, \
	 open("Data/GENCODE/annotation.gtf", "r") as ALLTRANSCRIPTS:
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
#!/soft/devel/python-2.7/bin/python

import sys, os
from sh import *
from time import sleep
from fnmatch import filter

if(len(sys.argv) != 1):
	print("No arguments passed.")
	exit()

os.chdir(sys.argv[1])

pidQueue = []

for transcript in filter(os.listdir("input"), "ENST*"):
	for configFile in filter(os.listdir("input/" + transcript), "*net"):
		#Map the loops for the query sequence.
		cmd("/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py",
			"-f input/ExpressedTranscripts.fasta"
			"-q input/" + transcript + "/" + configFile,
			"-j output/" + configFile,
			"-x " + configFile + ".xml",
			"-g all",
			"-n 25",
			"-Q sbi",
			"-c 1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50",
			"-v",
			"-m &"
			)
		
	 	pidQueue.append(cmdOut("$!"))

waitPID(pidQueue)
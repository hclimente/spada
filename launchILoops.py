#!/soft/devel/python-2.7/bin/python

import sys, os
from sh import *
from time import sleep

if(len(sys.argv) != 1):
	print("No arguments passed.")
	exit()

os.chdir(sys.argv[1])

pidQueue = []

for transcript in os.listdir("input"): #/ | egrep ENST`
	for configFile in os.listdir("input/" + transcript): #| egrep *net`
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

for job in pidQueue:
	while True:
		if cmdOut("ps --pid", job, " | grep -v", job, "| wc -l") == "0":
			break
		else:
			print("Awaiting for completion of iLoops jobs.")
			sleep(900)
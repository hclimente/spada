#!/soft/devel/python-2.7/bin/python

from libsmartas import cmdOut, waitPID

pidQueue = []

with open("Data/Input/TCGA_tags.txt", "r") as tags:
	for rawTag in tags:
		tag = rawTag.strip()
		cmd("Pipeline/standarizeInput.py TCGA", tag)

		if len(pidQueue) >= 3:
			waitPID(pidQueue)

		pid = cmdOut("SmartAS.py -f Data/Input/TCGA/" + tag + "/config.cfg >" + tag + ".log &; echo $!")
		pidQueue.append(pid)
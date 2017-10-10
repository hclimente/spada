import subprocess
import logging

logger = logging.getLogger("utils")

def cmd(base, *args):
	command = base
	for arg in args:
		command += " " + str(arg)

	logger.debug(command)
	subprocess.call(command, shell=True)

def cmdOut(*args):
	command = [ str(x) for x in args ]

	strCommand = ''
	for arg in command:
		strCommand += " " + arg

	logger.debug(strCommand)
	return subprocess.Popen(command, stdout=subprocess.PIPE)

def readTable(path, sep="\t", header=True):
	"""Read a table in a file, and generate a list of strings per row.
	Skips rows starting with "#".

	sep (str): field separator.
	header (bool): presence of a header, to discard the first row.
	"""
	counter = 0
	with open(path) as FILE:
		for line in FILE:
			if line[0]=="#":
				continue
			elif header and counter is 0:
				counter = 1
				continue

			yield line.strip().split(sep)

def geneclusterLaunch(tag,base,*args):
	command = base
	for argument in args:
		command += " " + argument

	cmd("qsub","-b y","-V","-N {}".format(tag),"-q short-low","-cwd","-e /data/users/hector/e-{}.log".format(tag),
		"-o /data/users/hector/o-{}.log".format(tag),command)

def readGeneset(sSetFile):

	geneSets = {}
	geneSetFile = "data/Databases/{1}".format(sSetFile)

	for line in readTable(geneSetFile,header=False):
		geneSet = line[0]
		genes = line[2:]
		geneSets[geneSet] = genes

	return geneSets

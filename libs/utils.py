from libs import options

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

	with open(tag+".sh","w") as configFile:
		configFile.write('#!/bin/sh\n')
		configFile.write('#$ -q normal\n')
		configFile.write('#$ -cwd\n')
		configFile.write("#$ -e /data/users/hector/{0}.log\n".format(tag))
		configFile.write("#$ -o /data/users/hector/{0}.log\n".format(tag))

		configFile.write("source ~/.bashrc\n")
		configFile.write(command+"\n")

	cmd("qsub","-N",tag,tag+".sh")

def launchJobs(gnNetwork,task,q="short"):
	import drmaa

	s = drmaa.Session()
	s.initialize()

	natSpec = ""
	if q=="normal":
		natSpec += "-q normal -l 'qname=normal' "
	else:
		natSpec += "-q short-high -l 'qname=short-high' "
	natSpec += "-cwd "
	natSpec += "-V "
	natSpec += "-N {}.{}".format(options.Options().tag,options.Options().step,task)

	cfg = options.Options().getCommandLineParameters(onlyModels=False)
	cfg.extend(["--parallel-range","$SGE_TASK_ID"])

	jt = s.createJobTemplate()
		
	jt.remoteCommand = 'pipeline/smartas.py'
	jt.args = cfg
	jt.joinFiles=True

	s.runBulkJobs(jt,1,len(gnNetwork.nodes()),options.Options().step)
	s.deleteJobTemplate(jt)

def launchSingleJob(task,name=""):
	
	if not name:
		import string
		import random
		name = 'SmartAS_task'
		name += ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))

	with open("{0}/Input/{1}.sh".format(options.Options().qout,name[:-1]),'w') as OUT:
		OUT.write("#$ -q sbi -l 'qname=sbi'\n")
		OUT.write("#$ -cwd\n")
		OUT.write("#$ -V\n")
		OUT.write("#$ -N {0}\n".format(name[:-1]))
		OUT.write("#$ -e {0}logs/{1}.out.txt\n".format(options.Options().qout,name[:-1]))
		OUT.write("#$ -o {0}logs/{1}.err.txt\n".format(options.Options().qout,name[:-1]))

		OUT.write(task)

	cmd("qsub {0}/Input/{1}.sh".format(options.Options().qout,name[:-1]))

def readGeneset(sSetFile):

	geneSets = {}
	geneSetFile = "{0}data/Databases/{1}".format(options.Options().wd,sSetFile)

	for line in readTable(geneSetFile,header=False):
		geneSet = line[0]
		genes = line[2:]
		geneSets[geneSet] = genes

	return geneSets
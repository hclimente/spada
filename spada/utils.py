class SpadaError(Exception):
	def __init__(self, value, logger):
		self.logger = logger
		self.value = value
	def __str__(self):
		logger.exception(value)
		return repr(self.value)

def readTable(path, sep = "\t", header = True, skipCommented = True):
	"""Read a table in a file, and generate a list of strings per row.
	Skips rows starting with "#".

	sep (str): field separator.
	header (bool): presence of a header, to discard the first row.
	"""
	counter = 0
	with open(path) as FILE:
		for line in FILE:
			if line[0]=="#" and skipCommented:
				continue
			elif header and counter is 0:
				counter = 1
				continue

			yield line.strip().split(sep)

def readGTF(gtf):
	for line in readTable(gtf, header=False):
		# read common fields
		parsed = {}
		parsed["chromosome"] = line[0]
		parsed["source"] 	 = line[1]
		parsed["feature"] 	 = line[2]
		parsed["start"] 	 = line[3]
		parsed["end"] 		 = line[4]
		parsed["score"] 	 = line[5]
		parsed["strand"] 	 = line[6]
		parsed["phase"] 	 = line[7]

		# read additional fields
		for field in line[8].strip().split(";"):
			if len(field):
				field = field.strip().split(" ")
				parsed[field[0]] = field[1].strip('"')

		yield parsed


def readPSIMITAB(psimitab):

	def parseField(field):
		parsed = {}

		if (len(field) > 1):
			x = field.split(":")
			parsed["type"] = x[0]

			if "(" in x[1]:
				x = x[1].strip(")").split("(")
				parsed["id"] = x[0]
				parsed["extra"] = x[1]
			else:
				parsed["id"] = x[1]

		return(parsed)

	for line in readTable(psimitab, skipCommented = False):

		parsed = {}
		parsed["geneA"] 		= [ parseField(x) for x in line[0].split("|") ]
		parsed["geneB"] 		= [ parseField(x) for x in line[1].split("|") ]
		parsed["symbolA"] 		= [ parseField(x) for x in line[2].split("|") ]
		parsed["symbolB"] 		= [ parseField(x) for x in line[3].split("|") ]
		parsed["aliasA"] 		= [ parseField(x) for x in line[4].split("|") ]
		parsed["aliasB"] 		= [ parseField(x) for x in line[5].split("|") ]
		parsed["method"] 		= [ parseField(x) for x in line[6].split("|") ]
		parsed["author"] 		= line[7]
		parsed["publication"] 	= [ parseField(x) for x in line[8].split("|") ]
		parsed["organismA"] 	= [ parseField(x) for x in line[9].split("|") ]
		parsed["organismB"] 	= [ parseField(x) for x in line[10].split("|") ]
		parsed["interaction"] 	= [ parseField(x) for x in line[11].split("|") ]
		parsed["source"] 		= [ parseField(x) for x in line[12].split("|") ]
		parsed["interactionId"] = [ parseField(x) for x in line[13].split("|") ]
		parsed["score"] 		= [ parseField(x) for x in line[14].split("|") ]

		yield parsed

def readFasta(fasta):
	with open(fasta) as FASTA:
		protein = ""
		sequence = ""

		for line in FASTA:
			if ">" in line:
				if sequence != "":
					yield protein, sequence
				protein  = line[1:].strip()
				sequence = ""

			else:
				sequence += line.strip()

		yield protein, sequence

def readGeneset(sSetFile):

	geneSets = {}
	geneSetFile = "data/Databases/{1}".format(sSetFile)

	for line in readTable(geneSetFile,header=False):
		geneSet = line[0]
		genes = line[2:]
		geneSets[geneSet] = genes

	return geneSets

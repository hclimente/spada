from spada.io import io

import logging

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
		parsed["tags"] 	 	 = set()

		# read additional fields
		for field in line[8].strip().split(";"):
			if len(field):
				k,v = field.strip().split(" ")
				v = v.strip('"')
				if k == 'tag':
					parsed['tags'].add(v)
				else:
					parsed[k] = v

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

def printSwitches(genes, txs, filename = "switches_spada.tsv"):

	logging.info("Writing switch information.")
	with open(filename, "w") as cList:
		cList.write("GeneId\tSymbol\tNormal_transcript\tTumor_transcript\t")
		cList.write("Normal_protein\tTumor_protein\tAnnotation\tDriverAnnotation\t")
		cList.write("NotNoise\tIsModel\tIsFunctional\tDriver\tSpecificDriver\tDruggable\t")
		cList.write("CDS_Normal\tCDS_Tumor\tCDS_change\tUTR_change\tPatients_affected\n")

		hallmarksDict = io.readGeneset("h.all.v5.0.entrez.gmt")
		bpDict = io.readGeneset("c5.bp.v4.0.entrez.gmt")

		for gene,info,switchDict,switch in genes.iterate_switches_byPatientNumber(txs,partialCreation=True,removeNoise=False):
			nIso = switch.nTranscript
			tIso = switch.tTranscript

			nUniprot 	= txs._net.node[switch.nTx]["Uniprot"]
			tUniprot 	= txs._net.node[switch.tTx]["Uniprot"]
			cdsChange 	= False
			utrChange 	= False

			if switch.cds_diff: 	 	cdsChange 	= True
			if switch.utr_diff: 		utrChange 	= True

			try:
				relevance = int(switch.is_functional)
			except Exception:
				relevance = None

			annotation,driverAnnotation = genes.getGeneAnnotation(gene,hallmarksDict,bpDict)

			cList.write("{}\t{}\t".format( gene, info["symbol"] ))
			cList.write("{}\t{}\t".format( nIso.name, tIso.name ))
			cList.write("{}\t{}\t".format( nUniprot, tUniprot ))
			cList.write("{}\t{}\t".format( annotation,driverAnnotation ))
			cList.write("{}\t{}\t".format( int(not switchDict["noise"]), int(switchDict["model"]) ))
			cList.write("{}\t{}\t".format( relevance, int(info["driver"]) ))
			cList.write("{}\t".format( int(info["specificDriver"]) ))
			cList.write("{}\t{}\t".format( int(info["druggable"]), int(bool(nIso.cds)) ))
			cList.write("{}\t{}\t".format( int(bool(tIso.cds)), int(cdsChange), ))
			cList.write("{}\t{}\n".format( int(utrChange), ",".join(switch.patients) ))

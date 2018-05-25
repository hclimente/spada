from spada.io import io

import logging
import numpy as np

class SpadaError(Exception):
	def __init__(self, value):
		self.logger = logging.getLogger('spada error')
		self.value = value
	def __str__(self):
		self.logger.exception(value)
		return repr(self.value)

def readTable(path, sep = "\t", header = True, skipCommented = True, keys = []):
	"""Read a table in a file, and generate a list of strings per row.
	Skips rows starting with "#".

	sep (str): field separator.
	header (bool): presence of a header, to discard the first row.
	"""
	headerRead = False
	with open(path) as FILE:
		for line in FILE:
			if line[0] == "#" and skipCommented:
				continue
			elif header and not headerRead:
				if not keys:
					keys = line.strip().split(sep)
				headerRead = True
				continue

			values = line.strip('\n').split(sep)

			if keys:
				if len(values) < len(keys):
					raise Exception("Parsing error in file {}, line {}.".format(path, line))
				yield dict(zip(keys, values))
			else:
				yield values

def readGTF(gtf):
	fields = ["chromosome","source","feature","start",
			  "end", "score", "strand", "phase", "tags"]
	for line in readTable(gtf, header = False, keys = fields):

		# read additional fields
		tags = set()
		for field in line['tags'].strip().split(";"):
			if len(field):
				k,v = field.strip().split(" ")
				v = v.strip('"')
				if k == 'tag':
					tags.add(v)
				else:
					line[k] = v
		line['tags'] = tags

		yield line


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

	fields = ["geneA","geneB","symbolA","symbolB","aliasA","aliasB","method",
			  "author","publication","organismA","organismB","interaction",
			  "source","interactionId","score"]

	for line in readTable(psimitab, skipCommented = False, keys = fields):

		parsed = {}

		for field in fields:
			if field == 'author':
				parsed[field] = line[field]
			else:
				parsed[field] = [ parseField(x) for x in line[field].split("|") ]

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

	with open(filename, "w") as OUT:
		OUT.write("Experiment\tGeneId\tSymbol\tNormal_transcript\t")
		OUT.write("Tumor_transcript\tDriver\tIs_main\tNot_noise\tIs_functional\t")
		OUT.write("CDS_normal\tCDS_tumor\tCDS_change\tUTR_change\tSamples\n")

		for gene,info,thisSwitch in genes.iterate_switches_byPatientNumber(txs, removeNoise = False):

			cdsChange = bool(thisSwitch.cds_diff)
			utrChange = bool(thisSwitch.utr_diff)

			driver = genes.isDriver(gene)

			OUT.write("{}\t".format( genes.tumor ))
			OUT.write("{}\t{}\t".format( gene, info["symbol"] ))
			OUT.write("{}\t{}\t".format( thisSwitch.nTx, thisSwitch.tTx ))
			OUT.write("%s\t%i\t" % ( driver, not thisSwitch.isNoise ))
			OUT.write("%i\t%i\t" % ( thisSwitch.isMain, thisSwitch.isFunctional))
			OUT.write("%i\t%i\t" % ( bool(thisSwitch.nTranscript.cds), bool(thisSwitch.tTranscript.cds) ))
			OUT.write("%i\t%i\t" % ( cdsChange, utrChange))
			OUT.write("{}\n".format( ",".join(thisSwitch.samples) ))

def parseExpression(FILE, header = False):

	for line in FILE:

		if not header:
			header = True
			continue

		xpr = line.strip().split('\t')
		tx = xpr.pop(0)
		xpr = np.array([xpr])
		xpr = xpr.astype(np.float)
		xpr = np.exp2(xpr) - .001

		yield tx,xpr

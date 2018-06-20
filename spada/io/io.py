from spada.io.error import SpadaError
from spada.biological_entities.gene_expression import GeneExpression

import logging
import numpy as np

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
		OUT.write("Experiment\tGeneId\tSymbol\tControl_transcript\t")
		OUT.write("Case_transcript\t")
		OUT.write("CDS_control\tCDS_case\tCDS_change\tUTR_change\tSamples\n")

		for gene,info,thisSwitch in genes.switches(txs):

			cdsChange = bool(thisSwitch.cds_diff)
			utrChange = bool(thisSwitch.utr_diff)

			OUT.write("{}\t{}\t{}\t".format( genes._name, gene, info["symbol"] ))
			OUT.write("{}\t{}\t".format( thisSwitch.ctrl, thisSwitch.case ))
			OUT.write("%i\t%i\t" % ( bool(thisSwitch.nTranscript.cds), bool(thisSwitch.tTranscript.cds) ))
			OUT.write("%i\t%i\t" % ( cdsChange, utrChange))
			OUT.write("{}\n".format( ",".join(thisSwitch.samples) ))

def parseExpression(ctrlFile, caseFile, genes, txs):

	gene2tx = getGene2Tx(txs)

	with open(ctrlFile, "r") as CTRL, open(caseFile, "r") as CASE:

		idsCtrl = readSamples(CTRL)
		idsCase = readSamples(CASE)

		# gene -> xpr
		expression = { }

		for (tx,ctrl),(tx2,case) in zip(parseExpressionLine(CTRL), parseExpressionLine(CASE)):

			if tx != tx2:
				raise SpadaError("Case and control expresion files mismatch: {} vs. {}.".format(tx, tx2))

			try:
				gene = txs._net.node[tx]["gene_id"]
			except KeyError:
				continue

			expression.setdefault(gene, GeneExpression(gene2tx[gene], idsCtrl, idsCase))
			expression[gene].addTx(tx, ctrl, case)

			if expression[gene].isComplete:
				yield gene, expression.pop(gene)

def parseExpressionLine(FILE, header = False):

	for line in FILE:

		if header:
			header = False
			continue

		xpr = line.strip().split('\t')
		tx = xpr.pop(0)

		xpr = np.array([xpr])
		xpr = xpr.astype(np.float)
		xpr = np.exp2(xpr) - .001
		xpr[xpr < 0] = 0

		yield tx,xpr

def readSamples(FILE):

	ids = FILE.readline().strip().split('\t')
	ids.pop(0)

	return ids

def getGene2Tx(txs):

	gene2tx = {}

	pairs = [ (t,i["gene_id"]) for t,i in txs.nodes(data=True) ]

	for tx,gene in pairs:
		gene2tx.setdefault(gene, set())
		gene2tx[gene].add(tx)

	return(gene2tx)

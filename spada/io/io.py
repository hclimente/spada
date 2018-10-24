from spada.io.error import SpadaError
from spada.bio.gene_expression import GeneExpression

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

def readGTF(gtf, tag_sep = ' '):
	fields = ["chromosome","source","feature","start",
			  "end", "score", "strand", "phase", "tags"]
	for line in readTable(gtf, header = False, keys = fields):

		line['start'] = int(line['start'])
		line['end'] = int(line['end'])

		# read additional fields
		tags = set()
		for field in line['tags'].strip().split(";"):
			if len(field):
				k,v = field.strip().split(tag_sep)
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
		OUT.write("Case_transcript\tCDS_control\tCDS_case\tCDS_change\t")
		OUT.write("5_UTR_change\t3_UTR_change\tSamples\n")

		for gene,info,thisSwitch in genes.switches(txs):

			cdsChange = bool(thisSwitch.cds_diff)
			utr5Change = bool(thisSwitch.utr_diff("5'"))
			utr3Change = bool(thisSwitch.utr_diff("3'"))

			OUT.write("{}\t{}\t{}\t".format( genes._name, gene, info["symbol"] ))
			OUT.write("{}\t{}\t".format( thisSwitch.ctrl, thisSwitch.case ))
			OUT.write("%i\t%i\t" % ( bool(thisSwitch.ctrlTranscript.cds), bool(thisSwitch.caseTranscript.cds) ))
			OUT.write("%i\t%i\t" % ( cdsChange, utr5Change))
			OUT.write("%i\t%s\n" % ( utr3Change, ",".join(sorted(thisSwitch.samples)) ))

def printSwitchesToGff(genes, txs, filename = "switches_spada.gff"):

	printed_txs = set()

	with open(filename, "w") as GTF:

		GTF.write('##gff-version 3\n')

		for _,_, thisSwitch in genes.switches(txs):
			for tx,isoform in [(thisSwitch.ctrlTranscript, thisSwitch.ctrlIsoform),
							   (thisSwitch.caseTranscript, thisSwitch.caseIsoform)]:
				name = tx._name

				if name in printed_txs:
					continue
				
				printed_txs.add(name)

				chromosome = txs[name]['chr']
				geneId = txs[name]['gene_id']
				symbol = genes[geneId]['symbol']

				gtfLine(GTF, chromosome, 'mRNA', 
						tx._tx_coordinates[0], tx._tx_coordinates[1], 
						tx._strand, 'ID={};Alias={}'.format(name, symbol))
				tags = 'Parent={}'.format(name)

				for start,end in tx._exons:
					gtfLine(GTF, chromosome, 'exon', start, end, tx._strand, tags)

				utr5 = list(tx.utr5.keys())
				utr3 = list(tx.utr3.keys())

				if utr5:
					if tx._strand == '+':
						gtfLine(GTF, chromosome, 'five_prime_UTR', utr5[0], utr5[-1], tx._strand, tags)
					elif tx._strand == '-':
						gtfLine(GTF, chromosome, 'five_prime_UTR', utr5[-1], utr5[0], tx._strand, tags)


				if utr3:
					if tx._strand == '+':
						gtfLine(GTF, chromosome, 'three_prime_UTR', utr3[0], utr3[-1], tx._strand, tags)
					elif tx._strand == '-':
						gtfLine(GTF, chromosome, 'three_prime_UTR', utr3[-1], utr3[0], tx._strand, tags)

				if isoform:

					cds = list(tx.cds.keys())

					if tx._strand == '+':
						gtfLine(GTF, chromosome, 'CDS', cds[0], cds[-1], tx._strand, tags)
						gtfLine(GTF, chromosome, 'start_codon', cds[0], cds[2], tx._strand, 
								'{};ID=cds_{}'.format(tags, name))
					elif tx._strand == '-':
						gtfLine(GTF, chromosome, 'CDS', cds[-1], cds[0], tx._strand, tags)
						gtfLine(GTF, chromosome, 'start_codon', cds[2], cds[0], tx._strand, 
								'{};ID=cds_{}'.format(tags, name))

					featureTypes = {'Pfam': isoform._pfam, 
									'Prosite': isoform._prosite, 
									'IDR': isoform._idr }
					for db, features in featureTypes.items():
						for feature_id, ranges in features.items():
							for s,e in ranges:
								start = isoform.structure[s - 1].genomicPosition
								end = isoform.structure[e - 1].genomicPosition
								gtfLine(GTF, txs[name]['chr'], 'translated_nucleotide_match', min(start,end), max(start,end), 
										txs[name]['strand'], 'CDS_matches=cds_{};Alias={};Note={};Target={} {} {}'.format(name, feature_id, db, txs[name]['chr'], s, e) )

			# isoform specific

def gtfLine(GTF, chromosome, feature, start, end, strand, tags):

	if start & (start > end):
		start2 = end
		end = start
		start = start2

	line  = '{}\tspada\t'.format(chromosome)
	line += '{}\t{}\t'.format(feature, start)
	line += '{}\t.\t'.format(end)
	line += '{}\t.\t'.format(strand)
	line += '{}\n'.format(tags)

	GTF.write(line)

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
				gene = txs[tx]["gene_id"]
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

def getGene2Tx(txs, with_aberrant = False):

	gene2tx = {}

	pairs = [ (t,i["gene_id"]) for t,i in txs.transcripts() if (with_aberrant or i['canonical']) ]

	for tx,gene in pairs:
		gene2tx.setdefault(gene, set())
		gene2tx[gene].add(tx)

	return(gene2tx)

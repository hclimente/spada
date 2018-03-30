from spada.methods import method

class StructuralAnalysis(method.Method):
	def __init__(self, gn_network, tx_network, isRandom = False):

		method.Method.__init__(self, __name__,gn_network,tx_network)

		self._tag = "_random" if isRandom else ""

	def run(self):

		self.featureAnalysis()
		self.ppiAnalysis()

		self.proteomeStatistics()

	def featureAnalysis(self):

		self.logger.info("Feature analysis.")

		with open("pfam_analysis{}.tsv".format(self._tag), "w") as PFAM, \
			 open("prosite_analysis{}.tsv".format(self._tag), "w") as PROSITE, \
			 open("idr_analysis{}.tsv".format(self._tag), "w") as IDR:

			self.writeDomainsHeader(PFAM)
			self.writeDomainsHeader(PROSITE)
			self.writeIDRHeader(IDR)

			for gene,info,thisSwitch in self._genes.iterate_switches_byPatientNumber(self._txs, removeNoise = False):
				pfam_change		= thisSwitch.analyzeDomains("Pfam")
				prosite_change	= thisSwitch.analyzeDomains("Prosite")
				idr_change   	= thisSwitch.analyzeIDR(0.2)

				self.writeDomains(PFAM, 'Pfam', gene, thisSwitch, pfam_change)
				self.writeDomains(PROSITE, 'Prosite', gene, thisSwitch, prosite_change)
				self.writeIDR(IDR, gene, thisSwitch, idr_change)

	def ppiAnalysis(self):

		self.logger.info("PPI analysis.")

		with open("ppi_analysis{}.tsv".format(self._tag), "w") as PPI:

			self.writePPIHeader(PPI)

			for gene,info,thisSwitch in self._genes.iterate_switches_byPatientNumber(self._txs, removeNoise = False):

				DDIchanges = self.analyzeDDIs(thisSwitch)
				self.writePPI(PPI, gene, thisSwitch, DDIchanges)

	def proteomeStatistics(self):

		with open("proteome_features{}.tsv".format(self._tag), "w") as OUT:

			self.writeProteomeHeader(OUT)

			genes = [ (x,i['expressedTxsN']) for x,i in self._genes.nodes(data=True) ]

			for gene, txs in genes:

				if not txs:
					continue

				expression = [ (tx,self._txs._net.node[tx]['median_TPM_N']) for tx in txs ]
				expression = sorted(expression, key = lambda t: -t[1])
				tx = expression[0][0]
				txInfo = self._txs._net.node[tx]

				for featureType in ['Pfam','Prosite']:
					for feature in txInfo[featureType]:
						i = 1
						for start,end in txInfo[featureType][feature]:
							OUT.write("{}\t{}\t".format(self._genes.tumor, txInfo['gene_id']))
							OUT.write("{}\t{}\t".format(tx, txInfo['median_TPM_N']))
							OUT.write("{}\t{}\t".format(featureType, feature))
							OUT.write("{}\t{}\t".format(i, end - start))
							OUT.write("{}\t{}\n".format(start, end))

						i += 1

	def analyzeDDIs(self, thisSwitch):

		DDIchanges = {}
		PPIs = { thisSwitch.nTx: {}, thisSwitch.tTx: {} }

		for tx in [thisSwitch.nTx, thisSwitch.tTx]:
			ppis = self._txs._net.edges(tx, data=True)

			for tx1, tx2, eDict in ppis:
				pairs = [ (next(iter(x)), next(iter(x))) for x in eDict['ddi'] if len(x) == 1 ] # case for homodomain interactions
				pairs.extend([ tuple(p) for p in eDict['ddi'] if len(p) > 1 ])					# case for heterodomain interactions
				domainsTx1 = self._txs._net.node[tx1]['Pfam'].keys()
				domainsTx2 = self._txs._net.node[tx2]['Pfam'].keys()

				PPIs[tx1][tx2] = set()
				for dA, dB in pairs:
					if dA in domainsTx1 and dB in domainsTx2:
						PPIs[tx1][tx2].add(dA + "@" + dB)
					if dA in domainsTx2 and dB in domainsTx1:
						PPIs[tx1][tx2].add(dB + "@" + dA)

		partners = set(PPIs[thisSwitch.nTx].keys()) | set(PPIs[thisSwitch.tTx].keys())

		for p in partners:

			nDDI = PPIs[thisSwitch.nTx].get(p, set())
			tDDI = PPIs[thisSwitch.tTx].get(p, set())
			DDIchanges[p] = { "nDDIs": nDDI - tDDI,
							  "tDDIs": tDDI - nDDI,
							  "bothDDIs": nDDI & tDDI,
							  "what": "Unaffected" }

			if DDIchanges[p]["bothDDIs"]:
				if DDIchanges[p]["nDDIs"] or DDIchanges[p]["tDDIs"]:
					DDIchanges[p]["what"] = "Affected"
			elif DDIchanges[p]["nDDIs"] and not DDIchanges[p]["tDDIs"]:
				DDIchanges[p]["what"] = "Lost_in_tumor"
			elif not DDIchanges[p]["nDDIs"] and DDIchanges[p]["tDDIs"]:
				DDIchanges[p]["what"] = "Gained_in_tumor"

		return DDIchanges

	def writeDomainsHeader(self, OUT):
		OUT.write("Experiment\tGeneId\tNormal_transcript\tTumor_transcript\t")
		OUT.write("Feature_type\tFeature\tWhat\tIndex\tNormal_start\tNormal_end\t")
		OUT.write("Tumor_start\tTumor_end\tNormal_MacroScore\tNormal_MicroScore\t")
		OUT.write("Normal_Jaccard\tTumor_MacroScore\tTumor_MicroScore\tTumor_Jaccard\n")

	def writeDomains(self, OUT, featureType, gene, thisSwitch, changes):
		for c in changes:
			OUT.write("{}\t{}\t{}\t".format(self._genes.tumor, gene, thisSwitch.nTx))
			OUT.write("{}\t{}\t{}\t".format(thisSwitch.tTx, featureType, c["feature"]))
			OUT.write("{}\t{}\t{}\t".format(c["what"], c["index"], c['nStart']))
			OUT.write("{}\t{}\t{}\t".format(c["nEnd"], c["tStart"], c['tEnd']))
			OUT.write("{}\t{}\t{}\t".format(c["nM"], c["nm"], c["nJ"]))
			OUT.write("{}\t{}\t{}\n".format(c["tM"], c["tm"], c["tJ"]))

	def writeIDRHeader(self, OUT):
		OUT.write("Experiment\tGeneId\tNormal_transcript\tTumor_transcript\t")
		OUT.write("What\tSequence\tStartPos\tEndPos\t")
		OUT.write("microScore\tmacroScore\tJaccard\n")

	def writeIDR(self, OUT, gene, thisSwitch, changes):
		for c in changes:
			OUT.write("{}\t".format( self._genes.tumor ))
			OUT.write("{}\t{}\t{}\t".format(gene, thisSwitch.nTx, thisSwitch.tTx))
			OUT.write("{}\t{}\t{}\t".format(c["what"], c["feature"], c["start"]))
			OUT.write("{}\t{}\t{}\t{}\n".format(c["end"], c["M"], c["m"], c["J"]))

	def writePPIHeader(self, OUT):
		OUT.write("Experiment\tGeneId\tNormal_transcript\tTumor_transcript\t")
		OUT.write("Other_gene\tOther_symbol\tOther_transcript\tWhat\t")
		OUT.write("#nDDIs\t#tDDIs\t#BothDDIs\tnDDIs\ttDDIs\tBothDDIs\n")

	def writePPI(self, OUT, gene, thisSwitch, DDIchanges):
		for tx, ddis in DDIchanges.items():

			other_gene = self._txs._net.node[tx]["gene_id"]

			OUT.write("{}\t{}\t".format( self._genes.tumor, gene ))
			OUT.write("{}\t{}\t".format( thisSwitch.nTx, thisSwitch.tTx ))
			OUT.write("{}\t{}\t".format( other_gene, self._genes._net.node[other_gene]["symbol"] ))
			OUT.write("{}\t{}\t".format( tx, ddis["what"] ))
			OUT.write("{}\t{}\t".format( len(ddis["nDDIs"]), len(ddis["tDDIs"]) ))
			OUT.write("{}\t{}\t".format( len(ddis["bothDDIs"]), ";".join(ddis["nDDIs"]) ))
			OUT.write("{}\t{}\n".format( ";".join(ddis["tDDIs"]), ";".join(ddis["bothDDIs"]) ))

	def writeProteomeHeader(self, OUT):
		OUT.write("Experiment\tGeneId\tTranscript\tExpression\t")
		OUT.write("Feature_type\tFeature\tIndex\tLength\tStart\tEnd\n")

if __name__ == '__main__':
	pass

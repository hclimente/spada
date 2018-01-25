from spada.methods import method

class StructuralAnalysis(method.Method):
	def __init__(self, gn_network, tx_network, isRandom = False):

		method.Method.__init__(self, __name__,gn_network,tx_network)

		self._tag = "_random" if isRandom else ""

	def run(self):

		self.logger.info("Feature analysis.")
		self.featureAnalysis()

		self.logger.info("PPI analysis.")
		self.ppiAnalysis()

	def featureAnalysis(self):

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

				self.writeDomains(PFAM, gene, thisSwitch, pfam_change)
				self.writeDomains(PROSITE, gene, thisSwitch, prosite_change)
				self.writeIDR(IDR, gene, thisSwitch, idr_change)

	def ppiAnalysis(self):

		with open("ppi_analysis{}.tsv".format(self._tag), "w") as PPI:

			self.writePPIHeader(PPI)

			for gene,info,thisSwitch in self._genes.iterate_switches_byPatientNumber(self._txs, removeNoise = False):

				DDIchanges = self.analyzeDDIs(thisSwitch)
				self.writePPI(PPI, gene, thisSwitch, DDIchanges)

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
		OUT.write("Gene\tNormalTranscript\tTumorTranscript\t")
		OUT.write("What\tFeature\tIndex\tnMacroScore\t")
		OUT.write("nMicroScore\tnJaccard\ttMacroScore\ttMicroScore\ttJaccard\n")

	def writeDomains(self, OUT, gene, thisSwitch, changes):
		for c in changes:
			OUT.write("{}\t{}\t{}\t".format(gene, thisSwitch.nTx, thisSwitch.tTx))
			OUT.write("{}\t{}\t{}\t".format(c["what"], c["feature"], c["index"]))
			OUT.write("{}\t{}\t{}\t".format(c["nM"], c["nm"], c["nJ"]))
			OUT.write("{}\t{}\t{}\n".format(c["tM"], c["tm"], c["tJ"]))

	def writeIDRHeader(self, OUT):
		OUT.write("Gene\tNormalTranscript\tTumorTranscript\t")
		OUT.write("What\tSequence\tStartPos\tEndPos\t")
		OUT.write("microScore\tmacroScore\tJaccard\n")

	def writeIDR(self, OUT, gene, thisSwitch, changes):
		for c in changes:
			OUT.write("{}\t{}\t{}\t".format(gene, thisSwitch.nTx, thisSwitch.tTx))
			OUT.write("{}\t{}\t{}\t".format(c["what"], c["feature"], c["start"]))
			OUT.write("{}\t{}\t{}\t{}\n".format(c["end"], c["M"], c["m"], c["J"]))

	def writePPIHeader(self, OUT):
		OUT.write("Switched_gene\tNormalTranscript\tTumorTranscript\t")
		OUT.write("Other_gene\tOther_transcript\tWhat\t")
		OUT.write("#nDDIs\t#tDDIs\t#BothDDIs\n")
		OUT.write("nDDIs\ttDDIs\tBothDDIs\n")

	def writePPI(self, OUT, gene, thisSwitch, DDIchanges):
		for tx, ddis in DDIchanges.items():
			OUT.write("{}\t{}\t{}\t".format(gene, thisSwitch.nTx, thisSwitch.tTx))
			OUT.write("{}\t{}\t".format(self._txs._net.node[tx]["gene_id"], tx))
			OUT.write("{}\t".format(ddis["what"]))
			OUT.write("{}\t{}\t".format(len(ddis["nDDIs"]), len(ddis["tDDIs"])))
			OUT.write("{}\t".format(len(ddis["bothDDIs"])))
			OUT.write("{}\t{}\t".format(";".join(ddis["nDDIs"]), ";".join(ddis["tDDIs"]) ))
			OUT.write("{}\n".format(";".join(ddis["bothDDIs"]) ))

if __name__ == '__main__':
	pass

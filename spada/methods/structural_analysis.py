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

				ppi_change = thisSwitch.analyzePPIs()

				self.writePPI(PPI, gene, thisSwitch, ppi_change)

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
		OUT.write("Other_gene\tOther_transcript\tWhat\tInvolved_domains\n")

	def writePPI(self, PPI, gene, thisSwitch, ppi_change):
		pass

if __name__ == '__main__':
	pass

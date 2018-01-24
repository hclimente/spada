from spada import utils
from spada.methods import method

import os

class StructuralAnalysis(method.Method):
	def __init__(self, gn_network, tx_network, isRandom = False):

		if not os.path.exists("structural_analysis"):
			utils.cmd("mkdir","structural_analysis")

		method.Method.__init__(self, __name__,gn_network,tx_network)

		tag = "_random" if isRandom else ""

		self._pfam_file = "structural_analysis/pfam_analysis{}.tsv".format(tag)
		self.PFAM = open(self._pfam_file,"w")
		self.writeDomainsHeader(self.PFAM)

		self._idr_file = "structural_analysis/idr_analysis{}.tsv".format(tag)
		self.IDR = open(self._idr_file,"w")
		self.writeIDRHeader(self.IDR)

		self._prosite_file = "structural_analysis/prosite_analysis{}.tsv".format(tag)
		self.PROSITE = open(self._prosite_file,"w")
		self.writeDomainsHeader(self.PROSITE)

	def run(self):
		self.logger.info("Structural analysis.")

		for gene,info,thisSwitch in self._genes.iterate_switches_byPatientNumber(self._txs, removeNoise = False):
			#if gene == "ENSG00.5": import pdb; pdb.set_trace()
			pfam_change		= thisSwitch.analyzeDomains("Pfam")
			prosite_change	= thisSwitch.analyzeDomains("Prosite")
			idr_change   	= thisSwitch.analyzeIDR(0.2)

			self.writeDomains(self.PFAM, gene, thisSwitch, pfam_change)
			self.writeDomains(self.PROSITE, gene, thisSwitch, prosite_change)
			self.writeIDR(self.IDR, gene, thisSwitch, idr_change)

		self.PFAM.close()
		self.IDR.close()
		self.PROSITE.close()

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

if __name__ == '__main__':
	pass

from spada.methods import method

import itertools
import operator
import random

class SimulateSwitches(method.Method):
	def __init__(self, gn_network, tx_network):

		method.Method.__init__(self, __name__, gn_network, tx_network)
		self.MAX_SWITCHES = 5

	def run(self):
		self.logger.info("Generating random switches.")

		# random: two random transcripts among those expressed in tumor or normal
		self.sampleTranscripts_random()

		# random: most expressed isoform is normal
		self.sampleTranscripts_fixNormal()

	def sampleTranscripts_fixNormal(self):

		switches = []

		for gene,info in self._genes.iterate_genes_byPatientNumber():

			if len(set(info["expressedTxsN"]) & set(info["expressedTxsT"])) < 2:
				next

			# set normal isoform as the most expressed in normal, shuffle the rest for tumor
			txExpression = [ (x,self._txs._net.node[x]["median_TPM_N"]) for x in info["expressedTxsN"] ]
			nTx = max(txExpression, key=operator.itemgetter(1))[0]
			nTxExpression = max(txExpression, key=operator.itemgetter(1))[1]

			if nTx not in info["expressedTxsN"]:
				self.logger.warning("Median most expressed transcript from gene {} is not considered expressed. \
					Probably due to the threshold applied. TPM={}. Will be skipped. ".format(gene,txExpression))
				continue

			txs = list(info["expressedTxsT"])
			if nTx in txs:
				txs.remove(nTx)

			numSwitches = self.MAX_SWITCHES
			if len(txs) < self.MAX_SWITCHES:
				numSwitches = len(txs)

			random.shuffle(txs)

			for i in range(numSwitches):
				tTx = txs[i]
				switches.append((gene, nTx, tTx))

		return(switches)

	def sampleTranscripts_random(self):

		switches = []

		for gene,info in self._genes.iterate_genes_byPatientNumber():

			allExpressedTxs = set(info["expressedTxsN"]) & set(info["expressedTxsT"])
			if len(allExpressedTxs) < 2:
				next

			allSwitches = [ x for x in itertools.combinations(allExpressedTxs,2) ]
			random.shuffle(allSwitches)

			allSwitches = allSwitches[0:self.MAX_SWITCHES]

			for oneSwitch in allSwitches:
				nTx = oneSwitch[0]
				tTx = oneSwitch[1]
				switches.append((gene, nTx, tTx))

		return(switches)

if __name__ == '__main__':
	pass

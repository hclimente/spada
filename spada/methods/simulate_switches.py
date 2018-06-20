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

		# random: two random transcripts among those expressed in case or control
		self.sampleTranscripts_random()

		# random: most expressed isoform is control
		self.sampleTranscripts_fixControl()

	def sampleTranscripts_fixControl(self):

		switches = []

		for gene,info in self._genes.genes():

			if len(set(info["expressedTxsN"]) & set(info["expressedTxsT"])) < 2:
				next

			# set control isoform as the most expressed in control, shuffle the rest for case
			txExpression = [ (x,self._txs._net.node[x]["median_TPM_N"]) for x in info["expressedTxsN"] ]
			ctrl = max(txExpression, key=operator.itemgetter(1))[0]
			ctrlExpression = max(txExpression, key=operator.itemgetter(1))[1]

			if ctrl not in info["expressedTxsN"]:
				self.logger.warning("Median most expressed transcript from gene {} is not considered expressed. \
					Probably due to the threshold applied. TPM={}. Will be skipped. ".format(gene,txExpression))
				continue

			txs = list(info["expressedTxsT"])
			if ctrl in txs:
				txs.remove(ctrl)

			numSwitches = self.MAX_SWITCHES
			if len(txs) < self.MAX_SWITCHES:
				numSwitches = len(txs)

			random.shuffle(txs)

			for i in range(numSwitches):
				case = txs[i]
				switches.append((gene, ctrl, case))

		return(switches)

	def sampleTranscripts_random(self):

		switches = []

		for gene,info in self._genes.genes():

			allExpressedTxs = set(info["expressedTxsN"]) & set(info["expressedTxsT"])
			if len(allExpressedTxs) < 2:
				next

			allSwitches = [ x for x in itertools.combinations(allExpressedTxs,2) ]
			random.shuffle(allSwitches)

			allSwitches = allSwitches[0:self.MAX_SWITCHES]

			for oneSwitch in allSwitches:
				ctrl = oneSwitch[0]
				case = oneSwitch[1]
				switches.append((gene, ctrl, case))

		return(switches)

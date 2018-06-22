from spada.biological_entities.switch import LiteSwitch
from spada.io import io
from spada.methods import method

import itertools
import numpy as np
import operator
import random

class SimulateSwitches(method.Method):
	def __init__(self, annotation = 'annotation.pklz'):

		method.Method.__init__(self, __name__, annotation)
		self.max = 5

	def run(self, ctrlFile, caseFile, method, threshold):

		self.logger.info("Generating random switches.")

		for gene,expression in io.parseExpression(ctrlFile, caseFile, self._genes, self._txs):

			expressed = np.median(expression._expressionCtrl, axis = 1) > threshold
			txs = [ tx for tx,e in zip(expression._storedTxs, expressed) if e ]

			if method == 'random':
				# random: two random transcripts among those expressed in case or control
				switches = self.sampleTranscripts_random(txs)
			elif method == 'fix_expressed':
				# random: most expressed isoform is control
				tx,tpm = expression._top_ctrl
				switches = self.sampleTranscripts_fixControl(tx, txs.remove(tx))
			elif method == 'fix_main':
				# random: most expressed isoform is control
				tx = [ t for t in txs if self._txs._net.node[t]["main"] ]
				if tx:
					tx = tx[0]
				else:
					continue
				switches = self.sampleTranscripts_fixControl(tx, txs.remove(tx))

			for ctrl,case in switches:
				if self._genes.valid_switch(gene, ctrl, case, self._txs):
					thisSwitch = LiteSwitch(ctrl, case, [])
					self._genes.update_node("switches", thisSwitch, full_name = gene)

		io.printSwitches(self._genes, self._txs, "switches_simulated_{}.tsv".format(method))

	def sampleTranscripts_fixControl(self, ctrl, cases):

		if not cases:
			return []

		switches = []
		random.shuffle(cases)

		for i in range(min(len(cases), self.max)):
			switches.append((ctrl, cases[i]))

		return switches

	def sampleTranscripts_random(self, txs):

		switches = [ x for x in itertools.combinations(txs, 2) ]
		random.shuffle(switches)

		switches = switches[0:self.max]

		return switches

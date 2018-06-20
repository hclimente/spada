from spada.biological_entities.switch import LiteSwitch
from spada.io import io
from spada.methods import method

class ComputeSwitches(method.Method):

	def __init__(self, gn_network, tx_network):

		method.Method.__init__(self, __name__,gn_network,tx_network)

		self._genes.flushSwitches()

	def run(self, ctrlFile, caseFile, minExpression):

		self.findSwitches(ctrlFile, caseFile, minExpression)
		io.printSwitches(self._genes, self._txs)
		self._genes.saveNetwork("genes.pkl")

	def findSwitches(self, ctrlFile, caseFile, minExpression):

		for gene, geneExpression in io.parseExpression(ctrlFile, caseFile, self._genes, self._txs):

			if len(geneExpression._allTxs) > 1:

				switches = geneExpression.detectSwitches(minExpression)

				for (ctrl,case),samples in switches.items():
					thisSwitch = LiteSwitch(ctrl, case, samples)
					self._genes.update_node("switches", thisSwitch, full_name = gene)

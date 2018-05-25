from spada.methods import method

class GetSwitches(method.Method):
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self, __name__, gn_network, tx_network)
		self._genes.flushSwitches()

	def run(self, switchesFile):

		self._genes.readSwitches(switchesFile, self._txs)
		self._genes.calculateCompatibilityTable()
		self._genes.saveNetwork("genes.pkl")

if __name__ == '__main__':
	pass

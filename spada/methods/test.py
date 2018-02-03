from spada.methods import method

class Test(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork=False):
		method.Method.__init__(self,__name__,gn_network,tx_network)

	def run(self):

		for u,v in self._gene_network._net.edges():
			self._gene_network._net.remove_edge(u, v)


		self._gene_network.importEduardInteractions()
		self._gene_network.saveNetwork("geneNetwork.updatedIntx.pkl")

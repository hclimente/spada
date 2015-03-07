from libs import options
from methods import method

import cPickle
import pdb

class Test(method.Method):
	def __init__(self, gn_network, tx_network,gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		for x in self._gene_network.nodes():
			self._gene_network._net.node[x]["DriverType"] = None
			self._gene_network._net.node[x]["ASDriver"] = False
		self._gene_network.readGeneInfo()
		self._gene_network.saveNetwork("geneNetwork.pkl")
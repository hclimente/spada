from libs import options
from methods import method

import cPickle
import pdb

class Test(method.Method):
	def __init__(self, gn_network, tx_network,gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		self._gene_network.calculateCompatibilityTable()
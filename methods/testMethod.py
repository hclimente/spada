from libs import options
from libs import utils
from methods import method

class Test(method.Method):
	def __init__(self,gn_network,tx_network,gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		from interface import out_network

		out_network.outCandidateList(self._gene_network,self._transcript_network)
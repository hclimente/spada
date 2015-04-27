from libs import options
from libs import utils
from methods import method

class Test(method.Method):
	def __init__(self,gn_network,tx_network,gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,partialCreation=False,removeNoise=True,only_models=True):
			thisSwitch.analyzeSplicing()
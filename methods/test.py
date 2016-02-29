from libs import options
from libs import utils
from methods import method

class Test(method.Method):
	def __init__(self,gn_network,tx_network,gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		from network import ucsc_isoform_network
		self._transcript_network = ucsc_isoform_network.UCSCIsoformNetwork()
		self.logger.info("Creating transcript network.")
		self._transcript_network.importTranscriptome()

		self._transcript_network.saveNetwork("txNetwork.pkl")
#!/soft/devel/python-2.7/bin/python

from libs.utils import *
from network import network
from methods import method
from libs import options

class NetworkAnalysis(method.Method):
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self, __name__, gn_network, tx_network)

	def run(self):
		self.logger.info("GUILD analysis.")

if __name__ == '__main__':
	pass
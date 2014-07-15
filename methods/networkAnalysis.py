#!/soft/devel/python-2.7/bin/python

from libs.utils import *
from network import network
from methods import method
from libs import options

class NetworkAnalysis(method.Method):
	def __init__(self, gn_network, tx_network):
		self._gene_network			= gn_network
		self._transcript_network	= tx_network
		
		method.Method.__init__(self)

	def run(self):
		logging.info("GUILD analysis.")

if __name__ == '__main__':
	pass
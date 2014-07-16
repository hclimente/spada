#!/soft/devel/python-2.7/bin/python

import logging

class Method:
	def __init__(self, name, gn_network, tx_network):
		self.logger = logging.getLogger(name)
		self._gene_network			= gn_network
		self._transcript_network	= tx_network